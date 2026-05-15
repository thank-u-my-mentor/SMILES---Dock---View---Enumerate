#!/usr/bin/env python
"""Mine ChEMBL for literature molecules similar to docked or expert SMILES."""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, QED, rdMolDescriptors
from rdkit.ML.Cluster import Butina

SMILES_SKILL = Path(__file__).resolve().parents[2] / "smiles-to-vina-docking" / "scripts"
sys.path.insert(0, str(SMILES_SKILL))
import dock_utils as du  # noqa: E402


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def clean_float(value: object, default: float = 0.0) -> float:
    try:
        if value in ("", None):
            return default
        return float(value)
    except Exception:
        return default


def canonical(smiles: str | None) -> str:
    return du.canonical(smiles) or ""


def mol(smiles: str | None) -> Chem.Mol | None:
    return du.mol_from_smiles(smiles)


def fp(smiles: str) -> object | None:
    m = mol(smiles)
    if m is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)


def properties(smiles: str) -> dict[str, object]:
    m = mol(smiles)
    if m is None:
        return {}
    return {
        "qed": round(float(QED.qed(m)), 4),
        "mw": round(float(Descriptors.MolWt(m)), 3),
        "logp": round(float(Crippen.MolLogP(m)), 3),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(m)), 3),
        "heavy_atoms": int(m.GetNumHeavyAtoms()),
    }


def contact_score(row: dict[str, str]) -> float:
    hydrophobic = clean_float(row.get("hydrophobic_count") or row.get("hydrophobic_count ç–æ°´"))
    return (
        3.0 * clean_float(row.get("hbond_count"))
        + 0.7 * hydrophobic
        + 0.05 * clean_float(row.get("vdw_contact_count"))
        + 2.0 * clean_float(row.get("pi_contact_count"))
        + 2.0 * clean_float(row.get("ch_pi_count"))
    )


def load_query_rows(args: argparse.Namespace) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if args.history_dir:
        history = du.read_csv(args.history_dir.expanduser() / du.LEDGER_FILENAME)
        for row in history:
            smiles = canonical(row.get("canonical_smiles") or row.get("smiles") or row.get("input_smiles"))
            if not smiles:
                continue
            props = properties(smiles)
            heavy = max(1.0, float(props.get("heavy_atoms", 1)))
            rows.append(
                {
                    "query_id": row.get("seq_id", ""),
                    "query_nickname": row.get("nickname", ""),
                    "query_smiles": smiles,
                    "query_affinity": row.get("affinity_kcal_mol", ""),
                    "query_contact_score": round(contact_score(row), 4),
                    "query_contact_efficiency": round(contact_score(row) / heavy, 4),
                    **{f"query_{key}": value for key, value in props.items()},
                }
            )
    for i, smiles in enumerate(args.expert_smiles or [], start=1):
        can = canonical(smiles)
        if can:
            props = properties(can)
            rows.append(
                {
                    "query_id": f"EXPERT_{i:02d}",
                    "query_nickname": "expert",
                    "query_smiles": can,
                    "query_affinity": "",
                    "query_contact_score": "",
                    "query_contact_efficiency": "",
                    **{f"query_{key}": value for key, value in props.items()},
                }
            )
    dedup: dict[str, dict[str, object]] = {}
    for row in rows:
        dedup.setdefault(str(row["query_smiles"]), row)
    return list(dedup.values())


def cluster_queries(rows: list[dict[str, object]], cutoff: float) -> list[dict[str, object]]:
    fps = [fp(str(row["query_smiles"])) for row in rows]
    valid = [(i, f) for i, f in enumerate(fps) if f is not None]
    if not valid:
        return []
    dists: list[float] = []
    valid_fps = [item[1] for item in valid]
    for i in range(1, len(valid_fps)):
        sims = DataStructs.BulkTanimotoSimilarity(valid_fps[i], valid_fps[:i])
        dists.extend([1.0 - sim for sim in sims])
    clusters = Butina.ClusterData(dists, len(valid_fps), 1.0 - cutoff, isDistData=True)
    representatives: list[dict[str, object]] = []
    for cluster_id, cluster in enumerate(clusters, start=1):
        members = [rows[valid[index][0]] for index in cluster]
        members.sort(
            key=lambda row: (
                clean_float(row.get("query_contact_efficiency")),
                clean_float(row.get("query_qed")),
            ),
            reverse=True,
        )
        rep = dict(members[0])
        rep["cluster_id"] = f"C{cluster_id:04d}"
        rep["cluster_size"] = len(members)
        representatives.append(rep)
    representatives.sort(key=lambda row: int(row["cluster_size"]), reverse=True)
    return representatives


def chembl_similarity_hits(smiles: str, threshold: int, limit: int) -> list[dict[str, object]]:
    from chembl_webresource_client.new_client import new_client

    hits = new_client.similarity.filter(smiles=smiles, similarity=threshold)
    out: list[dict[str, object]] = []
    for hit in hits[:limit]:
        structures = hit.get("molecule_structures") or {}
        out.append(
            {
                "chembl_id": hit.get("molecule_chembl_id", ""),
                "chembl_similarity": hit.get("similarity", ""),
                "chembl_pref_name": hit.get("pref_name", ""),
                "chembl_smiles": structures.get("canonical_smiles", ""),
                "max_phase": hit.get("max_phase", ""),
            }
        )
    return out


def chembl_activity_rows(molecule_id: str, target_keyword: str, min_pchembl: float, limit: int) -> list[dict[str, object]]:
    from chembl_webresource_client.new_client import new_client

    acts = new_client.activity.filter(molecule_chembl_id=molecule_id, pchembl_value__gte=min_pchembl).only(
        "molecule_chembl_id",
        "target_chembl_id",
        "target_pref_name",
        "standard_type",
        "standard_relation",
        "standard_value",
        "standard_units",
        "pchembl_value",
        "ligand_efficiency",
        "document_chembl_id",
        "assay_chembl_id",
    )
    out: list[dict[str, object]] = []
    keyword = target_keyword.lower().strip()
    for act in acts:
        pref = str(act.get("target_pref_name", ""))
        if keyword and keyword not in pref.lower():
            continue
        le = act.get("ligand_efficiency") or {}
        out.append(
            {
                "activity_target_chembl_id": act.get("target_chembl_id", ""),
                "activity_target_pref_name": pref,
                "standard_type": act.get("standard_type", ""),
                "standard_relation": act.get("standard_relation", ""),
                "standard_value": act.get("standard_value", ""),
                "standard_units": act.get("standard_units", ""),
                "pchembl_value": act.get("pchembl_value", ""),
                "ligand_efficiency_LE": le.get("le", ""),
                "ligand_efficiency_LLE": le.get("lle", ""),
                "ligand_efficiency_BEI": le.get("bei", ""),
                "ligand_efficiency_SEI": le.get("sei", ""),
                "document_chembl_id": act.get("document_chembl_id", ""),
                "assay_chembl_id": act.get("assay_chembl_id", ""),
            }
        )
        if len(out) >= limit:
            break
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history-dir", type=Path)
    parser.add_argument("--expert-smiles", action="append", default=[])
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--cluster-threshold", type=float, default=0.72)
    parser.add_argument("--top-clusters", type=int, default=50)
    parser.add_argument("--similarity-threshold", type=int, default=70)
    parser.add_argument("--hits-per-query", type=int, default=10)
    parser.add_argument("--activities-per-hit", type=int, default=20)
    parser.add_argument("--target-keyword", default="kinase", help="Filter activity target name; use empty string to keep all targets")
    parser.add_argument("--min-pchembl", type=float, default=5.0)
    parser.add_argument("--sleep", type=float, default=0.05)
    args = parser.parse_args()

    try:
        import chembl_webresource_client  # noqa: F401
    except Exception as exc:
        raise SystemExit("Install first: pip install chembl_webresource_client") from exc

    outdir = args.outdir.expanduser().resolve()
    queries = load_query_rows(args)
    representatives = cluster_queries(queries, args.cluster_threshold)
    representatives = representatives[: args.top_clusters]
    rep_fields = [
        "cluster_id", "cluster_size", "query_id", "query_nickname", "query_smiles",
        "query_affinity", "query_contact_score", "query_contact_efficiency",
        "query_qed", "query_mw", "query_logp", "query_tpsa", "query_heavy_atoms",
    ]
    write_csv(outdir / "dock_history_structure_clusters.csv", representatives, rep_fields)

    hit_rows: list[dict[str, object]] = []
    activity_rows: list[dict[str, object]] = []
    for rep in representatives:
        hits = chembl_similarity_hits(str(rep["query_smiles"]), args.similarity_threshold, args.hits_per_query)
        for hit in hits:
            hit_row = {**rep, **hit}
            hit_rows.append(hit_row)
            activities = chembl_activity_rows(
                str(hit.get("chembl_id", "")),
                args.target_keyword,
                args.min_pchembl,
                args.activities_per_hit,
            )
            for activity in activities:
                activity_rows.append({**hit_row, **activity})
            time.sleep(args.sleep)
    hit_fields = rep_fields + ["chembl_id", "chembl_similarity", "chembl_pref_name", "chembl_smiles", "max_phase"]
    activity_fields = hit_fields + [
        "activity_target_chembl_id", "activity_target_pref_name", "standard_type",
        "standard_relation", "standard_value", "standard_units", "pchembl_value",
        "ligand_efficiency_LE", "ligand_efficiency_LLE", "ligand_efficiency_BEI",
        "ligand_efficiency_SEI", "document_chembl_id", "assay_chembl_id",
    ]
    write_csv(outdir / "chembl_similarity_hits.csv", hit_rows, hit_fields)
    write_csv(outdir / "chembl_similarity_activities.csv", activity_rows, activity_fields)
    print(
        f"queries={len(queries)} clusters={len(representatives)} hits={len(hit_rows)} "
        f"activities={len(activity_rows)} wrote={outdir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
