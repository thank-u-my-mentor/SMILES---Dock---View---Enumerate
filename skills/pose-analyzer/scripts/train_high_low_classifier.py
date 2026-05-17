#!/usr/bin/env python
"""Train high-vs-low molecule classifiers from score_space_model outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")
rdBase.DisableLog("rdApp.*")


def read_csv(path: Path) -> list[dict[str, str]]:
    for encoding in ("utf-8-sig", "utf-8", "gb18030", "cp936", "latin1"):
        try:
            with path.open("r", encoding=encoding, newline="") as handle:
                return list(csv.DictReader(handle))
        except UnicodeDecodeError:
            continue
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def safe_float(value: object) -> float:
    try:
        if value in ("", None):
            return math.nan
        return float(value)
    except Exception:
        return math.nan


def mol_from_smiles(smiles: str) -> Chem.Mol | None:
    try:
        return Chem.MolFromSmiles(str(smiles or "").strip())
    except Exception:
        return None


def canonical(smiles: str) -> str:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def morgan_vector(smiles: str, n_bits: int = 2048) -> np.ndarray:
    arr = np.zeros((n_bits,), dtype=np.float32)
    mol = mol_from_smiles(smiles)
    if mol is None:
        return arr
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def tanimoto_to_set(smiles: str, fps: list[object]) -> float:
    mol = mol_from_smiles(smiles)
    if mol is None or not fps:
        return math.nan
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    return float(max(DataStructs.BulkTanimotoSimilarity(fp, fps)))


def fingerprint(smiles: str):
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def label_rows(rows: list[dict[str, str]], high_threshold: float, low_max: float) -> tuple[list[int], list[int]]:
    indexes: list[int] = []
    labels: list[int] = []
    for i, row in enumerate(rows):
        if str(row.get("smiles_sanity_status", "")).lower() == "reject":
            continue
        score = safe_float(row.get("official_binding_score"))
        if math.isnan(score):
            continue
        if score >= high_threshold:
            indexes.append(i)
            labels.append(1)
        elif 0.0 <= score <= low_max:
            indexes.append(i)
            labels.append(0)
    return indexes, labels


def classification_metrics(y_true: list[int], probs: list[float]) -> dict[str, float | str]:
    if not y_true:
        return {}
    preds = [1 if p >= 0.5 else 0 for p in probs]
    out: dict[str, float | str] = {
        "accuracy": round(float(accuracy_score(y_true, preds)), 6),
        "balanced_accuracy": round(float(balanced_accuracy_score(y_true, preds)), 6),
    }
    try:
        if len(set(y_true)) == 2:
            out["roc_auc"] = round(float(roc_auc_score(y_true, probs)), 6)
    except Exception:
        out["roc_auc"] = ""
    return out


def train_morgan_models(rows: list[dict[str, str]], labeled_indexes: list[int], labels: list[int], seed: int) -> tuple[dict[str, object], list[float]]:
    x_label = np.vstack([morgan_vector(rows[i].get("canonical_smiles", "")) for i in labeled_indexes])
    y = np.array(labels, dtype=int)
    x_all = np.vstack([morgan_vector(row.get("canonical_smiles", "")) for row in rows])
    min_class = min(int(np.sum(y == 0)), int(np.sum(y == 1))) if len(set(y)) == 2 else 0
    n_splits = min(5, min_class)
    oof = np.full((len(y),), np.nan, dtype=float)
    fold_metrics: list[dict[str, object]] = []
    if n_splits >= 2:
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        for fold, (train_idx, test_idx) in enumerate(skf.split(x_label, y), start=1):
            rf = RandomForestClassifier(n_estimators=300, min_samples_leaf=2, class_weight="balanced", random_state=seed + fold)
            rf.fit(x_label[train_idx], y[train_idx])
            probs = rf.predict_proba(x_label[test_idx])[:, 1]
            oof[test_idx] = probs
            metric = classification_metrics(y[test_idx].tolist(), probs.tolist())
            metric["fold"] = fold
            fold_metrics.append(metric)
    rf_final = RandomForestClassifier(n_estimators=500, min_samples_leaf=2, class_weight="balanced", random_state=seed)
    rf_final.fit(x_label, y)
    rf_probs = rf_final.predict_proba(x_all)[:, 1].tolist()
    logit = make_pipeline(
        StandardScaler(with_mean=False),
        LogisticRegression(max_iter=2000, class_weight="balanced", random_state=seed),
    )
    logit.fit(x_label, y)
    logit_probs = logit.predict_proba(x_all)[:, 1].tolist()
    probs = [0.7 * a + 0.3 * b for a, b in zip(rf_probs, logit_probs)]
    metrics: dict[str, object] = {
        "model": "morgan_rf_plus_logistic",
        "labeled_rows": len(labels),
        "high_rows": int(np.sum(y == 1)),
        "low_rows": int(np.sum(y == 0)),
        "cv_folds": n_splits,
        "fold_metrics": fold_metrics,
    }
    if np.isfinite(oof).all():
        metrics["oof_metrics"] = classification_metrics(labels, oof.tolist())
    return metrics, probs


ATOM_NUMS = [6, 7, 8, 9, 15, 16, 17, 35, 53]


def atom_features(atom: Chem.Atom) -> list[float]:
    atomic = [1.0 if atom.GetAtomicNum() == num else 0.0 for num in ATOM_NUMS]
    other = [0.0 if atom.GetAtomicNum() in ATOM_NUMS else 1.0]
    return [
        *atomic,
        *other,
        atom.GetDegree() / 4.0,
        atom.GetFormalCharge(),
        1.0 if atom.GetIsAromatic() else 0.0,
        1.0 if atom.IsInRing() else 0.0,
        atom.GetTotalNumHs() / 4.0,
    ]


def graph_from_smiles(smiles: str) -> tuple[np.ndarray, np.ndarray] | None:
    mol = mol_from_smiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    x = np.array([atom_features(atom) for atom in mol.GetAtoms()], dtype=np.float32)
    n = mol.GetNumAtoms()
    adj = np.eye(n, dtype=np.float32)
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        adj[i, j] = 1.0
        adj[j, i] = 1.0
    deg = np.sum(adj, axis=1)
    inv = np.diag(1.0 / np.sqrt(np.maximum(deg, 1e-6)))
    adj = inv @ adj @ inv
    return x, adj.astype(np.float32)


def train_graph_model(rows: list[dict[str, str]], labeled_indexes: list[int], labels: list[int], seed: int, epochs: int) -> tuple[dict[str, object], list[float]]:
    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
    except Exception as exc:
        return {"model": "graph_gcn", "status": f"unavailable:{exc.__class__.__name__}"}, [math.nan for _ in rows]

    graphs_all = [graph_from_smiles(row.get("canonical_smiles", "")) for row in rows]
    usable = [(idx, labels[pos]) for pos, idx in enumerate(labeled_indexes) if graphs_all[idx] is not None]
    if len(usable) < 8 or len({label for _, label in usable}) < 2:
        return {"model": "graph_gcn", "status": "not_enough_labeled_graphs"}, [math.nan for _ in rows]

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    in_dim = len(atom_features(Chem.MolFromSmiles("C").GetAtomWithIdx(0)))

    class TinyGCN(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.lin1 = nn.Linear(in_dim, 48)
            self.lin2 = nn.Linear(48, 48)
            self.out = nn.Linear(48, 1)

        def forward(self, x, adj):
            h = F.relu(adj @ self.lin1(x))
            h = F.relu(adj @ self.lin2(h))
            pooled = h.mean(dim=0)
            return self.out(pooled).squeeze()

    def tensor_graph(idx: int):
        graph = graphs_all[idx]
        assert graph is not None
        x, adj = graph
        return torch.tensor(x), torch.tensor(adj)

    model = TinyGCN()
    labels_arr = np.array([label for _, label in usable], dtype=int)
    pos_weight = torch.tensor([max(1.0, float(np.sum(labels_arr == 0)) / max(1, float(np.sum(labels_arr == 1))))])
    opt = torch.optim.Adam(model.parameters(), lr=0.006, weight_decay=1e-4)
    train_pairs = usable
    for _ in range(epochs):
        random.shuffle(train_pairs)
        for idx, label in train_pairs:
            x, adj = tensor_graph(idx)
            y = torch.tensor(float(label))
            loss = nn.BCEWithLogitsLoss(pos_weight=pos_weight)(model(x, adj).view(1), y.view(1))
            opt.zero_grad()
            loss.backward()
            opt.step()
    probs: list[float] = []
    model.eval()
    with torch.no_grad():
        for graph in graphs_all:
            if graph is None:
                probs.append(math.nan)
                continue
            x_np, adj_np = graph
            prob = torch.sigmoid(model(torch.tensor(x_np), torch.tensor(adj_np))).item()
            probs.append(float(prob))
    metrics = {
        "model": "tiny_graph_gcn",
        "status": "ok",
        "labeled_graphs": len(usable),
        "note": "Small-data neural baseline; use as a structural hypothesis signal, not as a final QSAR model.",
    }
    return metrics, probs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-matrix", type=Path, default=Path("~/vina_task2/score_space_model/binding_score_feature_matrix.csv"))
    parser.add_argument("--outdir", type=Path, default=Path("~/vina_task2/high_low_classifier"))
    parser.add_argument("--high-score-threshold", type=float, default=0.4)
    parser.add_argument("--low-score-max", type=float, default=0.25)
    parser.add_argument("--gnn-epochs", type=int, default=180)
    parser.add_argument("--seed", type=int, default=13)
    args = parser.parse_args()

    rows = read_csv(args.feature_matrix.expanduser())
    labeled_indexes, labels = label_rows(rows, args.high_score_threshold, args.low_score_max)
    if len(set(labels)) < 2:
        raise SystemExit("Need both high and low labeled molecules.")
    morgan_metrics, morgan_probs = train_morgan_models(rows, labeled_indexes, labels, args.seed)
    gnn_metrics, gnn_probs = train_graph_model(rows, labeled_indexes, labels, args.seed, args.gnn_epochs)

    high_fps = [fingerprint(rows[i].get("canonical_smiles", "")) for i, label in zip(labeled_indexes, labels) if label == 1]
    low_fps = [fingerprint(rows[i].get("canonical_smiles", "")) for i, label in zip(labeled_indexes, labels) if label == 0]
    high_fps = [fp for fp in high_fps if fp is not None]
    low_fps = [fp for fp in low_fps if fp is not None]

    pred_rows: list[dict[str, object]] = []
    for row, p_morgan, p_gnn in zip(rows, morgan_probs, gnn_probs):
        score = safe_float(row.get("official_binding_score"))
        label = ""
        if not math.isnan(score):
            if score >= args.high_score_threshold:
                label = "high"
            elif 0.0 <= score <= args.low_score_max:
                label = "low"
            else:
                label = "middle"
        p_final = p_morgan if math.isnan(p_gnn) else 0.65 * p_morgan + 0.35 * p_gnn
        smiles = row.get("canonical_smiles", "")
        pred_rows.append(
            {
                "seq_id": row.get("seq_id", ""),
                "nickname": row.get("nickname", ""),
                "label": label,
                "official_binding_score": row.get("official_binding_score", ""),
                "prob_high_morgan": round(float(p_morgan), 6),
                "prob_high_gnn": round(float(p_gnn), 6) if not math.isnan(p_gnn) else "",
                "prob_high_ensemble": round(float(p_final), 6),
                "nearest_high_tanimoto": round(tanimoto_to_set(smiles, high_fps), 6) if high_fps else "",
                "nearest_low_tanimoto": round(tanimoto_to_set(smiles, low_fps), 6) if low_fps else "",
                "smiles_sanity_status": row.get("smiles_sanity_status", ""),
                "smiles_sanity_reasons": row.get("smiles_sanity_reasons", ""),
                "canonical_smiles": smiles,
            }
        )
    pred_rows.sort(key=lambda row: safe_float(row.get("prob_high_ensemble")), reverse=True)
    outdir = args.outdir.expanduser().resolve()
    write_csv(
        outdir / "high_low_classifier_predictions.csv",
        pred_rows,
        [
            "seq_id", "nickname", "label", "official_binding_score",
            "prob_high_morgan", "prob_high_gnn", "prob_high_ensemble",
            "nearest_high_tanimoto", "nearest_low_tanimoto",
            "smiles_sanity_status", "smiles_sanity_reasons", "canonical_smiles",
        ],
    )
    metrics = {
        "high_threshold": args.high_score_threshold,
        "low_max": args.low_score_max,
        "morgan": morgan_metrics,
        "gnn": gnn_metrics,
        "output_note": "Use prob_high_ensemble to prioritize generated/docked unscored molecules that structurally resemble the high-score group.",
    }
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "high_low_classifier_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
    print(
        f"labeled={len(labels)} high={sum(labels)} low={len(labels)-sum(labels)} "
        f"wrote={outdir / 'high_low_classifier_predictions.csv'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
