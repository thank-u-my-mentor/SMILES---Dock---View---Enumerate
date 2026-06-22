#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path

import pandas as pd


MOTIFS = ["GTTDDS", "APNNGLL", "FADAG"]
try:
    from Bio.Align import PairwiseAligner
except Exception:
    PairwiseAligner = None


def clean_seq(value: str) -> str:
    return re.sub(r"[^A-Za-z]", "", str(value or "")).upper()


def clean_text(value: str) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def write_fasta(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            acc = clean_text(row.get("accession"))
            seq = clean_seq(row.get("sequence"))
            if not acc or not seq:
                continue
            level = clean_text(row.get("level")).replace(" ", "_").replace("/", "-")
            org = clean_text(row.get("organism")) or "unknown"
            name = clean_text(row.get("protein_name"))[:120].replace("|", ";")
            rank = clean_text(row.get("rank"))
            handle.write(f">{acc}|rank={rank}|level={level}|organism={org}|name={name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def read_fasta(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    header = ""
    parts: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                add_record(records, header, "".join(parts))
            header = line[1:]
            parts = []
        else:
            parts.append(line)
    if header:
        add_record(records, header, "".join(parts))
    return records


def clean_accession(header: str) -> str:
    first = str(header or "").split()[0]
    if "|" in first:
        for part in first.split("|"):
            if re.match(r"^[A-Z0-9_]+(?:\.\d+)?(?:-\d+)?$", part):
                return part
    return first


def add_record(records: dict[str, dict[str, str]], header: str, sequence: str) -> None:
    acc = clean_accession(header)
    records[acc] = {"accession": acc, "header": header, "sequence": clean_seq(sequence)}


def motif_positions(sequence: str) -> dict[str, int]:
    seq = clean_seq(sequence)
    return {motif: seq.find(motif) for motif in MOTIFS}


def motif_score(sequence: str) -> float:
    pos = motif_positions(sequence)
    return round(100.0 * sum(1 for value in pos.values() if value >= 0) / len(MOTIFS), 2)


def positional_identity(a: str, b: str) -> float:
    a = clean_seq(a)
    b = clean_seq(b)
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    matches = sum(x == y for x, y in zip(a[:n], b[:n]))
    return matches / max(len(a), len(b))


def aligned_identity(a: str, b: str) -> float:
    a = clean_seq(a)
    b = clean_seq(b)
    if not a or not b:
        return 0.0
    if PairwiseAligner is None:
        return positional_identity(a, b)
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -0.1
    aligner.extend_gap_score = -0.01
    aln = aligner.align(a, b)[0]
    top = str(aln[0])
    bottom = str(aln[1])
    matches = 0
    aligned_pairs = 0
    for x, y in zip(top, bottom):
        if x == "-" or y == "-":
            continue
        aligned_pairs += 1
        if x == y:
            matches += 1
    return matches / aligned_pairs if aligned_pairs else 0.0


def best_seed_identity(sequence: str, seeds: dict[str, dict[str, str]]) -> tuple[str, float]:
    best_acc = ""
    best_ident = -1.0
    for acc, row in seeds.items():
        ident = aligned_identity(sequence, row["sequence"])
        if ident > best_ident:
            best_ident = ident
            best_acc = acc
    return best_acc, round(best_ident, 6)


def infer_level(row: dict) -> str:
    length = int(row.get("length") or 0)
    motif = float(row.get("motif_active_site_score") or 0)
    ident = float(row.get("nearest_seed_identity") or 0)
    evalue = row.get("hmmsearch_evalue")
    try:
        evalue_f = float(evalue)
    except Exception:
        evalue_f = 1.0
    if motif >= 99 and 280 <= length <= 305 and evalue_f <= 1e-80 and ident >= 0.55:
        return "Level 2 high-confidence candidate"
    if motif >= 66 and 260 <= length <= 330 and evalue_f <= 1e-30:
        return "Level 3 motif-partial candidate"
    return "Review/low-confidence homolog"


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "_No rows._"
    cols = list(frame.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for _, row in frame.iterrows():
        values = [clean_text(row.get(col, "")).replace("|", "\\|") for col in cols]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def command_prepare(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.input_csv).fillna("")
    rows = frame.to_dict("records")
    write_fasta(rows, outdir / "seed.fasta")
    frame.to_csv(outdir / "seed_metadata.csv", index=False)
    print(f"seed_records={len(rows)}")
    print(f"seed_fasta={outdir / 'seed.fasta'}")
    print(f"seed_metadata={outdir / 'seed_metadata.csv'}")


def command_score(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    seeds = read_fasta(Path(args.seed_fasta))
    seqs = read_fasta(Path(args.fasta))
    hmmer = {}
    if args.hmmer_hits and Path(args.hmmer_hits).exists():
        with open(args.hmmer_hits, newline="", encoding="utf-8", errors="replace") as handle:
            hmmer = {row.get("accession", ""): row for row in csv.DictReader(handle)}
    metadata = {}
    if args.metadata and Path(args.metadata).exists():
        with open(args.metadata, newline="", encoding="utf-8-sig", errors="replace") as handle:
            metadata = {row.get("accession", ""): row for row in csv.DictReader(handle)}

    rows = []
    for rank, (acc, rec) in enumerate(seqs.items(), 1):
        seq = rec["sequence"]
        nearest, ident = best_seed_identity(seq, seeds)
        hit = hmmer.get(acc, {})
        meta = metadata.get(acc, {})
        pos = motif_positions(seq)
        row = {
            "rank": rank,
            "accession": acc,
            "protein_name": hit.get("description", "") or meta.get("protein_name", "") or meta.get("description", "") or meta.get("header", ""),
            "organism": meta.get("organism", "") or hit.get("organism", ""),
            "kingdom": meta.get("kingdom", ""),
            "length": len(seq),
            "nearest_seed": nearest,
            "nearest_seed_identity": ident,
            "hmmsearch_evalue": hit.get("evalue", "") or meta.get("hmmsearch_evalue", ""),
            "hmmsearch_bitscore": hit.get("score", "") or meta.get("hmmsearch_bitscore", ""),
            "motif_active_site_score": motif_score(seq),
            "motif_positions": pos,
            "md5": hashlib.md5(seq.encode()).hexdigest(),
            "sha256": hashlib.sha256(seq.encode()).hexdigest(),
            "sequence": seq,
        }
        row["level"] = "Seed/from input CSV" if acc in seeds else infer_level(row)
        rows.append(row)

    rows.sort(
        key=lambda row: (
            0 if str(row["level"]).startswith("Seed") else 1,
            -float(row["motif_active_site_score"]),
            -float(row["nearest_seed_identity"]),
            row["accession"],
        )
    )
    for idx, row in enumerate(rows, 1):
        row["rank"] = idx

    fields = [
        "rank",
        "accession",
        "level",
        "protein_name",
        "organism",
        "kingdom",
        "length",
        "nearest_seed",
        "nearest_seed_identity",
        "hmmsearch_evalue",
        "hmmsearch_bitscore",
        "motif_active_site_score",
        "motif_positions",
        "md5",
        "sha256",
        "sequence",
    ]
    out_csv = outdir / "fluorinase_hmmer_candidates_scored.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    high = [row for row in rows if str(row["level"]).startswith(("Seed", "Level 2"))]
    write_fasta(high, outdir / "fluorinase_high_confidence_plus_seed.fasta")
    print(f"scored_records={len(rows)}")
    print(f"high_confidence_plus_seed={len(high)}")
    print(f"scored_csv={out_csv}")
    print(f"high_confidence_fasta={outdir / 'fluorinase_high_confidence_plus_seed.fasta'}")


def command_summarize(args: argparse.Namespace) -> None:
    base = Path(args.outdir)
    seed = pd.read_csv(base / "seed_metadata.csv").fillna("")
    scored = pd.read_csv(base / "fluorinase_hmmer_candidates_scored.csv").fillna("")
    strict_hits = pd.read_csv(base / "hmmer" / "hmmer_hits.csv").fillna("")
    pairs = pd.read_csv(base / "tree_analysis" / "pairwise_identity_pairs_both_nongap.csv").fillna("")
    pairs_ge90 = pd.read_csv(base / "tree_analysis" / "pairwise_identity_pairs_ge90_both_nongap.csv").fillna("")
    broad_hits_path = base / "hmmer_e20" / "hmmer_hits.csv"
    broad_hits = pd.read_csv(broad_hits_path).fillna("") if broad_hits_path.exists() else pd.DataFrame()

    seed_accessions = set(seed["accession"].astype(str))
    scored["is_original_seed"] = scored["accession"].astype(str).isin(seed_accessions)
    new = scored[~scored["is_original_seed"]].copy()
    strict_new = base / "strict_hmmer_new_nonseed_candidates.csv"
    new_cols = [
        "accession",
        "level",
        "protein_name",
        "organism",
        "kingdom",
        "length",
        "nearest_seed",
        "nearest_seed_identity",
        "motif_active_site_score",
        "hmmsearch_evalue",
        "sequence",
    ]
    new[[col for col in new_cols if col in new.columns]].to_csv(strict_new, index=False)

    high_conf = scored[
        scored["level"].astype(str).str.startswith("Seed")
        | scored["level"].astype(str).str.startswith("Level 2")
    ]

    lines = [
        "# Nature Fluorine Enzyme HMMER Summary",
        "",
        f"- Original seed rows: {len(seed_accessions)}",
        f"- Strict HMMER 1e-50 parsed hits: {len(strict_hits)}",
        f"- Strict HMMER scored sequences: {len(scored)}",
        f"- New non-seed sequences in strict scored set: {len(new)}",
        f"- High-confidence plus seed FASTA records: {len(high_conf)}",
        "- Identity matrix records: 22",
        f"- Pairwise comparisons: {len(pairs)}",
        f"- Pairs with identity >=90%: {len(pairs_ge90)}",
    ]
    if not pairs.empty:
        vals = pd.to_numeric(pairs["identity_percent"], errors="coerce").dropna()
        lines.append(f"- Min/median/max pairwise identity: {vals.min():.2f}% / {vals.median():.2f}% / {vals.max():.2f}%")
    lines += ["", "## Strict HMMER New Non-Seed Sequences", ""]
    show_cols = [col for col in new_cols[:-1] if col in new.columns]
    lines.append(markdown_table(new[show_cols]))
    if not broad_hits.empty:
        lines += [
            "",
            "## Broad HMMER 1e-20 Warning",
            "",
            "The broader 1e-20 search produced 1000 parsed hits, but the top 100 already contains RNA polymerase delta, molybdate/tungstate binding protein, adenosyl-chloride synthase, and other distant homolog/noise. Use motif+length filtering before interpreting these as C-F enzymes.",
            "",
            "Top description counts among first 100 broad hits:",
            "",
            markdown_table(
                broad_hits.head(100)["description"]
                .value_counts()
                .head(20)
                .rename_axis("description")
                .reset_index(name="count")
            ),
        ]
    summary = base / "README_fluorinase_hmmer_summary.md"
    summary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"summary={summary}")
    print(f"strict_new_candidates={strict_new}")


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)
    p_prepare = sub.add_parser("prepare")
    p_prepare.add_argument("--input-csv", required=True)
    p_prepare.add_argument("--outdir", required=True)
    p_prepare.set_defaults(func=command_prepare)

    p_score = sub.add_parser("score")
    p_score.add_argument("--seed-fasta", required=True)
    p_score.add_argument("--fasta", required=True)
    p_score.add_argument("--hmmer-hits", default="")
    p_score.add_argument("--metadata", default="")
    p_score.add_argument("--outdir", required=True)
    p_score.set_defaults(func=command_score)

    p_summary = sub.add_parser("summarize")
    p_summary.add_argument("--outdir", required=True)
    p_summary.set_defaults(func=command_summarize)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
