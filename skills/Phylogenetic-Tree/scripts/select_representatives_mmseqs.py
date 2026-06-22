#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


CONDA_BIN = Path("/mnt/l/WSL/softwares/conda_envs/md/bin")

KINGDOM_KEYWORDS = [
    ("Bacteria", ["bacteria", "bacterium", "escherichia", "pseudomonas", "staphylococcus", "xanthomonas", "bacillus", "streptomyces", "clostridium", "salmonella", "klebsiella", "mycobacterium", "lactobacillus", "propionibacterium", "sporosarcina", "rhizobium"]),
    ("Archaea", ["archaea", "archaeon", "methano", "halobacter"]),
    ("Fungi", ["fungi", "fungus", "saccharomyces", "aspergillus", "candida", "neurospora", "fusarium", "penicillium", "steccherinum"]),
    ("Plant", ["viridiplantae", "plantae", "arabidopsis", "oryza sativa", "zea mays", "canavalia", "phaseolus", "artemisia"]),
    ("Animal", ["metazoa", "animalia", "homo sapiens", "rattus", "mus musculus", "drosophila", "danio", "xenopus", "bos taurus"]),
    ("Virus", ["virus", "viruses", "phage"]),
    ("Metagenome", ["metagenome", "environmental sample"]),
]


def infer_kingdom(text: str) -> str:
    lower = (text or "").lower()
    for kingdom, keys in KINGDOM_KEYWORDS:
        if any(key in lower for key in keys):
            return kingdom
    return "Unknown"


def normalize_kingdom(kingdom: str, organism: str = "", phylum: str = "") -> str:
    text = " ".join([kingdom or "", organism or "", phylum or ""]).lower()
    if any(key in text for key in ["streptophyta", "viridiplantae", "plantae", "embryophyta", "tracheophyta"]):
        return "Plant"
    if any(key in text for key in ["fungi", "ascomycota", "basidiomycota", "saccharomyces", "fusarium", "steccherinum"]):
        return "Fungi"
    if any(key in text for key in ["metazoa", "animalia", "chordata", "arthropoda", "mammalia", "homo sapiens", "rattus"]):
        return "Animal"
    if "bacteria" in text:
        return "Bacteria"
    if "archaea" in text:
        return "Archaea"
    if "virus" in text or "viruses" in text:
        return "Virus"
    return kingdom or infer_kingdom(" ".join([organism or "", phylum or ""]))


def clean_accession(raw: str) -> str:
    text = str(raw or "").strip()
    if not text:
        return ""
    if "|" in text:
        parts = [p for p in text.split("|") if p]
        for part in parts:
            if re.match(r"^[A-Z0-9]{6,10}(?:-\d+)?$", part):
                return part
    text = re.sub(r"^UniRef\d+_", "", text)
    text = re.sub(r"\.\d+$", "", text)
    match = re.search(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})", text)
    return match.group(1) if match else text.split()[0]


def sanitize_id(raw: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.:-]", "_", raw or "")
    return "_" + safe if safe and safe[0].isdigit() else safe


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header = ""
    seq: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header and seq:
                    records.append((header, re.sub(r"[^A-Za-z]", "", "".join(seq)).upper()))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
    if header and seq:
        records.append((header, re.sub(r"[^A-Za-z]", "", "".join(seq)).upper()))
    return records


def write_fasta(records: Sequence[Tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx:idx + 80] + "\n")


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def metadata_by_accession(paths: Sequence[Path]) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    for path in paths:
        for row in read_csv(path):
            keys = [
                row.get("accession", ""),
                row.get("seed_id", ""),
                row.get("uniprot_ids", "").split(";")[0].strip(),
            ]
            for key in keys:
                key = clean_accession(key)
                if key:
                    out.setdefault(key, row)
    return out


def parse_header_metadata(header: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if " OS=" in header:
        out["organism"] = header.split(" OS=", 1)[1].split(" OX=", 1)[0].strip()
    if " OX=" in header:
        out["taxid"] = header.split(" OX=", 1)[1].split()[0].strip()
    for key in ("kingdom", "organism", "length", "pdb", "entity", "chains", "uniprot"):
        match = re.search(rf"(?:^|\s){key}=([^=]+?)(?=\s+\w+=|$)", header)
        if match:
            out[key] = match.group(1).strip()
    return out


def dedupe_records(records: Sequence[Tuple[str, str]]) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    seen = set()
    for header, seq in records:
        acc = clean_accession(header.split()[0])
        key = acc or seq
        if key in seen:
            continue
        seen.add(key)
        out.append((header, seq))
    return out


def find_mmseqs() -> Optional[str]:
    return shutil.which("mmseqs") or str(CONDA_BIN / "mmseqs") if (CONDA_BIN / "mmseqs").exists() else None


def cluster_mmseqs(input_fasta: Path, outdir: Path, min_seq_id: float, coverage: float, threads: int) -> Dict[str, List[str]]:
    mmseqs = find_mmseqs()
    if not mmseqs:
        return {}
    output_base = outdir / "mmseqs_cluster"
    tmp_dir = outdir / "mmseqs_tmp"
    cmd = [
        mmseqs,
        "easy-linclust",
        str(input_fasta),
        str(output_base),
        str(tmp_dir),
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(coverage),
        "--cov-mode",
        "1",
        "--threads",
        str(threads),
        "--remove-tmp-files",
        "1",
    ]
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"MMseqs failed:\n{result.stderr[-4000:]}")
    tsv = outdir / "mmseqs_cluster_cluster.tsv"
    if not tsv.exists():
        tsv = outdir / "mmseqs_cluster.tsv"
    clusters: Dict[str, List[str]] = {}
    if tsv.exists():
        with tsv.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    clusters.setdefault(parts[0], []).append(parts[1])
    return clusters


def select_representatives(
    records: Sequence[Tuple[str, str]],
    core_ids: set,
    clusters: Dict[str, List[str]],
    max_rep: int,
) -> List[Tuple[str, str]]:
    id_to_record = {header.split()[0]: (header, seq) for header, seq in records}
    selected: List[Tuple[str, str]] = []
    selected_ids = set()

    for header, seq in records:
        seq_id = header.split()[0]
        if seq_id in core_ids or clean_accession(seq_id) in core_ids:
            selected.append((header, seq))
            selected_ids.add(seq_id)

    remaining = max(0, max_rep - len(selected))
    if clusters:
        sorted_clusters = sorted(clusters.items(), key=lambda item: len(item[1]), reverse=True)
        total_members = sum(len(members) for _, members in sorted_clusters) or 1
        for rep_id, members in sorted_clusters:
            if remaining <= 0:
                break
            candidates = [id_to_record[m] for m in members if m in id_to_record and m not in selected_ids]
            candidates.sort(key=lambda item: len(item[1]), reverse=True)
            quota = max(1, round(remaining * len(members) / total_members))
            quota = min(quota, remaining, len(candidates))
            if len(members) >= 100:
                quota = min(max(quota, 20), remaining, len(candidates))
            elif len(members) >= 20:
                quota = min(max(quota, 5), remaining, len(candidates))
            for header, seq in candidates[:quota]:
                if remaining <= 0:
                    break
                selected.append((header, seq))
                selected_ids.add(header.split()[0])
                remaining -= 1
        if remaining > 0:
            for header, seq in sorted(records, key=lambda item: len(item[1]), reverse=True):
                seq_id = header.split()[0]
                if seq_id in selected_ids:
                    continue
                selected.append((header, seq))
                selected_ids.add(seq_id)
                remaining -= 1
                if remaining <= 0:
                    break
    else:
        bins: Dict[int, List[Tuple[str, str]]] = {}
        for header, seq in records:
            seq_id = header.split()[0]
            if seq_id in selected_ids:
                continue
            bins.setdefault(len(seq) // 50, []).append((header, seq))
        while remaining > 0 and any(bins.values()):
            for key in sorted(bins):
                if remaining <= 0:
                    break
                if bins[key]:
                    header, seq = bins[key].pop(0)
                    selected.append((header, seq))
                    remaining -= 1
    return selected[:max_rep]


def write_nodes(
    path: Path,
    records: Sequence[Tuple[str, str]],
    core_ids: set,
    metadata: Dict[str, Dict[str, str]],
) -> None:
    fields = ["id", "uniprot_id", "species_code", "kingdom", "organism", "phylum_class", "length", "is_representative", "is_core", "community", "cluster_id"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for header, seq in records:
            seq_id = header.split()[0]
            acc = clean_accession(seq_id)
            md = metadata.get(acc, {})
            hm = parse_header_metadata(header)
            organism = md.get("organism", "") or md.get("Organism", "") or hm.get("organism", "")
            phylum = md.get("phylum", "") or md.get("phylum_class", "") or "Unknown"
            kingdom = normalize_kingdom(md.get("kingdom", "") or hm.get("kingdom", ""), organism, phylum)
            writer.writerow({
                "id": seq_id,
                "uniprot_id": acc,
                "species_code": "UNKNOWN",
                "kingdom": kingdom or "Unknown",
                "organism": organism or "Unknown",
                "phylum_class": phylum,
                "length": str(len(seq)),
                "is_representative": "TRUE",
                "is_core": "TRUE" if (seq_id in core_ids or acc in core_ids) else "FALSE",
                "community": "-1",
                "cluster_id": "",
            })


def main() -> None:
    parser = argparse.ArgumentParser(description="General MMseqs representative selector without hydrolase length filters.")
    parser.add_argument("input_fasta")
    parser.add_argument("--core-fasta", default="")
    parser.add_argument("--output", default="/mnt/e/Tree-Metal-F/ssn")
    parser.add_argument("--metadata", nargs="*", default=["/mnt/e/Tree-Metal-F/hmmer/homolog_metadata.csv", "/mnt/e/Tree-Metal-F/tree_core_metadata.csv"])
    parser.add_argument("--max-rep", type=int, default=200)
    parser.add_argument("--min-seq-id", type=float, default=0.35)
    parser.add_argument("--coverage", type=float, default=0.75)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    records = dedupe_records(read_fasta(Path(args.input_fasta)))
    core_records = read_fasta(Path(args.core_fasta)) if args.core_fasta else []
    core_ids = {header.split()[0] for header, _ in core_records} | {clean_accession(header.split()[0]) for header, _ in core_records}
    cluster_input = outdir / "cluster_input.fasta"
    write_fasta(records, cluster_input)
    clusters = cluster_mmseqs(cluster_input, outdir, args.min_seq_id, args.coverage, args.threads)
    reps = select_representatives(records, core_ids, clusters, args.max_rep)
    write_fasta(reps, outdir / "representatives.fasta")
    md = metadata_by_accession([Path(p) for p in args.metadata])
    write_nodes(outdir / "ssn_nodes.csv", reps, core_ids, md)
    print(f"Input records: {len(records)}")
    print(f"Clusters: {len(clusters) if clusters else 0}")
    print(f"Representatives: {len(reps)}")
    print(f"Core kept: {sum(1 for h, _ in reps if h.split()[0] in core_ids or clean_accession(h.split()[0]) in core_ids)}/{len(core_records)}")
    print(f"Output: {outdir}")


if __name__ == "__main__":
    main()
