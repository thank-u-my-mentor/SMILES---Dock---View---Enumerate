#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence
from urllib.parse import urljoin


DEFAULT_API_BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
USER_AGENT = "Phylogenetic-Tree/0.1"

KINGDOM_KEYWORDS = [
    ("Bacteria", ["bacteria", "bacterium", "escherichia", "pseudomonas", "staphylococcus", "xanthomonas", "bacillus", "streptomyces", "clostridium", "salmonella", "klebsiella", "mycobacterium", "lactobacillus", "propionibacterium", "sporosarcina", "rhizobium"]),
    ("Archaea", ["archaea", "archaeon", "methano", "halobacter"]),
    ("Fungi", ["fungi", "fungus", "saccharomyces", "aspergillus", "candida", "neurospora", "fusarium", "penicillium", "steccherinum"]),
    ("Plant", ["viridiplantae", "plantae", "arabidopsis", "oryza sativa", "zea mays", "canavalia", "phaseolus"]),
    ("Animal", ["metazoa", "animalia", "homo sapiens", "rattus", "mus musculus", "drosophila", "danio", "xenopus", "bos taurus"]),
    ("Virus", ["virus", "viruses", "phage"]),
    ("Metagenome", ["metagenome", "environmental sample"]),
]


@dataclass
class FastaRecord:
    header: str
    sequence: str

    @property
    def id(self) -> str:
        return self.header.split()[0]


@dataclass
class Hit:
    accession: str
    evalue: Optional[float] = None
    score: Optional[float] = None
    description: str = ""


def infer_kingdom(text: str) -> str:
    lower = (text or "").lower()
    for kingdom, keys in KINGDOM_KEYWORDS:
        if any(key in lower for key in keys):
            return kingdom
    return "Unknown"


def clean_accession(raw: str) -> str:
    text = str(raw or "").strip()
    if "|" in text:
        parts = [p for p in text.split("|") if p]
        for part in parts:
            if re.match(r"^[A-Z0-9]{6,10}(?:-\d+)?$", part):
                return part
    text = re.sub(r"^UniRef\d+_", "", text)
    text = re.sub(r"\.\d+$", "", text)
    match = re.search(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})", text)
    return match.group(1) if match else text.split()[0]


def maybe_float(value) -> Optional[float]:
    try:
        if value in (None, ""):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def read_fasta(path: Path, strip_gaps: bool = True) -> List[FastaRecord]:
    records: List[FastaRecord] = []
    header: Optional[str] = None
    seq: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header and seq:
                    sequence = "".join(seq)
                    if strip_gaps:
                        sequence = re.sub(r"[-.\s]", "", sequence)
                    records.append(FastaRecord(header, re.sub(r"[^A-Za-z]", "", sequence).upper()))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
    if header and seq:
        sequence = "".join(seq)
        if strip_gaps:
            sequence = re.sub(r"[-.\s]", "", sequence)
        records.append(FastaRecord(header, re.sub(r"[^A-Za-z]", "", sequence).upper()))
    return records


def write_fasta(records: Iterable[FastaRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for rec in records:
            handle.write(f">{rec.header}\n")
            seq = rec.sequence.replace(" ", "").replace("\n", "")
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx:idx + 80] + "\n")


def run_alignment(input_fasta: Path, aligned_fasta: Path, threads: int) -> None:
    clustalo = shutil.which("clustalo") or "/mnt/l/WSL/softwares/conda_envs/md/bin/clustalo"
    if not Path(clustalo).exists():
        raise RuntimeError("clustalo was not found; provide --aligned-input or add clustalo to PATH.")
    aligned_fasta.parent.mkdir(parents=True, exist_ok=True)
    cmd = [clustalo, "-i", str(input_fasta), "-o", str(aligned_fasta), "--outfmt", "fasta", "--force", "--threads", str(threads)]
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"clustalo failed:\n{result.stderr[-3000:]}")


def submit_hmmsearch(
    alignment_text: str,
    api_base: str,
    database: str,
    evalue: str,
    domain_evalue: str,
    report_evalue: str,
    report_domain_evalue: str,
    timeout: int,
) -> str:
    import requests

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})
    payload = {
        "database": database,
        "input": alignment_text,
        "incE": str(evalue),
        "incdomE": str(domain_evalue),
        "E": str(report_evalue),
        "domE": str(report_domain_evalue),
    }
    url = urljoin(api_base.rstrip("/") + "/", "search/hmmsearch")
    response = session.post(url, json=payload, timeout=timeout)
    if response.status_code not in {200, 201, 202}:
        raise RuntimeError(f"EBI HMMER hmmsearch failed: HTTP {response.status_code}\n{response.text[:1000]}")
    data = response.json()
    job_id = data.get("id") or data.get("jobId") or data.get("uuid")
    if not job_id:
        raise RuntimeError(f"Could not find HMMER job id in response: {data}")
    return str(job_id)


def poll_result(api_base: str, job_id: str, timeout: int, poll_seconds: int, max_polls: int) -> Dict:
    import requests

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})
    url = urljoin(api_base.rstrip("/") + "/", f"result/{job_id}")
    last_text = ""
    for _ in range(max_polls):
        response = session.get(url, timeout=timeout)
        last_text = response.text[:1000]
        if response.status_code in {202, 204}:
            time.sleep(poll_seconds)
            continue
        if response.status_code != 200:
            raise RuntimeError(f"HMMER result failed: HTTP {response.status_code}\n{last_text}")
        data = response.json()
        status = str(data.get("status", "")).upper()
        if status in {"PENDING", "RUNNING", "QUEUED"}:
            time.sleep(poll_seconds)
            continue
        return data
    raise TimeoutError(f"Timed out waiting for HMMER job {job_id}: {last_text}")


def extract_hits(payload: object) -> List[Hit]:
    hits: Dict[str, Hit] = {}

    def walk(obj: object) -> None:
        if isinstance(obj, dict):
            keys = {str(k).lower(): k for k in obj.keys()}
            acc_key = next((keys[k] for k in ("acc", "accession", "target", "name") if k in keys), None)
            if acc_key is not None:
                acc = clean_accession(str(obj.get(acc_key, "")))
                if acc and acc.lower() not in {"none", "unknown", "query"}:
                    evalue = maybe_float(obj.get(keys.get("evalue", "")) or obj.get(keys.get("e-value", "")))
                    score = maybe_float(obj.get(keys.get("score", "")) or obj.get(keys.get("bits", "")) or obj.get(keys.get("bit_score", "")))
                    desc = str(obj.get(keys.get("desc", ""), "") or obj.get(keys.get("description", ""), ""))
                    old = hits.get(acc)
                    if old is None or ((evalue if evalue is not None else 1e99) < (old.evalue if old.evalue is not None else 1e99)):
                        hits[acc] = Hit(acc, evalue, score, desc)
            for value in obj.values():
                walk(value)
        elif isinstance(obj, list):
            for value in obj:
                walk(value)

    walk(payload)
    return sorted(hits.values(), key=lambda h: (h.evalue if h.evalue is not None else 1e99, -(h.score or 0)))


def download_uniprot_fastas(accessions: Sequence[str], retries: int = 2) -> List[FastaRecord]:
    import requests

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    records: List[FastaRecord] = []
    for idx, accession in enumerate(accessions, 1):
        if idx % 100 == 0:
            print(f"Downloaded {len(records)}/{idx - 1} UniProt FASTA records...")
        for attempt in range(retries + 1):
            response = session.get(UNIPROT_FASTA.format(accession=accession), timeout=30)
            if response.status_code == 200 and response.text.startswith(">"):
                tmp = Path("__tmp_uniprot_fetch.fasta")
                tmp.write_text(response.text, encoding="utf-8")
                try:
                    records.extend(read_fasta(tmp))
                finally:
                    tmp.unlink(missing_ok=True)
                break
            if response.status_code in {400, 404}:
                break
            time.sleep(1.0 * (attempt + 1))
    return records


def fasta_metadata(records: Sequence[FastaRecord]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for rec in records:
        header = rec.header
        accession = clean_accession(rec.id)
        organism = ""
        taxid = ""
        if " OS=" in header:
            organism = header.split(" OS=", 1)[1].split(" OX=", 1)[0].strip()
        if " OX=" in header:
            taxid = header.split(" OX=", 1)[1].split()[0].strip()
        rows.append({
            "accession": accession,
            "organism": organism,
            "taxid": taxid,
            "kingdom": infer_kingdom(organism),
            "length": str(len(rec.sequence)),
            "header": header,
        })
    return rows


def dedupe_records(records: Sequence[FastaRecord], max_total: int) -> List[FastaRecord]:
    out: List[FastaRecord] = []
    seen = set()
    for rec in records:
        key = clean_accession(rec.id) or rec.sequence
        if key in seen:
            continue
        seen.add(key)
        out.append(rec)
        if len(out) >= max_total:
            break
    return out


def write_hits_csv(hits: Sequence[Hit], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["accession", "evalue", "score", "description"])
        for hit in hits:
            writer.writerow([hit.accession, hit.evalue if hit.evalue is not None else "", hit.score if hit.score is not None else "", hit.description])


def write_metadata_csv(rows: Sequence[Dict[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "organism", "taxid", "kingdom", "length", "header"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run EBI HMMER hmmsearch from an aligned protein FASTA/HMM input.")
    parser.add_argument("input_fasta", help="Seed FASTA. It will be aligned unless --aligned-input is set.")
    parser.add_argument("--outdir", default="/mnt/e/Tree-Metal-F/hmmer")
    parser.add_argument("--aligned-input", action="store_true", help="Treat input_fasta as an already aligned FASTA/HMMER-compatible MSA.")
    parser.add_argument("--database", default="uniprot", help="HMMER sequence database, e.g. uniprot, refprot, swissprot, pdb, rp15/rp35/rp55/rp75.")
    parser.add_argument("--evalue", default="1e-10", help="Sequence inclusion E-value; use 1e-10 to 1e-20 for strict searches.")
    parser.add_argument("--domain-evalue", default=None, help="Domain inclusion E-value; defaults to --evalue.")
    parser.add_argument("--report-evalue", default=None, help="Sequence reporting E-value; defaults to --evalue.")
    parser.add_argument("--report-domain-evalue", default=None, help="Domain reporting E-value; defaults to --domain-evalue.")
    parser.add_argument("--max-hits", type=int, default=1000)
    parser.add_argument("--max-total", type=int, default=1000)
    parser.add_argument("--api-base", default=DEFAULT_API_BASE)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--poll-seconds", type=int, default=5)
    parser.add_argument("--max-polls", type=int, default=720)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    input_fasta = Path(args.input_fasta)
    aligned = input_fasta if args.aligned_input else outdir / "query_alignment.fasta"
    if not args.aligned_input:
        run_alignment(input_fasta, aligned, args.threads)

    domain_evalue = args.domain_evalue or args.evalue
    report_evalue = args.report_evalue or args.evalue
    report_domain_evalue = args.report_domain_evalue or domain_evalue
    alignment_text = aligned.read_text(encoding="utf-8", errors="replace")
    job_id = submit_hmmsearch(
        alignment_text,
        args.api_base,
        args.database,
        args.evalue,
        domain_evalue,
        report_evalue,
        report_domain_evalue,
        args.timeout,
    )
    (outdir / "hmmer_job_id.txt").write_text(job_id + "\n", encoding="utf-8")
    result = poll_result(args.api_base, job_id, args.timeout, args.poll_seconds, args.max_polls)
    (outdir / "hmmer_result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")

    hits = extract_hits(result)[:args.max_hits]
    write_hits_csv(hits, outdir / "hmmer_hits.csv")
    homologs = download_uniprot_fastas([hit.accession for hit in hits])
    seeds = read_fasta(input_fasta)
    homolog_only = dedupe_records(homologs, args.max_hits)
    combined = dedupe_records(seeds + homologs, args.max_total)
    write_fasta(seeds, outdir / "core.fasta")
    write_fasta(homolog_only, outdir / "homologs.fasta")
    write_fasta(combined, outdir / "homologs_plus_core.fasta")
    write_metadata_csv(fasta_metadata(combined), outdir / "homolog_metadata.csv")

    print(f"HMMER job id: {job_id}")
    print(f"Hits parsed: {len(hits)}")
    print(f"UniProt FASTA records downloaded: {len(homologs)}")
    print(f"Combined FASTA records: {len(combined)}")
    print(f"Output: {outdir}")


if __name__ == "__main__":
    main()
