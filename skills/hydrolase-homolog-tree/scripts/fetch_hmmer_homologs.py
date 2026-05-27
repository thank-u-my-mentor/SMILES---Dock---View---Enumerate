#!/usr/bin/env python3
"""
Fetch a bounded homolog FASTA set from a small curated seed FASTA.

Default mode uses the EMBL-EBI HMMER web service with phmmer for each seed
sequence, merges unique UniProt hits, and downloads FASTA entries from UniProt.
If --db-fasta is provided, local profile-HMM search is used instead:
seed FASTA -> ClustalO/MAFFT alignment -> hmmbuild -> hmmsearch.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import urljoin

HMMER_BASE = "https://www.ebi.ac.uk/Tools/hmmer"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
USER_AGENT = "hydrolase-homolog-tree/1.0"


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
    source_query: str
    evalue: Optional[float] = None
    score: Optional[float] = None
    description: str = ""


def read_fasta(path: Path) -> List[FastaRecord]:
    records: List[FastaRecord] = []
    header: Optional[str] = None
    seq: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header and seq:
                    records.append(FastaRecord(header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
    if header and seq:
        records.append(FastaRecord(header, "".join(seq)))
    return records


def write_fasta(records: Iterable[FastaRecord], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for rec in records:
            handle.write(f">{rec.header}\n")
            seq = rec.sequence.replace(" ", "").replace("\n", "")
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx:idx + 80] + "\n")


def run_command(cmd: Sequence[str], cwd: Optional[Path] = None) -> None:
    print(" ".join(str(x) for x in cmd))
    result = subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(cmd)}\n"
            f"STDOUT:\n{result.stdout[-2000:]}\nSTDERR:\n{result.stderr[-4000:]}"
        )


def first_tool(*names: str) -> Optional[str]:
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    return None


def clean_accession(raw: str) -> str:
    text = str(raw).strip()
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


def extract_hits_from_hmmer_json(payload: object, source_query: str) -> List[Hit]:
    hits: List[Hit] = []

    def walk(obj: object) -> None:
        if isinstance(obj, dict):
            keys = {str(k).lower(): k for k in obj.keys()}
            accession_key = next((keys[k] for k in ("acc", "accession", "target", "name") if k in keys), None)
            if accession_key is not None:
                accession = clean_accession(str(obj.get(accession_key, "")))
                if accession and accession.lower() not in {"none", "unknown"}:
                    hits.append(
                        Hit(
                            accession=accession,
                            source_query=source_query,
                            evalue=maybe_float(obj.get(keys.get("evalue", ""))),
                            score=maybe_float(obj.get(keys.get("score", "")) or obj.get(keys.get("bits", ""))),
                            description=str(obj.get(keys.get("desc", ""), "") or obj.get(keys.get("description", ""), "")),
                        )
                    )
            for value in obj.values():
                walk(value)
        elif isinstance(obj, list):
            for value in obj:
                walk(value)

    walk(payload)
    return hits


def hmmer_web_phmmer(query: FastaRecord, seqdb: str, max_hits: int, evalue: float, timeout: int) -> List[Hit]:
    import requests

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    endpoint = urljoin(HMMER_BASE + "/", "search/phmmer")
    fasta_text = f">{query.header}\n{query.sequence}\n"
    data = {
        "seq": fasta_text,
        "seqdb": seqdb,
        "E": str(evalue),
        "domE": str(evalue),
    }
    response = session.post(endpoint, data=data, allow_redirects=False, timeout=timeout)
    if response.status_code not in {200, 201, 202, 302, 303}:
        raise RuntimeError(f"HMMER phmmer failed for {query.id}: HTTP {response.status_code} {response.text[:300]}")

    result_url = response.headers.get("Location") or response.url
    if result_url.startswith("/"):
        result_url = urljoin(HMMER_BASE, result_url)

    last_text = ""
    for _ in range(90):
        params = {"output": "json", "range": f"1,{max_hits}"}
        result = session.get(result_url, params=params, timeout=timeout)
        last_text = result.text[:500]
        if result.status_code in {202, 204}:
            time.sleep(3)
            continue
        if result.status_code != 200:
            raise RuntimeError(f"HMMER result failed for {query.id}: HTTP {result.status_code} {last_text}")
        try:
            payload = result.json()
        except json.JSONDecodeError as exc:
            raise RuntimeError(f"HMMER result for {query.id} was not JSON: {last_text}") from exc
        return extract_hits_from_hmmer_json(payload, query.id)[:max_hits]
    raise TimeoutError(f"Timed out waiting for HMMER result for {query.id}: {last_text}")


def download_uniprot_fasta(accessions: Sequence[str], retries: int = 2) -> List[FastaRecord]:
    import requests

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    records: List[FastaRecord] = []
    for idx, accession in enumerate(accessions, 1):
        if idx % 50 == 0:
            print(f"Downloaded {len(records)}/{idx - 1} UniProt FASTA entries...")
        url = UNIPROT_ENTRY.format(accession=accession)
        for attempt in range(retries + 1):
            response = session.get(url, timeout=30)
            if response.status_code == 200 and response.text.startswith(">"):
                parsed = read_fasta_from_text(response.text)
                records.extend(parsed)
                break
            if response.status_code in {404, 400}:
                break
            time.sleep(1.5 * (attempt + 1))
    return records


def read_fasta_from_text(text: str) -> List[FastaRecord]:
    tmp_records: List[FastaRecord] = []
    header: Optional[str] = None
    seq: List[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header and seq:
                tmp_records.append(FastaRecord(header, "".join(seq)))
            header = line[1:].strip()
            seq = []
        else:
            seq.append(line)
    if header and seq:
        tmp_records.append(FastaRecord(header, "".join(seq)))
    return tmp_records


def build_seed_alignment(seed: Path, outdir: Path) -> Path:
    aln = outdir / "seed_aln.fasta"
    clustalo = first_tool("clustalo")
    mafft = first_tool("mafft")
    if clustalo:
        run_command([clustalo, "-i", str(seed), "-o", str(aln), "--outfmt", "fasta", "--force"])
    elif mafft:
        with aln.open("w", encoding="utf-8") as handle:
            result = subprocess.run([mafft, "--auto", str(seed)], text=True, capture_output=True)
            if result.returncode != 0:
                raise RuntimeError(result.stderr[-4000:])
            handle.write(result.stdout)
    else:
        raise RuntimeError("Need clustalo or mafft to build a seed alignment for local HMM search.")
    return aln


def parse_tblout(tblout: Path, max_hits: int) -> List[Hit]:
    hits: List[Hit] = []
    seen = set()
    with tblout.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            accession = parts[0]
            if accession in seen:
                continue
            seen.add(accession)
            evalue = maybe_float(parts[4])
            score = maybe_float(parts[5])
            desc = " ".join(parts[18:]) if len(parts) > 18 else ""
            hits.append(Hit(clean_accession(accession), "profile_hmm", evalue, score, desc))
            if len(hits) >= max_hits:
                break
    return hits


def fasta_index(records: Iterable[FastaRecord]) -> Dict[str, FastaRecord]:
    index: Dict[str, FastaRecord] = {}
    for rec in records:
        index[rec.id] = rec
        cleaned = clean_accession(rec.id)
        index.setdefault(cleaned, rec)
    return index


def local_hmmsearch(seed_fasta: Path, db_fasta: Path, outdir: Path, max_hits: int, evalue: float, threads: int) -> Tuple[List[Hit], List[FastaRecord]]:
    hmmbuild = first_tool("hmmbuild")
    hmmsearch = first_tool("hmmsearch")
    if not hmmbuild or not hmmsearch:
        raise RuntimeError("Need hmmbuild and hmmsearch on PATH for --mode local-hmmsearch.")
    aln = build_seed_alignment(seed_fasta, outdir)
    hmm = outdir / "seed_profile.hmm"
    tblout = outdir / "hmmsearch.tblout"
    run_command([hmmbuild, str(hmm), str(aln)])
    run_command([
        hmmsearch,
        "--cpu", str(threads),
        "-E", str(evalue),
        "--tblout", str(tblout),
        str(hmm),
        str(db_fasta),
    ])
    hits = parse_tblout(tblout, max_hits)
    db_records = fasta_index(read_fasta(db_fasta))
    selected = [db_records[h.accession] for h in hits if h.accession in db_records]
    return hits, selected


def deduplicate_records(seed_records: Sequence[FastaRecord], homolog_records: Sequence[FastaRecord], max_total: int) -> List[FastaRecord]:
    output: List[FastaRecord] = []
    seen = set()
    for rec in list(seed_records) + list(homolog_records):
        key = clean_accession(rec.id)
        if key in seen:
            continue
        seen.add(key)
        output.append(rec)
        if len(output) >= max_total:
            break
    return output


def write_hits_csv(hits: Sequence[Hit], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["accession", "source_query", "evalue", "score", "description"])
        for hit in hits:
            writer.writerow([hit.accession, hit.source_query, hit.evalue or "", hit.score or "", hit.description])


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch bounded HMMER homologs for hydrolase seed FASTA.")
    parser.add_argument("seed_fasta", help="Curated seed FASTA, usually 5-25 characterized enzymes.")
    parser.add_argument("--outdir", default="hmmer_homologs")
    parser.add_argument("--mode", choices=["auto", "web-phmmer", "local-hmmsearch"], default="auto")
    parser.add_argument("--db-fasta", default=None, help="Local protein FASTA database for profile hmmsearch.")
    parser.add_argument("--seqdb", default="uniprotkb", help="HMMER web sequence database.")
    parser.add_argument("--max-hits", type=int, default=2500, help="Homolog cap before adding seeds.")
    parser.add_argument("--max-total", type=int, default=3000, help="Final FASTA cap including seeds.")
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=60)
    args = parser.parse_args()

    seed_fasta = Path(args.seed_fasta).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    seed_records = read_fasta(seed_fasta)
    if not 1 <= len(seed_records) <= 100:
        raise ValueError(f"Unexpected seed count {len(seed_records)}. This workflow expects a small curated FASTA.")

    core_fasta = outdir / "core.fasta"
    write_fasta(seed_records, core_fasta)

    mode = args.mode
    if mode == "auto":
        mode = "local-hmmsearch" if args.db_fasta else "web-phmmer"

    all_hits: List[Hit] = []
    homolog_records: List[FastaRecord] = []
    if mode == "local-hmmsearch":
        if not args.db_fasta:
            raise ValueError("--db-fasta is required for local-hmmsearch.")
        all_hits, homolog_records = local_hmmsearch(
            seed_fasta=seed_fasta,
            db_fasta=Path(args.db_fasta).resolve(),
            outdir=outdir,
            max_hits=args.max_hits,
            evalue=args.evalue,
            threads=args.threads,
        )
    else:
        per_query = max(50, args.max_hits // max(1, len(seed_records)) + 25)
        hit_by_acc: Dict[str, Hit] = {}
        for query in seed_records:
            print(f"HMMER web phmmer: {query.id}")
            hits = hmmer_web_phmmer(query, args.seqdb, per_query, args.evalue, args.timeout)
            for hit in hits:
                existing = hit_by_acc.get(hit.accession)
                if existing is None or ((hit.evalue or 1e99) < (existing.evalue or 1e99)):
                    hit_by_acc[hit.accession] = hit
            if len(hit_by_acc) >= args.max_hits:
                break
        all_hits = sorted(hit_by_acc.values(), key=lambda h: (h.evalue if h.evalue is not None else 1e99, -(h.score or 0)))[:args.max_hits]
        homolog_records = download_uniprot_fasta([h.accession for h in all_hits])

    write_hits_csv(all_hits, outdir / "hmmer_hits.csv")
    combined = deduplicate_records(seed_records, homolog_records, args.max_total)
    homolog_only = deduplicate_records([], homolog_records, args.max_hits)
    write_fasta(homolog_only, outdir / "homologs.fasta")
    write_fasta(combined, outdir / "homologs_plus_core.fasta")
    print(f"Seed records: {len(seed_records)}")
    print(f"HMMER hits: {len(all_hits)}")
    print(f"Downloaded homolog FASTA records: {len(homolog_records)}")
    print(f"Combined FASTA records: {len(combined)}")
    print(f"Output: {outdir}")


if __name__ == "__main__":
    main()
