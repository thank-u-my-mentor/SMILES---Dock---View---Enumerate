#!/usr/bin/env python3
"""
Enrich seed provenance rows with UniProt literature/PDB candidates.

This is deliberately a first-pass discovery layer, not a final evidence parser:
it uses UniProt accessions from a seed CSV, retrieves UniProtKB JSON, extracts
publication/citation IDs and PDB cross-references, and writes reviewable CSVs.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPARC_ENTRY = "https://rest.uniprot.org/uniparc/{accession}.json"
USER_AGENT = "hydrolase-homolog-tree/1.0"


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: List[Dict[str, str]], headers: List[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_accession(raw: str) -> str:
    text = str(raw or "").strip()
    if "|" in text:
        for part in text.split("|"):
            if re.match(r"^[A-Z0-9]{6,10}(?:-\d+)?$", part.strip()):
                return part.strip()
    text = re.sub(r"^UniRef\d+_", "", text)
    text = re.sub(r"\.\d+$", "", text)
    match = re.search(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})", text)
    return match.group(1) if match else text.split()[0] if text else ""


def first_value(row: Dict[str, str], names: Iterable[str]) -> str:
    lower = {k.lower(): k for k in row}
    for name in names:
        key = lower.get(name.lower())
        if key and str(row.get(key, "")).strip():
            return str(row[key]).strip()
    return ""


def get_json(url: str, timeout: int, retries: int) -> Optional[dict]:
    for attempt in range(retries + 1):
        req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
        try:
            with urlopen(req, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except HTTPError as exc:
            if exc.code == 404:
                return None
            if attempt == retries:
                print(f"  warning: HTTP {exc.code} for {url}")
                return None
        except Exception as exc:
            if attempt == retries:
                print(f"  warning: {type(exc).__name__} for {url}: {exc}")
                return None
        except URLError:
            if attempt == retries:
                raise
        time.sleep(1.5 * (attempt + 1))
    return None


def citation_ids(citation: dict) -> Dict[str, str]:
    ids = {"doi": "", "pubmed_id": "", "title": "", "journal": "", "year": ""}
    ids["title"] = str(citation.get("title", "") or "")
    ids["journal"] = str(citation.get("journal", "") or citation.get("citationCrossReferences", "") or "")
    ids["year"] = str(citation.get("publicationDate", "") or citation.get("year", "") or "")
    for cref in citation.get("citationCrossReferences", []) or []:
        db = str(cref.get("database", "")).lower()
        value = str(cref.get("id", "") or "")
        if db == "doi":
            ids["doi"] = value
        elif db in {"pubmed", "pubmedid"}:
            ids["pubmed_id"] = value
    return ids


def extract_publications(payload: dict, seed: Dict[str, str], accession: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    refs = payload.get("references", []) or []
    for idx, ref in enumerate(refs, 1):
        citation = ref.get("citation", {}) or {}
        ids = citation_ids(citation)
        scopes = []
        for scope in ref.get("referencePositions", []) or []:
            scopes.append(str(scope))
        rows.append({
            "seed_uniprot_id": accession,
            "seed_enzyme_name": first_value(seed, ["enzyme_name", "enzyme"]),
            "seed_source_organism": first_value(seed, ["source_organism", "organism"]),
            "source": "UniProt reference",
            "candidate_rank": str(idx),
            "candidate_doi": ids["doi"],
            "candidate_pubmed_id": ids["pubmed_id"],
            "candidate_title": ids["title"],
            "candidate_journal": ids["journal"],
            "candidate_year": ids["year"],
            "candidate_note": "; ".join(scopes),
            "candidate_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
        })
    return rows


def extract_pdb(payload: dict, seed: Dict[str, str], accession: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    xrefs = payload.get("uniProtKBCrossReferences", []) or []
    rank = 0
    for xref in xrefs:
        if str(xref.get("database", "")).upper() != "PDB":
            continue
        rank += 1
        pdb_id = str(xref.get("id", "") or "")
        props = {p.get("key", ""): p.get("value", "") for p in xref.get("properties", []) or []}
        rows.append({
            "seed_uniprot_id": accession,
            "seed_enzyme_name": first_value(seed, ["enzyme_name", "enzyme"]),
            "seed_source_organism": first_value(seed, ["source_organism", "organism"]),
            "source": "UniProt PDB cross-reference",
            "candidate_rank": str(rank),
            "pdb_id": pdb_id,
            "method": props.get("Method", ""),
            "resolution": props.get("Resolution", ""),
            "chains": props.get("Chains", ""),
            "candidate_url": f"https://www.rcsb.org/structure/{pdb_id}" if pdb_id else "",
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Discover seed literature/PDB candidates from UniProt accessions.")
    parser.add_argument("seed_csv", help="CSV from extract_seed_metadata_from_xlsx.py.")
    parser.add_argument("--outdir", default="seed_literature_uniprot")
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--sleep", type=float, default=0.2)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    seeds = read_csv(Path(args.seed_csv))

    publication_rows: List[Dict[str, str]] = []
    pdb_rows: List[Dict[str, str]] = []
    enriched_rows: List[Dict[str, str]] = []
    for idx, seed in enumerate(seeds, 1):
        accession = clean_accession(first_value(seed, ["uniprot_id", "Uniprot id", "UniProt ID"]))
        out = dict(seed)
        out.setdefault("uniprot_literature_status", "")
        out.setdefault("uniprot_candidate_doi", "")
        out.setdefault("uniprot_candidate_pubmed_id", "")
        out.setdefault("uniprot_candidate_pdb", "")
        if not accession:
            out["uniprot_literature_status"] = "missing_uniprot_id"
            enriched_rows.append(out)
            continue
        print(f"[{idx}/{len(seeds)}] UniProt {accession}")
        if accession.upper().startswith("UPI"):
            payload = get_json(UNIPARC_ENTRY.format(accession=accession), timeout=args.timeout, retries=args.retries)
            out["uniprot_literature_status"] = "uniparc_checked"
            enriched_rows.append(out)
            time.sleep(args.sleep)
            continue
        payload = get_json(UNIPROT_ENTRY.format(accession=accession), timeout=args.timeout, retries=args.retries)
        if not payload:
            out["uniprot_literature_status"] = "not_found"
            enriched_rows.append(out)
            continue
        pubs = extract_publications(payload, seed, accession)
        pdbs = extract_pdb(payload, seed, accession)
        publication_rows.extend(pubs)
        pdb_rows.extend(pdbs)
        first_doi = next((r["candidate_doi"] for r in pubs if r.get("candidate_doi")), "")
        first_pm = next((r["candidate_pubmed_id"] for r in pubs if r.get("candidate_pubmed_id")), "")
        first_pdb = next((r["pdb_id"] for r in pdbs if r.get("pdb_id")), "")
        if not first_value(seed, ["key_doi", "doi", "DOI"]) and first_doi:
            out["key_doi"] = first_doi
            out["evidence_source"] = f"UniProt {accession}"
            out["provenance_confidence"] = out.get("provenance_confidence") or "low"
            out["notes_for_cloning"] = (out.get("notes_for_cloning", "") + " DOI candidate filled from UniProt; verify paper methods.").strip()
        out["uniprot_literature_status"] = f"references={len(pubs)};pdb={len(pdbs)}"
        out["uniprot_candidate_doi"] = first_doi
        out["uniprot_candidate_pubmed_id"] = first_pm
        out["uniprot_candidate_pdb"] = first_pdb
        enriched_rows.append(out)
        time.sleep(args.sleep)

    seed_headers = list(enriched_rows[0].keys()) if enriched_rows else []
    for col in ["uniprot_literature_status", "uniprot_candidate_doi", "uniprot_candidate_pubmed_id", "uniprot_candidate_pdb"]:
        if col not in seed_headers:
            seed_headers.append(col)
    write_csv(outdir / "seed_provenance_uniprot_enriched.csv", enriched_rows, seed_headers)

    pub_headers = [
        "seed_uniprot_id", "seed_enzyme_name", "seed_source_organism", "source",
        "candidate_rank", "candidate_doi", "candidate_pubmed_id", "candidate_title",
        "candidate_journal", "candidate_year", "candidate_note", "candidate_url",
    ]
    pdb_headers = [
        "seed_uniprot_id", "seed_enzyme_name", "seed_source_organism", "source",
        "candidate_rank", "pdb_id", "method", "resolution", "chains", "candidate_url",
    ]
    write_csv(outdir / "uniprot_publication_candidates.csv", publication_rows, pub_headers)
    write_csv(outdir / "uniprot_pdb_candidates.csv", pdb_rows, pdb_headers)
    print(f"Wrote {outdir}")
    print(f"Publication candidates: {len(publication_rows)}")
    print(f"PDB candidates: {len(pdb_rows)}")


if __name__ == "__main__":
    main()
