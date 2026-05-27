#!/usr/bin/env python3
"""
Find tree UniProt IDs with likely enzyme literature that are missing from base.xlsx.

This is a fast triage layer: it queries UniProtKB JSON for tree nodes, extracts
annotation score, literature and PDB cross-references, and classifies references
as enzyme/structure/promiscuity-like versus genome/proteome-only.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}.json"
USER_AGENT = "hydrolase-homolog-tree/1.0"

GENOME_TERMS = [
    "genome", "genomic", "metagenome", "metagenomic", "chromosome",
    "transcriptome", "proteome", "proteomic", "whole-genome",
    "draft genome", "complete genome", "microbiome", "isolate genome",
]

ENZYME_TERMS = [
    "enzyme", "hydrolase", "lipase", "esterase", "dehalogenase",
    "epoxide hydrolase", "amidase", "amylase", "protease", "subtilisin",
    "catalytic", "catalysis", "activity", "substrate", "biochemical",
    "characterization", "characterisation", "cloning", "expression",
    "purification", "kinetic", "crystal", "structure", "active site",
    "mutant", "promiscu", "enantioselect", "transesterification",
    "aminolysis", "amidation", "michael", "mannich", "knoevenagel",
]


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
        except URLError as exc:
            if attempt == retries:
                print(f"  warning: URL error for {url}: {exc}")
                return None
        time.sleep(1.5 * (attempt + 1))
    return None


def citation_ids(citation: dict) -> Dict[str, str]:
    out = {"doi": "", "pubmed_id": "", "title": "", "journal": "", "year": ""}
    out["title"] = str(citation.get("title", "") or "")
    out["journal"] = str(citation.get("journal", "") or "")
    out["year"] = str(citation.get("publicationDate", "") or citation.get("year", "") or "")
    for cref in citation.get("citationCrossReferences", []) or []:
        db = str(cref.get("database", "")).lower()
        value = str(cref.get("id", "") or "")
        if db == "doi":
            out["doi"] = value
        elif db in {"pubmed", "pubmedid"}:
            out["pubmed_id"] = value
    return out


def term_hits(text: str, terms: List[str]) -> List[str]:
    low = text.lower()
    return [term for term in terms if term in low]


def annotation_score(payload: dict) -> str:
    score = payload.get("annotationScore")
    if isinstance(score, (int, float, str)):
        return str(score)
    if isinstance(score, dict):
        return str(score.get("score", "") or score.get("value", "") or "")
    return ""


def protein_name(payload: dict) -> str:
    protein = payload.get("proteinDescription", {}) or {}
    rec = protein.get("recommendedName", {}) or {}
    full = rec.get("fullName", {}) or {}
    return str(full.get("value", "") or "")


def organism_name(payload: dict) -> str:
    org = payload.get("organism", {}) or {}
    return str(org.get("scientificName", "") or org.get("commonName", "") or "")


def pdb_ids(payload: dict) -> List[str]:
    ids = []
    for xref in payload.get("uniProtKBCrossReferences", []) or []:
        if str(xref.get("database", "")).upper() == "PDB":
            ids.append(str(xref.get("id", "") or ""))
    return [x for x in ids if x]


def build_base_ids(base_rows: List[Dict[str, str]]) -> set:
    ids = set()
    for row in base_rows:
        uid = clean_accession(first_value(row, ["uniprot_id", "Uniprot id", "UniProt ID"]))
        if uid:
            ids.add(uid)
    return ids


def main() -> None:
    parser = argparse.ArgumentParser(description="Triage tree UniProt IDs for missed enzyme literature.")
    parser.add_argument("--nodes", required=True, help="ssn_nodes.csv from the tree/SSN run.")
    parser.add_argument("--base", required=True, help="seed_provenance_template_v2.csv or base-derived CSV.")
    parser.add_argument("--outdir", default="tree_literature_triage")
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--sleep", type=float, default=0.15)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--max", type=int, default=0, help="Debug limit; 0 means all.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    nodes = read_csv(Path(args.nodes))
    base_ids = build_base_ids(read_csv(Path(args.base)))

    node_by_uid: Dict[str, Dict[str, str]] = {}
    for row in nodes:
        uid = clean_accession(row.get("uniprot_id") or row.get("id"))
        if uid and not uid.upper().startswith("UPI"):
            node_by_uid.setdefault(uid, row)

    query_ids = [uid for uid in node_by_uid if uid not in base_ids]
    if args.max:
        query_ids = query_ids[:args.max]

    def process_uid(uid: str) -> tuple:
        payload = get_json(UNIPROT_ENTRY.format(accession=uid), args.timeout, args.retries)
        node = node_by_uid[uid]
        local_refs: List[Dict[str, str]] = []
        if not payload:
            return ({
                "uniprot_id": uid, "tree_id": node.get("id", ""), "status": "not_found",
                "review_priority": "skip", "reason": "UniProtKB not found",
            }, local_refs)
        enzyme_refs = []
        genome_refs = []
        other_refs = []
        for ridx, ref in enumerate(payload.get("references", []) or [], 1):
            citation = citation_ids(ref.get("citation", {}) or {})
            text = " ".join([citation["title"], citation["journal"], str(ref.get("referencePositions", ""))])
            e_hits = term_hits(text, ENZYME_TERMS)
            g_hits = term_hits(text, GENOME_TERMS)
            if e_hits:
                kind = "enzyme_or_structure"
                enzyme_refs.append(citation)
            elif g_hits:
                kind = "genome_or_proteome"
                genome_refs.append(citation)
            else:
                kind = "other_or_unclear"
                other_refs.append(citation)
            local_refs.append({
                "uniprot_id": uid,
                "tree_id": node.get("id", ""),
                "reference_rank": str(ridx),
                "reference_kind": kind,
                "enzyme_terms": ";".join(e_hits),
                "genome_terms": ";".join(g_hits),
                "doi": citation["doi"],
                "pubmed_id": citation["pubmed_id"],
                "title": citation["title"],
                "journal": citation["journal"],
                "year": citation["year"],
                "uniprot_url": f"https://www.uniprot.org/uniprotkb/{uid}/entry",
            })

        pdb = pdb_ids(payload)
        priority = "low"
        reason = "no enzyme-like publication title"
        if enzyme_refs:
            priority = "high"
            reason = "enzyme/structure/activity terms in UniProt references"
        elif pdb:
            priority = "medium"
            reason = "PDB cross-reference present; inspect primary citation"
        elif genome_refs and not enzyme_refs:
            priority = "skip"
            reason = "only genome/proteome-like references detected"

        best_ref = enzyme_refs[0] if enzyme_refs else other_refs[0] if other_refs else genome_refs[0] if genome_refs else {}
        summary = {
            "uniprot_id": uid,
            "tree_id": node.get("id", ""),
            "protein_name": protein_name(payload),
            "organism": organism_name(payload) or node.get("organism", ""),
            "kingdom": node.get("kingdom", ""),
            "annotation_score": annotation_score(payload),
            "pdb_ids": ";".join(pdb),
            "reference_count": str(len(payload.get("references", []) or [])),
            "enzyme_like_reference_count": str(len(enzyme_refs)),
            "genome_like_reference_count": str(len(genome_refs)),
            "review_priority": priority,
            "reason": reason,
            "best_doi": best_ref.get("doi", ""),
            "best_pubmed_id": best_ref.get("pubmed_id", ""),
            "best_title": best_ref.get("title", ""),
            "uniprot_url": f"https://www.uniprot.org/uniprotkb/{uid}/entry",
        }
        return summary, local_refs

    summary_rows: List[Dict[str, str]] = []
    ref_rows: List[Dict[str, str]] = []
    workers = max(1, args.workers)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_uid = {executor.submit(process_uid, uid): uid for uid in query_ids}
        for idx, future in enumerate(as_completed(future_to_uid), 1):
            uid = future_to_uid[future]
            try:
                summary, refs = future.result()
            except Exception as exc:
                node = node_by_uid[uid]
                summary, refs = ({
                    "uniprot_id": uid, "tree_id": node.get("id", ""), "status": "error",
                    "review_priority": "skip", "reason": f"{type(exc).__name__}: {exc}",
                }, [])
            summary_rows.append(summary)
            ref_rows.extend(refs)
            print(f"[{idx}/{len(query_ids)}] {uid} -> {summary.get('review_priority', '')}")
            if args.sleep:
                time.sleep(args.sleep)

    summary_headers = [
        "uniprot_id", "tree_id", "protein_name", "organism", "kingdom",
        "annotation_score", "pdb_ids", "reference_count",
        "enzyme_like_reference_count", "genome_like_reference_count",
        "review_priority", "reason", "best_doi", "best_pubmed_id",
        "best_title", "uniprot_url", "status",
    ]
    ref_headers = [
        "uniprot_id", "tree_id", "reference_rank", "reference_kind",
        "enzyme_terms", "genome_terms", "doi", "pubmed_id", "title",
        "journal", "year", "uniprot_url",
    ]
    write_csv(outdir / "tree_uniprot_literature_triage.csv", summary_rows, summary_headers)
    write_csv(outdir / "tree_uniprot_reference_details.csv", ref_rows, ref_headers)
    print(f"Wrote {outdir / 'tree_uniprot_literature_triage.csv'}")
    print(f"Wrote {outdir / 'tree_uniprot_reference_details.csv'}")


if __name__ == "__main__":
    main()
