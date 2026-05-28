#!/usr/bin/env python3
"""Build a human-readable sequence seed review table from literature-agent CSVs.

Default mode is offline: infer organism/source clues and concrete enzyme IDs from
existing CSV evidence, then write an HTML review dashboard plus an enhanced CSV.
Network sequence fetching is intentionally not enabled here yet; use the emitted
database query links for manual review or add explicit fetch modes later.
"""

from __future__ import annotations

import argparse
import csv
import html
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import quote_plus


CONTEXT_FIELDS = [
    "enzyme_name",
    "enzyme_name_or_target",
    "title",
    "abstract",
    "evidence_summary",
    "metadata_keywords",
    "enzyme_family",
    "flavin_cofactor",
    "catalyst_or_photosensitizer",
    "cofactor_or_photosensitizer",
    "mechanistic_evidence",
    "structure_similarity_hint",
    "search_track",
    "paper_domain",
    "reaction_type",
]

GENERIC_NAME_PATTERNS = [
    r"\benzymes?\b",
    r"\bfamily members?\b",
    r"\bpanel\b",
    r"\bspecific variants? not named\b",
    r"\bspecific enzyme unclear\b",
    r"\bspecific not named\b",
    r"\brounds? of protein engineering\b",
    r"\bprotein engineering\b",
    r"\bFMNsq\b",
    r"\bflavin semiquinone\b",
    r"\bTyr\s*[->→]\s*Phe\b",
    r"\bCys-incorporated\b",
    r"\bvariants?\b",
    r"\bmutants?\b",
    r"\bwild type\b",
    r"\bapo-",
    r"\bholo-",
    r"\bbiosynthesis enzyme\b",
    r"\bflavoenzymes?\b",
    r"\boxidases? and dehydrogenases?\b",
    r"\bdehydrogenases?\b",
    r"\breductases?\b",
]

REVIEW_HINTS = [
    "review", "strategies", "advances", "versatile catalysts", "biochemical implications",
    "accounts", "perspective", "concept", "on their way",
]

ORGANISM_ALIASES = {
    "CvFAP": "Chlorella variabilis",
    "Chlorella variabilis FAP": "Chlorella variabilis",
    "GluER": "Gluconobacter oxydans",
    "GsOYE": "Galdieria sulphuraria",
    "PaDADH": "Pseudomonas aeruginosa",
    "PfBAL": "Pseudomonas fluorescens",
    "Pseudomonas aeruginosa D-arginine dehydrogenase": "Pseudomonas aeruginosa",
    "avenolide biosynthetic flavoenzyme": "Streptomyces avermitilis",
}

KNOWN_ORGANISM_NAMES = [
    "Bacillus amyloliquefaciens",
    "Enterobacter cloacae",
    "Galdieria sulphuraria",
    "Gluconobacter oxydans",
    "Chlorella variabilis",
    "Gluconobacter",
    "Saccharomyces pastorianus",
    "Saccharomyces cerevisiae",
    "Nicotiana tabacum",
    "Escherichia coli",
    "Pseudomonas aeruginosa",
    "Pseudomonas fluorescens",
    "Streptomyces avermitilis",
]

NAME_PATTERNS = [
    r"\bPfBAL(?:-[A-Z0-9]+)*\b",
    r"\bbenzaldehyde lyase\b",
    r"\b[A-Z][a-z]{1,4}ER(?:-[A-Z0-9]+)*\b",
    r"\b[A-Z][a-z]OYE\d*(?:-[A-Z0-9]+)*\b",
    r"\bOYE\d+(?:-[A-Z0-9]+)*\b",
    r"\bCvFAP(?:-[A-Z][0-9]+[A-Z])*\b",
    r"\bFAP(?:-[A-Z][0-9]+[A-Z])*\b",
    r"\b[A-Z][a-z]{1,5}FDH(?:-[A-Z0-9]+)*\b",
    r"\b[A-Z][a-z]{1,5}NTR(?:-[A-Z0-9]+)*\b",
    r"\bPaDADH(?:-[A-Z0-9]+)*\b",
    r"\bPseudomonas aeruginosa D-arginine dehydrogenase\b",
    r"\bD-2-hydroxyglutarate dehydrogenase\b",
    r"\bovenolide biosynthetic flavoenzyme\b",
    r"\b[A-Z][A-Za-z0-9]{1,8}(?:ase|ER|OYE|FAP|FDH|NTR)(?:-[A-Z0-9]+)*\b",
]

MUTATION_RE = re.compile(r"\b[A-Z][0-9]{1,4}[A-Z]\b")
EC_RE = re.compile(r"\bEC\s*(?:number|no\.?)?\s*[:#]?\s*(\d+\.\d+\.\d+\.(?:\d+|-))\b", re.I)
FALSE_PDB_IDS = {"3CL2", "6H2O"}

KNOWN_PAPER_OVERRIDES = {
    "10.1038/s41586-023-06822-x": {
        "resolution_status": "non_flavin_exclude",
        "resolution_reason": "ThDP-dependent PfBAL paper; exclude from flavin-dependent FASTA seeds.",
        "specific_enzyme_names": "PfBAL | benzaldehyde lyase",
        "parent_enzyme_names": "benzaldehyde lyase",
        "canonical_enzyme": "PfBAL",
        "mutations_or_variants": "",
        "organism_candidates": "Pseudomonas fluorescens",
        "organism_evidence": "doi_override:benzaldehyde lyase from Pseudomonas fluorescens",
        "enzyme_family": "ThDP-dependent lyase",
        "enzyme_name_original": "benzaldehyde lyase from Pseudomonas fluorescens (PfBAL)",
        "flavin_cofactor": "no flavin enzyme cofactor; ThDP enzyme plus eosin Y photocatalyst",
        "cofactor_class": "thdp",
        "cofactor_detail": "ThDP | eosin Y external photocatalyst",
        "cofactor_evidence": "doi_override:user correction; article is ThDP-dependent radical acylation",
        "flavin_dependency_status": "non_flavin_thdp",
        "enzyme_function_class": "ThDP-dependent benzaldehyde lyase",
        "ec_number_candidates": "",
        "sequence_resolution_route": "exclude_non_flavin",
        "pdb_ids": "",
        "pdb_evidence": "",
    },
    "10.1038/s41929-023-01065-5": {
        "resolution_status": "family_only",
        "resolution_reason": "Uses flavin-dependent EREDs, but local metadata does not resolve a concrete enzyme/organism; do not infer OYE1.",
        "specific_enzyme_names": "",
        "parent_enzyme_names": "ERED/OYE family",
        "canonical_enzyme": "ERED/OYE family",
        "mutations_or_variants": "",
        "organism_candidates": "",
        "organism_evidence": "doi_override:ERED family, no specific organism in local metadata",
        "enzyme_family": "ERED/OYE",
        "enzyme_name_original": "flavin-dependent ene-reductases (EREDs)",
        "flavin_cofactor": "flavin unspecified; exogenous Ru(bpy)3 photosensitizer",
        "cofactor_class": "flavin",
        "cofactor_detail": "flavin | Ru(bpy)3 external photocatalyst",
        "cofactor_evidence": "doi_override:flavin-dependent ene-reductases",
        "flavin_dependency_status": "likely_flavin_dependent",
        "enzyme_function_class": "ERED/OYE ene-reductase",
        "ec_number_candidates": "",
        "sequence_resolution_route": "paper_supplement_then_database",
    },
    "10.1002/anie.202311762": {
        "cofactor_class": "flavin",
        "cofactor_detail": "FMN/flavin in GluER",
        "cofactor_evidence": "abstract:flavin-dependent ene-reductase GluER; PDB 6O08 contains FMN",
        "flavin_dependency_status": "confirmed_flavin_dependent",
        "enzyme_function_class": "ERED/OYE ene-reductase",
        "sequence_resolution_route": "pdb_first",
        "pdb_ids": "6O08",
        "uniprot_id": "A1E8I9",
        "organism_candidates": "Gluconobacter oxydans",
        "organism_evidence": "doi_override:PDB 6O08 organism Gluconobacter oxydans",
    },
    "10.1002/ange.202311762": {
        "cofactor_class": "flavin",
        "cofactor_detail": "FMN/flavin in GluER",
        "cofactor_evidence": "abstract:flavin-dependent ene-reductase GluER; PDB 6O08 contains FMN",
        "flavin_dependency_status": "confirmed_flavin_dependent",
        "enzyme_function_class": "ERED/OYE ene-reductase",
        "sequence_resolution_route": "pdb_first",
        "pdb_ids": "6O08",
        "uniprot_id": "A1E8I9",
        "organism_candidates": "Gluconobacter oxydans",
        "organism_evidence": "doi_override:PDB 6O08 organism Gluconobacter oxydans",
    },
    "10.1021/jacs.5c19848": {
        "resolution_status": "specific_with_organism",
        "resolution_reason": "Supporting information gives wild-type GsOYE accession M2XAQ9 and PDB 6S0G.",
        "specific_enzyme_names": "GsOYE | wild-type GsOYE",
        "parent_enzyme_names": "GsOYE",
        "canonical_enzyme": "GsOYE",
        "mutations_or_variants": "wild type",
        "organism_candidates": "Galdieria sulphuraria",
        "organism_evidence": "supporting_information:OYE from Galdieria sulphuraria (GsOYE)",
        "enzyme_family": "ERED/OYE",
        "enzyme_name_original": "OYE from Galdieria sulphuraria (GsOYE), wild type",
        "cofactor_class": "flavin",
        "cofactor_detail": "FMN/flavin",
        "cofactor_evidence": "supporting_information:flavin-dependent ERED/OYE GsOYE",
        "flavin_dependency_status": "confirmed_flavin_dependent",
        "enzyme_function_class": "ERED/OYE ene-reductase repurposed as isomerase",
        "ec_number_candidates": "",
        "sequence_resolution_route": "supplement_accession_first",
        "pdb_ids": "6S0G",
        "pdb_evidence": "supporting_information:PDB 6S0G; wild-type GsOYE",
        "uniprot_id": "M2XAQ9",
    },
}


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: Sequence[Dict[str, str]], fields: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def norm_space(value: str) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def row_context(row: Dict[str, str]) -> str:
    return " ".join(norm_space(row.get(k, "")) for k in CONTEXT_FIELDS)


def split_multi(value: str) -> List[str]:
    parts = re.split(r"\s*\|\s*|[;,]\s*", str(value or ""))
    return [norm_space(p) for p in parts if norm_space(p)]


def clean_candidate_name(value: str) -> str:
    value = norm_space(value)
    value = re.sub(r"^\(?\s*(engineered|evolved|directed evolution|mutant|variant)\s+", "", value, flags=re.I)
    value = re.sub(r"\s+from\s+[A-Z][a-z]+(?:\s+[a-z][a-z-]+){0,3}\b", "", value)
    value = re.sub(r"\bspecific\s+variants?\s+not\s+named(?:\s+in\s+abstract)?\b", "specific not named", value, flags=re.I)
    value = re.sub(r"\b(?:mutants?|variants?)\b", "", value, flags=re.I)
    value = re.sub(r"\b[A-Z][0-9]{1,4}[A-Z]\s+mutant\b", "", value)
    value = re.sub(r"\bTyr\s*[->→]\s*Phe/Trp\s+mutants?\b", "", value, flags=re.I)
    value = value.strip(" .,:;()[]")
    value = norm_space(value)
    return value


def is_generic_name(value: str, title: str = "") -> bool:
    low = f"{value} {title}".lower()
    if not value or value.lower() in {"unclear", "nan", "none", "n/a", "specific not named", "specific not named in abstract"}:
        return True
    if re.fullmatch(r"EC\s*\d+\.\d+\.\d+\.(?:\d+|-)", value, flags=re.I):
        return True
    if re.fullmatch(r"UniProt\s+[A-Z0-9]{4,12}", value, flags=re.I):
        return True
    if any(phrase in value.lower() for phrase in [
        "rounds of protein engineering", "specific not named", "specific enzyme unclear",
        "flavin semiquinone", "fmnsq", "cys-incorporated", "tyr→phe", "tyr->phe",
    ]):
        return True
    if re.fullmatch(r"[A-Z][0-9]{1,4}[A-Z](?:\s+mutant)?", value):
        return True
    if re.fullmatch(r"(?:specific\s+)?variants?(?:\s+not\s+named(?:\s+in\s+abstract)?)?", value, flags=re.I):
        return True
    if re.fullmatch(r"(?:mutants?|wild type|apo-[A-Za-z0-9]+|holo-[A-Za-z0-9]+)", value, flags=re.I):
        return True
    if any(re.search(pattern, value.lower()) for pattern in GENERIC_NAME_PATTERNS):
        return True
    if len(value) > 120 or value.count("|") >= 2:
        return True
    if any(hint in low for hint in REVIEW_HINTS):
        return True
    return False


def extract_specific_names(row: Dict[str, str]) -> List[str]:
    text = row_context(row)
    names: List[str] = []
    for field in ["enzyme_name", "enzyme_name_or_target"]:
        for part in split_multi(row.get(field, "")):
            cleaned = clean_candidate_name(part)
            if cleaned and not is_generic_name(cleaned):
                names.append(cleaned)
    for pattern in NAME_PATTERNS:
        for match in re.findall(pattern, text):
            cleaned = clean_candidate_name(match)
            if cleaned and not is_generic_name(cleaned):
                names.append(cleaned)
    return prefer_specific_names(dedupe(names))[:10]


def prefer_specific_names(names: Sequence[str]) -> List[str]:
    scored = []
    for name in names:
        low = name.lower()
        score = 0
        if re.search(r"\b(?:[A-Z][a-z]{1,4}ER|[A-Z][a-z]OYE\d*|OYE\d+|CvFAP|GsOYE|OaER)\b", name):
            score -= 4
        if any(alias.lower() in low for alias in ORGANISM_ALIASES):
            score -= 2
        if " from " in low:
            score -= 1
        if any(generic in low for generic in ["fatty acid photodecarboxylase", "photodecarboxylase", "ered/oye"]):
            score += 2
        if len(name) > 70:
            score += 3
        scored.append((score, name))
    return [name for _, name in sorted(scored, key=lambda item: (item[0], item[1].lower()))]


def parent_enzyme_name(name: str) -> str:
    value = norm_space(name)
    value = re.sub(r"\s*\([^)]*\)", "", value)
    value = re.sub(r"\s+from\s+[A-Z][a-z]+(?:\s+[a-z][a-z-]+){0,3}\b", "", value)
    value = re.sub(r"(?:-[A-Z][0-9]{1,4}[A-Z])+$", "", value)
    return value.strip(" -")


def extract_mutations(names: Sequence[str], row: Dict[str, str]) -> List[str]:
    text = " ".join(list(names) + [row.get("enzyme_name", ""), row.get("enzyme_name_or_target", ""), row.get("abstract", ""), row.get("evidence_summary", "")])
    return dedupe(MUTATION_RE.findall(text))


def extract_organisms(row: Dict[str, str], names: Sequence[str]) -> Tuple[List[str], List[str]]:
    text = row_context(row)
    organisms: List[str] = []
    evidence: List[str] = []
    for name in names:
        for key, organism in ORGANISM_ALIASES.items():
            if key.lower() in name.lower():
                organisms.append(organism)
                evidence.append(f"alias:{key}")
    for organism in KNOWN_ORGANISM_NAMES:
        if organism.lower() in text.lower():
            organisms.append(organism)
            evidence.append(f"text:{organism}")
    for match in re.findall(r"\bfrom\s+([A-Z][a-z]+(?:\s+[a-z][a-z-]+){0,2})", text):
        if not looks_like_false_organism(match):
            organisms.append(match)
            evidence.append(f"from:{match}")
    for match in re.findall(r"\b([A-Z][a-z]+(?:\s+[a-z][a-z-]+){0,2})\s+(?:FAP|OYE|ERED|ene-reductase|nitroreductase|monooxygenase|halogenase)\b", text):
        if not looks_like_false_organism(match):
            organisms.append(match)
            evidence.append(f"prefix:{match}")
    return dedupe(organisms), dedupe(evidence)


def looks_like_false_organism(value: str) -> bool:
    low = value.lower().strip()
    bad = {
        "the", "this", "these", "photoenzymatic", "flavin", "fatty acid", "old yellow",
        "single electron", "directed evolution", "natural redox", "functional group",
        "journal of", "nature communications",
        "enzyme", "mass spectra of", "photodecarboxylase", "fatty acid photodecarboxylase",
        "th", "abstract photoenzymatic", "light footprinting", "collision cross",
        "nitronate", "enterobacter cloacae enables", "saccharomyces cerevisiae and",
    }
    if len(low) < 5:
        return True
    if low in bad:
        return True
    if any(low.startswith(prefix) for prefix in ["the ", "this ", "these ", "mass ", "abstract "]):
        return True
    if any(word in low for word in ["spectra", "reaction", "substrate", "protein engineering"]):
        return True
    return False


def dedupe(values: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for value in values:
        value = norm_space(value)
        key = value.lower()
        if not value or key in seen:
            continue
        seen.add(key)
        out.append(value)
    return out


def extract_pdb_ids(row: Dict[str, str]) -> Tuple[List[str], List[str]]:
    text = " ".join([row.get("pdb_ids", ""), row_context(row)])
    ids: List[str] = []
    evidence: List[str] = []
    for raw in split_multi(row.get("pdb_ids", "")):
        candidate = raw.upper()
        if is_pdb_id(candidate):
            ids.append(candidate)
            evidence.append(f"csv:{candidate}")
    for match in re.findall(r"\b(?:PDB|RCSB|structure)\s*(?:ID|entry|code|accession)?s?\s*[:#]?\s*([0-9][A-Za-z0-9]{3})\b", text, flags=re.I):
        candidate = match.upper()
        if is_pdb_id(candidate):
            ids.append(candidate)
            evidence.append(f"text:{candidate}")
    return dedupe(ids), dedupe(evidence)


def is_pdb_id(value: str) -> bool:
    value = norm_space(value).upper()
    return bool(re.fullmatch(r"[0-9][A-Z0-9]{3}", value)) and value not in FALSE_PDB_IDS


def pdb_urls(pdb_ids: Sequence[str], query: str) -> Dict[str, str]:
    pdb_ids = [p.upper() for p in pdb_ids if is_pdb_id(p)]
    if pdb_ids:
        pdb_entry_url = " | ".join(f"https://www.rcsb.org/structure/{p}" for p in pdb_ids)
        rcsb_query = "https://www.rcsb.org/search?request=" + quote_plus(" ".join(pdb_ids))
    else:
        pdb_entry_url = ""
        rcsb_query = "https://www.rcsb.org/search?request=" + quote_plus(query)
    return {
        "pdb_url": pdb_entry_url,
        "rcsb_search_url": rcsb_query,
    }


def extract_ec_numbers(row: Dict[str, str]) -> List[str]:
    return dedupe(EC_RE.findall(row_context(row)))


def infer_cofactor_and_function(row: Dict[str, str], names: Sequence[str]) -> Dict[str, str]:
    text = " ".join([row_context(row), " ".join(names)]).lower()
    evidence: List[str] = []
    detail: List[str] = []

    flavin_hits = []
    for label, pattern in [
        ("FAD", r"\bfad\b|flavin adenine dinucleotide"),
        ("FMN", r"\bfmn\b|flavin mononucleotide"),
        ("flavin", r"\bflavin\b|\bflavoenzyme\b|flavin-dependent|flavin dependent"),
    ]:
        if re.search(pattern, text):
            flavin_hits.append(label)
    if flavin_hits:
        detail.extend(dedupe(flavin_hits))
        evidence.append("text:" + "/".join(dedupe(flavin_hits)))

    non_flavin_hits = []
    for label, pattern in [
        ("ThDP", r"\bthdp\b|thiamine diphosphate|thiamin diphosphate|thiamine pyrophosphate|\btpp\b"),
        ("PQQ", r"\bpqq\b|pyrroloquinoline quinone"),
        ("PLP", r"\bplp\b|pyridoxal phosphate|pyridoxal-5"),
        ("heme", r"\bheme\b|\bhaem\b|\bcytochrome p450\b|\bp450\b"),
        ("iron", r"nonheme iron|non-heme iron|\biron enzyme\b"),
        ("NAD(P)H", r"\bnadph\b|\bnadh\b|\bnad\(p\)h\b|\bnadp\+\b|\bnad\+\b"),
    ]:
        if re.search(pattern, text):
            non_flavin_hits.append(label)
    if non_flavin_hits:
        detail.extend(dedupe(non_flavin_hits))
        evidence.append("text:" + "/".join(dedupe(non_flavin_hits)))

    family, function = infer_enzyme_function(row, names)
    cofactor_class = "unknown"
    status = "unknown"

    if "ThDP" in non_flavin_hits and not flavin_hits:
        cofactor_class = "thdp"
        status = "non_flavin_thdp"
    elif "PQQ" in non_flavin_hits and not flavin_hits:
        cofactor_class = "pqq"
        status = "non_flavin_pqq"
    elif "PLP" in non_flavin_hits and not flavin_hits:
        cofactor_class = "plp"
        status = "non_flavin_plp"
    elif "heme" in non_flavin_hits and not flavin_hits:
        cofactor_class = "heme"
        status = "non_flavin_heme"
    elif flavin_hits and any(hit in non_flavin_hits for hit in ["ThDP", "PQQ", "PLP", "heme", "iron"]):
        cofactor_class = "mixed_or_external"
        status = "mixed_or_review"
    elif flavin_hits:
        cofactor_class = "flavin"
        if re.search(r"flavin-dependent|flavin dependent|flavoenzyme|fmn-dependent|fad-dependent", text):
            status = "confirmed_flavin_dependent"
        else:
            status = "likely_flavin_dependent"
    elif family in {"ERED/OYE ene-reductase", "fatty acid photodecarboxylase", "flavin-dependent halogenase", "flavin monooxygenase", "nitroreductase", "photolyase/cryptochrome"}:
        cofactor_class = "flavin"
        status = "likely_flavin_dependent"
        evidence.append(f"family:{family}")
    elif "NAD(P)H" in non_flavin_hits:
        cofactor_class = "nad(p)h_only"
        status = "unknown_or_nadph_only"

    return {
        "cofactor_class": cofactor_class,
        "cofactor_detail": " | ".join(dedupe(detail)),
        "cofactor_evidence": " | ".join(dedupe(evidence)),
        "flavin_dependency_status": status,
        "enzyme_function_class": function,
    }


def infer_enzyme_function(row: Dict[str, str], names: Sequence[str]) -> Tuple[str, str]:
    text = " ".join([row_context(row), " ".join(names)]).lower()
    checks = [
        ("ThDP-dependent benzaldehyde lyase", r"pfbal|benzaldehyde lyase|thdp|thiamine diphosphate"),
        ("ERED/OYE ene-reductase", r"\bered\b|ene[- ]reductase|old yellow enzyme|\boye\d*\b|gluer|gsoye|oaer"),
        ("fatty acid photodecarboxylase", r"\bfap\b|fatty acid photodecarboxylase|photodecarboxylase|cvfap"),
        ("flavin-dependent halogenase", r"flavin-dependent halogenase|flavin dependent halogenase|halogenase"),
        ("flavin monooxygenase", r"flavin monooxygenase|\bfmo\b|monooxygenase"),
        ("nitroreductase", r"nitroreductase|\bntr\b|nitronate monooxygenase"),
        ("photolyase/cryptochrome", r"photolyase|cryptochrome"),
        ("flavin reductase/dehydrogenase", r"flavin reductase|flavin-dependent dehydrogenase|flavin dependent dehydrogenase|d-arginine dehydrogenase|fdh\b|dehydrogenase"),
        ("PQQ enzyme", r"\bpqq\b|pyrroloquinoline quinone"),
        ("PLP enzyme", r"\bplp\b|pyridoxal phosphate"),
        ("heme/P450 enzyme", r"\bp450\b|\bheme\b|\bhaem\b"),
    ]
    for label, pattern in checks:
        if re.search(pattern, text):
            return label, label
    family = norm_space(row.get("enzyme_family", "")) or "unknown"
    return family, family


def choose_sequence_route(row: Dict[str, str], names: Sequence[str], organisms: Sequence[str], pdb_ids: Sequence[str], dependency: str) -> str:
    if dependency.startswith("non_flavin_"):
        return "exclude_non_flavin"
    if pdb_ids:
        return "pdb_first"
    if "pdb available" in row_context(row).lower():
        return "pdb_search_first"
    if names and organisms:
        return "uniprot_by_name_organism"
    if names:
        return "paper_supplement_then_database"
    return "manual_literature_review"


def apply_known_override(candidate: Dict[str, str], doi: str) -> None:
    override = KNOWN_PAPER_OVERRIDES.get(doi)
    if not override:
        return
    for key, value in override.items():
        if key in {"pdb_ids", "uniprot_id"} and candidate.get(key) and value:
            candidate[key] = " | ".join(dedupe(split_multi(candidate[key]) + split_multi(value)))
        elif key in {"organism_candidates", "organism_evidence"} and value:
            candidate[key] = value
        else:
            candidate[key] = value
    canonical = candidate.get("canonical_enzyme") or parent_enzyme_name(split_multi(candidate.get("specific_enzyme_names", ""))[0])
    organisms = split_multi(candidate.get("organism_candidates", ""))
    organism = organisms[0] if organisms else ""
    urls = query_urls(canonical, organism)
    purls = pdb_urls(split_multi(candidate.get("pdb_ids", "")), " ".join(t for t in [canonical, organism, candidate.get("title", "")] if t))
    candidate.update(urls)
    candidate.update(purls)
    candidate["uniprot_entry_url"] = uniprot_entry_url(candidate.get("uniprot_id", ""))
    candidate["protein_query"] = " ".join(t for t in [canonical, organism] if t)
    candidate["dna_query"] = " ".join(t for t in [canonical, organism, "gene"] if t)


def tree_ready_row(row: Dict[str, str]) -> Dict[str, str]:
    uniprot_ids = split_multi(row.get("uniprot_id", ""))
    pdb_ids = split_multi(row.get("pdb_ids", ""))
    seed_rank = score_tree_seed(row)
    return {
        "seed_rank": seed_rank,
        "include_in_core_fasta": "",
        "enzyme_short_name": row.get("canonical_enzyme", ""),
        "family": row.get("enzyme_family", ""),
        "function_class": row.get("enzyme_function_class", ""),
        "source_organism": row.get("organism_candidates", ""),
        "uniprot_id": uniprot_ids[0] if uniprot_ids else "",
        "pdb_id": pdb_ids[0] if pdb_ids else "",
        "key_doi": row.get("doi", ""),
        "paper_title": row.get("title", ""),
        "flavin_dependency_status": row.get("flavin_dependency_status", ""),
        "cofactor_class": row.get("cofactor_class", ""),
        "sequence_resolution_route": row.get("sequence_resolution_route", ""),
        "evidence_source": row.get("organism_evidence", ""),
        "evidence_note": row.get("resolution_reason", ""),
        "protein_query": row.get("protein_query", ""),
        "uniprot_url": row.get("uniprot_entry_url") or row.get("uniprot_url", ""),
        "pdb_url": row.get("pdb_url", ""),
        "reaction_type": row.get("reaction_type", ""),
        "new_to_nature_reaction": row.get("new_to_nature_reaction", ""),
        "manual_notes": "",
    }


def score_tree_seed(row: Dict[str, str]) -> str:
    if row.get("flavin_dependency_status") != "confirmed_flavin_dependent":
        return "review"
    if row.get("resolution_status") != "specific_with_organism":
        return "review"
    if row.get("uniprot_id") and row.get("pdb_ids"):
        return "A_pdb_uniprot"
    if row.get("uniprot_id"):
        return "B_uniprot"
    if row.get("pdb_ids"):
        return "B_pdb"
    return "C_name_organism"


def build_tree_ready_rows(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    candidates = [
        tree_ready_row(row)
        for row in rows
        if row.get("flavin_dependency_status") in {"confirmed_flavin_dependent", "likely_flavin_dependent"}
        and row.get("resolution_status") != "non_flavin_exclude"
    ]
    rank_order = {"A_pdb_uniprot": 0, "B_uniprot": 1, "B_pdb": 2, "C_name_organism": 3, "review": 4}
    return sorted(candidates, key=lambda r: (rank_order.get(r["seed_rank"], 9), r.get("enzyme_short_name", ""), r.get("key_doi", "")))


def classify_row(row: Dict[str, str], names: Sequence[str], organisms: Sequence[str]) -> Tuple[str, str]:
    title = row.get("title", "")
    base_name = row.get("enzyme_name") or row.get("enzyme_name_or_target") or ""
    dependency = row.get("flavin_dependency_status", "")
    if dependency.startswith("non_flavin_"):
        return "non_flavin_exclude", "Non-flavin enzyme/cofactor system; not a FASTA seed for the flavin-dependent goal."
    if any(hint in title.lower() for hint in REVIEW_HINTS):
        return "review_or_family_context", "Review/background paper; use for leads, not direct FASTA seed."
    if not names and is_generic_name(base_name, title):
        return "family_only", "No concrete enzyme symbol/name extracted."
    if names and organisms:
        return "specific_with_organism", "Concrete enzyme plus organism/source evidence."
    if names:
        return "specific_missing_organism", "Concrete enzyme candidate but organism/source needs database or paper check."
    return "needs_manual_review", "Insufficient sequence-resolution evidence."


def query_urls(name: str, organism: str) -> Dict[str, str]:
    terms = " ".join(t for t in [parent_enzyme_name(name), organism] if t)
    protein_terms = terms or name or organism
    dna_terms = " ".join(t for t in [parent_enzyme_name(name), organism, "gene"] if t)
    return {
        "uniprot_url": "https://www.uniprot.org/uniprotkb?query=" + quote_plus(protein_terms),
        "ncbi_protein_url": "https://www.ncbi.nlm.nih.gov/protein/?term=" + quote_plus(protein_terms),
        "ncbi_nucleotide_url": "https://www.ncbi.nlm.nih.gov/nuccore/?term=" + quote_plus(dna_terms),
        "ena_url": "https://www.ebi.ac.uk/ena/browser/text-search?query=" + quote_plus(dna_terms),
    }


def uniprot_entry_url(uniprot_id: str) -> str:
    ids = split_multi(uniprot_id)
    if not ids:
        return ""
    return " | ".join(f"https://www.uniprot.org/uniprotkb/{uid}/entry" for uid in ids)


def combine_rows(enzyme_rows: List[Dict[str, str]], homolog_rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    rows = []
    by_doi = defaultdict(dict)
    for row in homolog_rows:
        doi = norm_space(row.get("doi", "")).lower()
        if doi:
            by_doi[doi].update(row)
    for row in enzyme_rows:
        merged = dict(by_doi.get(norm_space(row.get("doi", "")).lower(), {}))
        merged.update(row)
        rows.append(merged)
    for row in homolog_rows:
        doi = norm_space(row.get("doi", "")).lower()
        if doi and not any(norm_space(r.get("doi", "")).lower() == doi for r in enzyme_rows):
            rows.append(dict(row))
    return rows


def build_candidates(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    out = []
    for row in rows:
        doi = norm_space(row.get("doi", "")).lower()
        names = extract_specific_names(row)
        mutations = extract_mutations(names, row)
        organisms, organism_evidence = extract_organisms(row, names)
        pdb_ids, pdb_evidence = extract_pdb_ids(row)
        cofactor = infer_cofactor_and_function(row, names)
        row_for_status = dict(row)
        row_for_status.update(cofactor)
        status, reason = classify_row(row_for_status, names, organisms)
        route = choose_sequence_route(row_for_status, names, organisms, pdb_ids, cofactor["flavin_dependency_status"])
        name_for_query = parent_enzyme_name(names[0]) if names else clean_candidate_name(row.get("enzyme_name") or row.get("enzyme_name_or_target", ""))
        organism_for_query = organisms[0] if organisms else ""
        urls = query_urls(name_for_query, organism_for_query)
        pdb_link_query = " ".join(t for t in [name_for_query, organism_for_query, row.get("title", "")] if t)
        purls = pdb_urls(pdb_ids, pdb_link_query)
        candidate = {
            "resolution_status": status,
            "resolution_reason": reason,
            "specific_enzyme_names": " | ".join(names),
            "parent_enzyme_names": " | ".join(dedupe(parent_enzyme_name(n) for n in names)),
            "mutations_or_variants": " | ".join(mutations),
            "organism_candidates": " | ".join(organisms),
            "organism_evidence": " | ".join(organism_evidence),
            "protein_query": " ".join(t for t in [name_for_query, organism_for_query] if t),
            "dna_query": " ".join(t for t in [name_for_query, organism_for_query, "gene"] if t),
            **urls,
            **purls,
            "enzyme_family": row.get("enzyme_family", ""),
            "enzyme_name_original": row.get("enzyme_name") or row.get("enzyme_name_or_target", ""),
            "doi": row.get("doi", ""),
            "title": row.get("title", ""),
            "year": row.get("year", ""),
            "journal": row.get("journal", ""),
            "search_track": row.get("search_track", ""),
            "flavin_cofactor": row.get("flavin_cofactor", ""),
            "reaction_type": row.get("reaction_type", ""),
            "new_to_nature_reaction": row.get("new_to_nature_reaction", ""),
            "characterized_enzyme_evidence": row.get("characterized_enzyme_evidence", ""),
            "structure_similarity_hint": row.get("structure_similarity_hint", ""),
            "criteria_status": row.get("criteria_status", ""),
            "manual_review_priority": row.get("manual_review_priority", row.get("seed_priority", "")),
            "pdb_ids": " | ".join(pdb_ids),
            "pdb_evidence": " | ".join(pdb_evidence),
            "uniprot_id": row.get("uniprot_id", ""),
            "uniprot_entry_url": uniprot_entry_url(row.get("uniprot_id", "")),
            "ec_number_candidates": " | ".join(extract_ec_numbers(row)),
            "sequence_resolution_route": route,
            **cofactor,
            "evidence_summary": row.get("evidence_summary") or row.get("abstract", ""),
        }
        apply_known_override(candidate, doi)
        out.append(candidate)
    return sorted(dedupe_candidates(out), key=sort_key)


def canonical_enzyme_name(row: Dict[str, str]) -> str:
    if row.get("canonical_enzyme"):
        return row["canonical_enzyme"]
    blob = " | ".join([
        row.get("specific_enzyme_names", ""),
        row.get("parent_enzyme_names", ""),
        row.get("enzyme_name_original", ""),
    ])
    priority = [
        r"\bPfBAL\b", r"\bCvFAP\b", r"\bGluER\b", r"\bGsOYE\b", r"\bOaER\b", r"\bOYE1\b", r"\bPaDADH\b",
        r"\bPqsL\b", r"\bD-2-hydroxyglutarate dehydrogenase\b", r"\b[A-Z][a-z]{1,4}ER\b", r"\b[A-Z][a-z]OYE\d*\b",
        r"\b[A-Z][a-z]{1,5}FDH\b", r"\b[A-Z][a-z]{1,5}NTR\b",
        r"\bbenzaldehyde lyase\b",
        r"\bPseudomonas aeruginosa D-arginine dehydrogenase\b",
        r"\bovenolide biosynthetic flavoenzyme\b",
    ]
    for pattern in priority:
        match = re.search(pattern, blob)
        if match:
            return match.group(0)
    parents = split_multi(row.get("parent_enzyme_names", ""))
    if parents:
        return parents[0]
    names = split_multi(row.get("specific_enzyme_names", ""))
    if names:
        return parent_enzyme_name(names[0])
    return parent_enzyme_name(row.get("enzyme_name_original", ""))


def dedupe_candidates(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    grouped: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    status_rank = {
        "specific_with_organism": 0,
        "specific_missing_organism": 1,
        "needs_manual_review": 2,
        "family_only": 3,
        "review_or_family_context": 4,
        "non_flavin_exclude": 5,
    }
    for row in rows:
        canonical = canonical_enzyme_name(row)
        organism = split_multi(row.get("organism_candidates", ""))
        organism_key = organism[0] if organism else ""
        key = (norm_space(row.get("doi", "")).lower(), canonical.lower(), organism_key.lower())
        if key not in grouped:
            merged = dict(row)
            merged["canonical_enzyme"] = canonical
            grouped[key] = merged
            continue
        current = grouped[key]
        if status_rank.get(row["resolution_status"], 9) < status_rank.get(current["resolution_status"], 9):
            current["resolution_status"] = row["resolution_status"]
            current["resolution_reason"] = row["resolution_reason"]
        for field in [
            "specific_enzyme_names", "parent_enzyme_names", "mutations_or_variants",
            "organism_candidates", "organism_evidence", "pdb_ids", "uniprot_id",
            "pdb_evidence", "ec_number_candidates", "cofactor_detail", "cofactor_evidence",
        ]:
            current[field] = " | ".join(dedupe(split_multi(current.get(field, "")) + split_multi(row.get(field, ""))))
        for field in [
            "cofactor_class", "flavin_dependency_status", "enzyme_function_class",
            "sequence_resolution_route", "pdb_url", "rcsb_search_url",
        ]:
            if not current.get(field) and row.get(field):
                current[field] = row[field]
        if len(norm_space(row.get("evidence_summary", ""))) > len(norm_space(current.get("evidence_summary", ""))):
            current["evidence_summary"] = row.get("evidence_summary", "")
        if not current.get("title") and row.get("title"):
            current["title"] = row["title"]
    for row in grouped.values():
        canonical = row.get("canonical_enzyme") or canonical_enzyme_name(row)
        organisms = split_multi(row.get("organism_candidates", ""))
        organism = organisms[0] if organisms else ""
        urls = query_urls(canonical, organism)
        purls = pdb_urls(split_multi(row.get("pdb_ids", "")), " ".join(t for t in [canonical, organism, row.get("title", "")] if t))
        row.update(urls)
        row.update(purls)
        row["protein_query"] = " ".join(t for t in [canonical, organism] if t)
        row["dna_query"] = " ".join(t for t in [canonical, organism, "gene"] if t)
    return list(grouped.values())


def sort_key(row: Dict[str, str]) -> Tuple[int, int, int, str]:
    status_score = {
        "specific_with_organism": 0,
        "specific_missing_organism": 1,
        "needs_manual_review": 2,
        "family_only": 3,
        "review_or_family_context": 4,
        "non_flavin_exclude": 5,
    }.get(row["resolution_status"], 5)
    dependency_score = {
        "confirmed_flavin_dependent": 0,
        "likely_flavin_dependent": 1,
        "mixed_or_review": 2,
        "unknown": 3,
        "unknown_or_nadph_only": 3,
    }.get(row.get("flavin_dependency_status", ""), 4)
    priority_score = 0 if row.get("manual_review_priority") == "high" else 1
    new_score = 0 if str(row.get("new_to_nature_reaction", "")).lower() == "true" else 1
    return (dependency_score, status_score, priority_score, new_score, row.get("title", ""))


def render_html(rows: Sequence[Dict[str, str]], output_path: Path) -> None:
    counts = defaultdict(int)
    dep_counts = defaultdict(int)
    cofactor_counts = defaultdict(int)
    function_counts = defaultdict(int)
    for row in rows:
        counts[row["resolution_status"]] += 1
        dep_counts[row.get("flavin_dependency_status", "unknown")] += 1
        cofactor_counts[row.get("cofactor_class", "unknown")] += 1
        function_counts[row.get("enzyme_function_class", "unknown")] += 1
    cards = "\n".join([
        "<div class='stat stat-total'><b>Total rows</b><span>{}</span></div>".format(len(rows)),
        *[
            f"<div class='stat'><b>{html.escape(k)}</b><span>{v}</span></div>"
            for k, v in sorted(dep_counts.items(), key=lambda item: (-item[1], item[0]))[:8]
        ],
    ])
    status_options = options_from(rows, "resolution_status", "All readiness")
    dependency_options = options_from(rows, "flavin_dependency_status", "All dependencies")
    cofactor_options = options_from(rows, "cofactor_class", "All cofactors")
    function_options = options_from(rows, "enzyme_function_class", "All functions")
    route_options = options_from(rows, "sequence_resolution_route", "All routes")
    table_rows = "\n".join(render_row(row) for row in rows)
    doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Sequence Seed Review</title>
<style>
:root {{ color-scheme: light; --ink:#17202a; --muted:#64748b; --line:#d6dee8; --bg:#f5f7fa; --panel:#ffffff; --soft:#eef4f7; --accent:#006d77; --good:#0f766e; --warn:#9a5b00; --bad:#9f1239; --blue:#285a8d; --violet:#6d4bb3; --amber:#b7791f; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; font:13px/1.42 system-ui, -apple-system, Segoe UI, sans-serif; color:var(--ink); background:var(--bg); }}
header {{ padding:18px 24px 12px; background:linear-gradient(180deg, #ffffff 0%, #f7fbfc 100%); border-bottom:1px solid var(--line); position:sticky; top:0; z-index:5; box-shadow:0 10px 28px rgba(22, 34, 51, .06); }}
h1 {{ margin:0; font-size:21px; letter-spacing:0; }}
.sub {{ color:var(--muted); max-width:1160px; margin-top:5px; }}
.stats {{ display:grid; grid-template-columns:repeat(auto-fit, minmax(178px, 1fr)); gap:8px; margin-top:12px; }}
.stat {{ border:1px solid var(--line); background:#fbfdff; padding:8px 10px; border-radius:6px; min-width:0; display:flex; justify-content:space-between; gap:12px; box-shadow:0 2px 10px rgba(15, 23, 42, .035); transition:transform .16s ease, box-shadow .16s ease, border-color .16s ease; }}
.stat:hover {{ transform:translateY(-1px); box-shadow:0 8px 18px rgba(15, 23, 42, .08); border-color:#9ccfca; }}
.stat b {{ overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }}
.stat span {{ font-weight:700; color:var(--accent); }}
.stat-total {{ background:#e8f4f3; border-color:#9ccfca; }}
.toolbar {{ display:grid; grid-template-columns:minmax(260px, 1.6fr) repeat(6, minmax(150px, 1fr)); gap:8px; padding:10px 24px; background:rgba(238,244,247,.94); backdrop-filter:blur(8px); border-bottom:1px solid var(--line); position:sticky; top:129px; z-index:4; }}
input, select {{ height:34px; min-width:0; border:1px solid var(--line); border-radius:6px; padding:0 9px; background:white; color:var(--ink); transition:border-color .16s ease, box-shadow .16s ease, transform .16s ease; }}
input:focus, select:focus {{ outline:none; border-color:#5aa7a7; box-shadow:0 0 0 3px rgba(0, 109, 119, .14); }}
select:hover, input:hover {{ border-color:#9bb4c8; }}
main {{ padding:16px 24px 30px; }}
table {{ border-collapse:collapse; width:100%; background:var(--panel); border:1px solid var(--line); table-layout:fixed; box-shadow:0 18px 45px rgba(15, 23, 42, .07); }}
th, td {{ border-bottom:1px solid var(--line); vertical-align:top; padding:8px 9px; word-break:break-word; }}
th {{ text-align:left; background:#f1f5f8; position:sticky; top:184px; z-index:3; font-size:12px; color:#334155; }}
tbody tr {{ animation:rowIn .26s ease both; border-left:4px solid transparent; transition:transform .12s ease, box-shadow .12s ease, background .12s ease, border-color .12s ease; }}
tbody tr:hover {{ transform:translateX(2px); box-shadow:inset 4px 0 0 #9ccfca; }}
tbody tr[data-dependency="confirmed_flavin_dependent"] {{ border-left-color:#38a169; }}
tbody tr[data-dependency="likely_flavin_dependent"] {{ border-left-color:#3182ce; }}
tbody tr[data-dependency="mixed_or_review"] {{ border-left-color:#d69e2e; }}
tbody tr[data-dependency^="non_flavin"] {{ border-left-color:#d53f8c; }}
tr:hover td {{ background:#fbfdff; }}
.pill {{ display:inline-block; padding:2px 7px; border-radius:999px; font-size:12px; border:1px solid var(--line); background:#f8fafc; white-space:nowrap; margin:0 4px 4px 0; transition:transform .14s ease; }}
.pill:hover {{ transform:translateY(-1px); }}
.specific_with_organism, .confirmed_flavin_dependent {{ color:var(--good); border-color:#86c5b9; background:#eaf8f5; }}
.likely_flavin_dependent, .specific_missing_organism {{ color:var(--blue); border-color:#9ebee1; background:#edf5ff; }}
.mixed_or_review, .unknown, .unknown_or_nadph_only, .family_only, .review_or_family_context {{ color:var(--warn); border-color:#e4bd75; background:#fff7e6; }}
.non_flavin_exclude, .non_flavin_thdp, .non_flavin_pqq, .non_flavin_plp, .non_flavin_heme {{ color:var(--bad); border-color:#e5a3af; background:#fff0f3; }}
.supplement_accession_first, .pdb_first {{ color:var(--violet); border-color:#c2b5ee; background:#f3f0ff; }}
.pdb_search_first {{ color:#5f6b7a; border-color:#cbd5e1; background:#f8fafc; }}
.thdp {{ color:#8a4b18; border-color:#f1bf98; background:#fff2e8; }}
.flavin {{ color:#26655a; border-color:#a4d7cb; background:#ecfbf7; }}
.small {{ color:var(--muted); font-size:12px; margin-top:3px; }}
.title {{ width:25%; }}
.enzyme {{ width:19%; }}
.dependency {{ width:17%; }}
.links {{ width:18%; }}
.evidence {{ max-height:9.2em; overflow:auto; }}
a {{ color:var(--accent); text-decoration:none; font-weight:500; }}
a:hover {{ text-decoration:underline; }}
.links a {{ display:block; margin-bottom:3px; }}
@keyframes rowIn {{ from {{ opacity:0; transform:translateY(4px); }} to {{ opacity:1; transform:translateY(0); }} }}
@media (max-width: 980px) {{
  header {{ position:static; }}
  .toolbar {{ position:static; grid-template-columns:1fr 1fr; }}
  th {{ position:static; }}
  table {{ table-layout:auto; min-width:980px; }}
  main {{ overflow:auto; }}
}}
</style>
</head>
<body>
<header>
  <h1>Sequence Seed Review</h1>
  <div class="sub">Offline evidence table for flavin-dependence, enzyme function, organism source, PDB-first routes, and FASTA seed readiness.</div>
  <div class="stats">{cards}</div>
</header>
<div class="toolbar">
  <input id="q" placeholder="Search title, enzyme, organism, DOI, reaction">
  <select id="dependency">{dependency_options}</select>
  <select id="cofactor">{cofactor_options}</select>
  <select id="functionClass">{function_options}</select>
  <select id="status">{status_options}</select>
  <select id="route">{route_options}</select>
  <select id="pdb"><option value="">All PDB</option><option value="yes">PDB present</option><option value="no">No PDB</option></select>
</div>
<main>
<table id="tbl">
<thead><tr>
<th class="dependency">Dependency</th><th class="enzyme">Specific Enzyme</th><th>Organism / Route</th><th class="links">Sequence Links</th><th class="title">Paper / Evidence</th><th>Reaction</th>
</tr></thead>
<tbody>
{table_rows}
</tbody>
</table>
</main>
<script>
const q = document.getElementById('q');
const status = document.getElementById('status');
const dependency = document.getElementById('dependency');
const cofactor = document.getElementById('cofactor');
const functionClass = document.getElementById('functionClass');
const route = document.getElementById('route');
const pdb = document.getElementById('pdb');
const rows = [...document.querySelectorAll('tbody tr')];
function apply() {{
  const needle = q.value.toLowerCase();
  const filters = {{
    status: status.value,
    dependency: dependency.value,
    cofactor: cofactor.value,
    functionClass: functionClass.value,
    route: route.value,
  }};
  rows.forEach(row => {{
    const okText = !needle || row.innerText.toLowerCase().includes(needle);
    const okStatus = !filters.status || row.dataset.status === filters.status;
    const okDep = !filters.dependency || row.dataset.dependency === filters.dependency;
    const okCofactor = !filters.cofactor || row.dataset.cofactor === filters.cofactor;
    const okFunction = !filters.functionClass || row.dataset.functionClass === filters.functionClass;
    const okRoute = !filters.route || row.dataset.route === filters.route;
    const okPdb = !pdb.value || row.dataset.pdb === pdb.value;
    row.style.display = okText && okStatus && okDep && okCofactor && okFunction && okRoute && okPdb ? '' : 'none';
  }});
}}
q.addEventListener('input', apply);
status.addEventListener('input', apply);
dependency.addEventListener('input', apply);
cofactor.addEventListener('input', apply);
functionClass.addEventListener('input', apply);
route.addEventListener('input', apply);
pdb.addEventListener('input', apply);
</script>
</body></html>"""
    output_path.write_text(doc, encoding="utf-8")


def options_from(rows: Sequence[Dict[str, str]], field: str, label: str) -> str:
    values = sorted({norm_space(row.get(field, "")) or "unknown" for row in rows})
    options = [f'<option value="">{html.escape(label)}</option>']
    options.extend(f'<option>{html.escape(value)}</option>' for value in values)
    return "".join(options)


def render_row(row: Dict[str, str]) -> str:
    doi = row.get("doi", "")
    doi_link = f"<a href='https://doi.org/{html.escape(doi)}'>{html.escape(doi)}</a>" if doi else ""
    evidence = norm_space(row.get("evidence_summary", ""))[:520]
    dep = row.get("flavin_dependency_status", "unknown")
    cofactor = row.get("cofactor_class", "unknown")
    function_class = row.get("enzyme_function_class", "unknown")
    route = row.get("sequence_resolution_route", "manual_literature_review")
    pdb_ids = split_multi(row.get("pdb_ids", ""))
    pdb_present = "yes" if pdb_ids else "no"
    pdb_links = "".join(
        f'<a href="https://www.rcsb.org/structure/{html.escape(pid)}">PDB {html.escape(pid)}</a>'
        for pid in pdb_ids
    )
    if not pdb_links and row.get("rcsb_search_url"):
        pdb_links = f'<a href="{html.escape(row["rcsb_search_url"])}">RCSB search</a>'
    uniprot_ids = split_multi(row.get("uniprot_id", ""))
    uniprot_id_links = "".join(
        f'<a href="https://www.uniprot.org/uniprotkb/{html.escape(uid)}/entry">UniProt {html.escape(uid)}</a>'
        for uid in uniprot_ids
    )
    return f"""<tr data-status="{html.escape(row['resolution_status'])}" data-dependency="{html.escape(dep)}" data-cofactor="{html.escape(cofactor)}" data-function-class="{html.escape(function_class)}" data-route="{html.escape(route)}" data-pdb="{pdb_present}">
<td><span class="pill {html.escape(dep)}">{html.escape(dep)}</span><span class="pill {html.escape(cofactor)}">{html.escape(cofactor)}</span><div class="small">{html.escape(row.get('cofactor_detail',''))}</div><div class="small">{html.escape(row.get('cofactor_evidence',''))}</div><div class="small">EC: {html.escape(row.get('ec_number_candidates',''))}</div></td>
<td class="enzyme"><b>{html.escape(row.get('canonical_enzyme') or row.get('specific_enzyme_names') or row.get('enzyme_name_original',''))}</b><div class="small">Names: {html.escape(row.get('specific_enzyme_names',''))}</div><div class="small">Parent: {html.escape(row.get('parent_enzyme_names',''))}</div><div class="small">Variants: {html.escape(row.get('mutations_or_variants',''))}</div><div class="small">Function: {html.escape(function_class)}</div><div class="small">Family: {html.escape(row.get('enzyme_family',''))}</div></td>
<td><span class="pill {html.escape(row['resolution_status'])}">{html.escape(row['resolution_status'])}</span><div class="small">{html.escape(row['resolution_reason'])}</div><b>{html.escape(row.get('organism_candidates',''))}</b><div class="small">{html.escape(row.get('organism_evidence',''))}</div><div class="small">Route: {html.escape(route)}</div></td>
<td class="links">{uniprot_id_links}<a href="{html.escape(row['uniprot_url'])}">UniProt search</a><a href="{html.escape(row['ncbi_protein_url'])}">NCBI Protein</a><a href="{html.escape(row['ncbi_nucleotide_url'])}">NCBI Nucleotide/gene</a><a href="{html.escape(row['ena_url'])}">ENA text search</a>{pdb_links}<div class="small">Protein: {html.escape(row.get('protein_query',''))}</div><div class="small">DNA: {html.escape(row.get('dna_query',''))}</div></td>
<td class="title"><b>{html.escape(row.get('title',''))}</b><div>{doi_link}</div><div class="small">{html.escape(str(row.get('year','')))} {html.escape(row.get('journal',''))}</div><div class="evidence small">{html.escape(evidence)}</div></td>
<td>{html.escape(row.get('reaction_type',''))}<div class="small">new-to-nature={html.escape(str(row.get('new_to_nature_reaction','')))}</div><div class="small">{html.escape(row.get('characterized_enzyme_evidence',''))}</div><div class="small">{html.escape(row.get('structure_similarity_hint',''))}</div></td>
</tr>"""


def write_xlsx(path: Path, rows: Sequence[Dict[str, str]], fields: Sequence[str]) -> None:
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font, PatternFill, Border, Side
        from openpyxl.formatting.rule import FormulaRule
        from openpyxl.utils import get_column_letter
    except ImportError:
        print("[SequenceSeed] openpyxl is unavailable; skipped XLSX export.")
        return

    wb = Workbook()
    default = wb.active
    wb.remove(default)
    tree_rows = build_tree_ready_rows(rows)
    manual_rows = [row for row in rows if row.get("resolution_status") in {"family_only", "specific_missing_organism", "review_or_family_context"}]
    excluded_rows = [row for row in rows if row.get("resolution_status") == "non_flavin_exclude" or row.get("flavin_dependency_status", "").startswith("non_flavin_")]

    tree_fields = [
        "seed_rank", "include_in_core_fasta", "enzyme_short_name", "family", "function_class",
        "source_organism", "uniprot_id", "pdb_id", "key_doi", "paper_title",
        "flavin_dependency_status", "cofactor_class", "sequence_resolution_route",
        "evidence_source", "evidence_note", "protein_query", "uniprot_url", "pdb_url",
        "reaction_type", "new_to_nature_reaction", "manual_notes",
    ]
    add_sheet(wb, "tree_ready_core", tree_rows, tree_fields, freeze="A2")
    add_sheet(wb, "all_candidates", rows, fields, freeze="A2")
    add_sheet(wb, "manual_review", manual_rows, fields, freeze="A2")
    add_sheet(wb, "excluded_non_flavin", excluded_rows, fields, freeze="A2")

    fills = {
        "header": PatternFill("solid", fgColor="1F2937"),
        "a": PatternFill("solid", fgColor="DDF7EC"),
        "b": PatternFill("solid", fgColor="E5F0FF"),
        "c": PatternFill("solid", fgColor="FFF3CF"),
        "exclude": PatternFill("solid", fgColor="FFE4EA"),
    }
    thin = Side(style="thin", color="D9E2EC")
    border = Border(bottom=thin)
    for ws in wb.worksheets:
        ws.sheet_view.showGridLines = False
        ws.auto_filter.ref = ws.dimensions
        for cell in ws[1]:
            cell.fill = fills["header"]
            cell.font = Font(color="FFFFFF", bold=True)
            cell.alignment = Alignment(vertical="center", wrap_text=True)
            cell.border = border
        for row_cells in ws.iter_rows(min_row=2):
            for cell in row_cells:
                cell.alignment = Alignment(vertical="top", wrap_text=True)
                cell.border = border
        for col_idx, width in enumerate(preferred_widths(ws), start=1):
            ws.column_dimensions[get_column_letter(col_idx)].width = width
        if ws.title == "tree_ready_core" and ws.max_row > 1:
            ws.conditional_formatting.add(
                f"A2:A{ws.max_row}",
                FormulaRule(formula=['LEFT($A2,1)="A"'], fill=fills["a"]),
            )
            ws.conditional_formatting.add(
                f"A2:A{ws.max_row}",
                FormulaRule(formula=['LEFT($A2,1)="B"'], fill=fills["b"]),
            )
            ws.conditional_formatting.add(
                f"A2:A{ws.max_row}",
                FormulaRule(formula=['LEFT($A2,1)="C"'], fill=fills["c"]),
            )
        if ws.title == "excluded_non_flavin" and ws.max_row > 1:
            ws.conditional_formatting.add(f"A2:A{ws.max_row}", FormulaRule(formula=["TRUE"], fill=fills["exclude"]))
    wb.save(path)


def add_sheet(wb, title: str, rows: Sequence[Dict[str, str]], fields: Sequence[str], freeze: str = "A2") -> None:
    ws = wb.create_sheet(title)
    ws.append(list(fields))
    for row in rows:
        ws.append([row.get(field, "") for field in fields])
    ws.freeze_panes = freeze


def preferred_widths(ws) -> List[int]:
    widths: List[int] = []
    for col in ws.iter_cols():
        header = str(col[0].value or "")
        max_len = max([len(str(cell.value or "")) for cell in col[:80]] + [len(header)])
        if header in {"paper_title", "evidence_summary", "evidence_note", "evidence_source"}:
            widths.append(min(max(max_len + 2, 28), 58))
        elif header.endswith("_url"):
            widths.append(32)
        elif header in {"include_in_core_fasta", "new_to_nature_reaction"}:
            widths.append(18)
        else:
            widths.append(min(max(max_len + 2, 12), 34))
    return widths


def main() -> None:
    parser = argparse.ArgumentParser(description="Build sequence seed review HTML from literature-agent CSV outputs.")
    parser.add_argument("--enzyme-csv", required=True, help="Path to enzyme_seed_candidates.csv")
    parser.add_argument("--homolog-csv", help="Optional path to homolog_candidate_seeds.csv")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--prefix", default="sequence_seed", help="Output file prefix")
    args = parser.parse_args()

    enzyme_rows = read_csv(Path(args.enzyme_csv))
    homolog_rows = read_csv(Path(args.homolog_csv)) if args.homolog_csv else []
    rows = build_candidates(combine_rows(enzyme_rows, homolog_rows))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fields = [
        "resolution_status", "resolution_reason", "specific_enzyme_names",
        "canonical_enzyme", "parent_enzyme_names", "mutations_or_variants", "organism_candidates",
        "organism_evidence", "protein_query", "dna_query", "uniprot_url",
        "ncbi_protein_url", "ncbi_nucleotide_url", "ena_url", "pdb_url", "rcsb_search_url", "uniprot_entry_url",
        "cofactor_class", "cofactor_detail", "cofactor_evidence", "flavin_dependency_status",
        "enzyme_function_class", "ec_number_candidates", "sequence_resolution_route", "enzyme_family",
        "enzyme_name_original", "doi", "title", "year", "journal", "search_track",
        "flavin_cofactor", "reaction_type", "new_to_nature_reaction",
        "characterized_enzyme_evidence", "structure_similarity_hint", "criteria_status",
        "manual_review_priority", "pdb_ids", "pdb_evidence", "uniprot_id", "evidence_summary",
    ]
    csv_path = outdir / f"{args.prefix}_candidates.csv"
    tree_csv_path = outdir / f"{args.prefix}_tree_ready.csv"
    xlsx_path = outdir / f"{args.prefix}_tree_ready.xlsx"
    html_path = outdir / f"{args.prefix}_review.html"
    write_csv(csv_path, rows, fields)
    tree_fields = [
        "seed_rank", "include_in_core_fasta", "enzyme_short_name", "family", "function_class",
        "source_organism", "uniprot_id", "pdb_id", "key_doi", "paper_title",
        "flavin_dependency_status", "cofactor_class", "sequence_resolution_route",
        "evidence_source", "evidence_note", "protein_query", "uniprot_url", "pdb_url",
        "reaction_type", "new_to_nature_reaction", "manual_notes",
    ]
    write_csv(tree_csv_path, build_tree_ready_rows(rows), tree_fields)
    write_xlsx(xlsx_path, rows, fields)
    render_html(rows, html_path)
    print(f"[SequenceSeed] Wrote {csv_path} ({len(rows)} rows)")
    print(f"[SequenceSeed] Wrote {tree_csv_path}")
    print(f"[SequenceSeed] Wrote {xlsx_path}")
    print(f"[SequenceSeed] Wrote {html_path}")


if __name__ == "__main__":
    main()
