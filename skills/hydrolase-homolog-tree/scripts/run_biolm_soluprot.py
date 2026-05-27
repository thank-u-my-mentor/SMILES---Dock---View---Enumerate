#!/usr/bin/env python3
"""
Run BioLM SoluProt API on a FASTA file and write CSV predictions.

Requires:
  export BIOLM_API_KEY="..."

The API accepts at most 100 sequences per request, so this batches inputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import time
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


API_URL = "https://biolm.ai/api/v3/soluprot/predict/"


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
                    records.append((header, clean_sequence("".join(seq))))
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if header and seq:
        records.append((header, clean_sequence("".join(seq))))
    return records


def clean_sequence(seq: str) -> str:
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())


def chunks(records: List[Tuple[str, str]], size: int) -> Iterable[List[Tuple[str, str]]]:
    for i in range(0, len(records), size):
        yield records[i:i + size]


def post_json(url: str, token: str, payload: dict, timeout: int, retries: int) -> dict:
    data = json.dumps(payload).encode("utf-8")
    for attempt in range(retries + 1):
        req = Request(
            url,
            data=data,
            headers={
                "Authorization": f"Token {token}",
                "Content-Type": "application/json",
                "Accept": "application/json",
            },
            method="POST",
        )
        try:
            with urlopen(req, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")[:500]
            if attempt == retries:
                raise RuntimeError(f"BioLM API HTTP {exc.code}: {body}") from exc
        except URLError as exc:
            if attempt == retries:
                raise RuntimeError(f"BioLM API connection error: {exc}") from exc
        time.sleep(2 * (attempt + 1))
    raise RuntimeError("BioLM API failed without a response")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run BioLM SoluProt API on FASTA and write CSV.")
    parser.add_argument("fasta")
    parser.add_argument("--output", default="soluprot_predictions.csv")
    parser.add_argument("--api-key", default=os.getenv("BIOLM_API_KEY"))
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--sleep", type=float, default=0.5)
    args = parser.parse_args()

    if not args.api_key:
        raise SystemExit("Missing API key. Set BIOLM_API_KEY or pass --api-key.")
    if not 1 <= args.batch_size <= 100:
        raise SystemExit("--batch-size must be 1..100 for BioLM SoluProt.")

    records = [(rid, seq) for rid, seq in read_fasta(Path(args.fasta)) if 20 <= len(seq) <= 5000]
    out_rows: List[Dict[str, str]] = []
    for batch_idx, batch in enumerate(chunks(records, args.batch_size), 1):
        print(f"BioLM SoluProt batch {batch_idx}: {len(batch)} sequences")
        payload = {
            "params": {"batch_size": args.batch_size, "min_sequence_len": 20, "max_sequence_len": 5000},
            "items": [{"sequence": seq} for _, seq in batch],
        }
        data = post_json(API_URL, args.api_key, payload, args.timeout, args.retries)
        results = data.get("results", [])
        if len(results) != len(batch):
            raise RuntimeError(f"Expected {len(batch)} results, got {len(results)}: {str(data)[:500]}")
        for (seq_id, seq), result in zip(batch, results):
            score = result.get("soluble", "")
            soluble = result.get("is_soluble", "")
            out_rows.append({
                "id": seq_id,
                "soluprot_probability": score,
                "prediction": "Soluble" if soluble is True else "Insoluble" if soluble is False else "",
                "sequence_length": len(seq),
            })
        time.sleep(args.sleep)

    with Path(args.output).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "soluprot_probability", "prediction", "sequence_length"])
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"Wrote {args.output} with {len(out_rows)} predictions.")


if __name__ == "__main__":
    main()
