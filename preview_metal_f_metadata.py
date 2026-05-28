import csv
import importlib.util
from pathlib import Path


module_path = Path("/mnt/e/Codex/pymolrc_py_updated.py")
source = module_path.read_text(encoding="utf-8")
prefix = source.split("@cmd.extend\ndef simple_hetero", 1)[0]
prefix = prefix.replace("from pymol import cmd, stored\n\n", "")
namespace = {}
exec(prefix, namespace, namespace)

pdb_ids = "1BS3 1E6A 1T5G 2AU6 2AU9 2QFP 3VRS 4CEX 4GOA 5VEJ 6QDY 6RI4 6RI6 6RI8 6RII 6RIK 6RIL 6YTT 7T0Y".split()
rows = []
for pdb_id in pdb_ids:
    meta = namespace["_rcsb_metadata_bundle"](pdb_id)
    rows.append(
        {
            "pdb_id": pdb_id,
            "doi": meta.get("doi", ""),
            "pubmed_id": meta.get("pubmed_id", ""),
            "citation_title": meta.get("citation_title", ""),
            "uniprot_ids": meta.get("uniprot_ids", ""),
            "protein_entity_ids": meta.get("protein_entity_ids", ""),
            "organisms": meta.get("organisms", ""),
        }
    )

out = Path("/home/qin/Metal-F_project/Metal-F_metadata_preview.csv")
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(out)
for row in rows[:5]:
    print(row)
