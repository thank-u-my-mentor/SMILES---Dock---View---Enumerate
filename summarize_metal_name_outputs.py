from pathlib import Path


base = Path("/home/qin/Metal-F_project")
for name in ["CF3", "C-F", "fluoroaryl"]:
    d = base / f"Metal-{name}"
    ids = (d / f"Metal-{name}_geometry_pass_pdb_ids.txt").read_text().split()
    fasta = d / f"Metal-{name}_hmmer_full_length_seed.fasta"
    text = fasta.read_text() if fasta.exists() else ""
    count = sum(1 for line in text.splitlines() if line.startswith(">"))
    print(f"Metal-{name}")
    print("  pass_ids:", " ".join(ids))
    print("  pass_count:", len(ids))
    print("  fasta_count:", count)
    print("  fasta_size:", fasta.stat().st_size if fasta.exists() else 0)
    for suffix in [
        "geometry_summary.csv",
        "literature_seed.csv",
        "protein_entities_kept.csv",
        "protein_entities_excluded.csv",
        "hmmer_full_length_seed.fasta",
    ]:
        path = d / f"Metal-{name}_{suffix}"
        print(" ", path, path.stat().st_size if path.exists() else "missing")
