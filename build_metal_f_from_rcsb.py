from pathlib import Path


RCSB_SEARCH_URL = "https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_nonpolymer_instance_annotation.comp_id%22%2C%22operator%22%3A%22exact_match%22%2C%22value%22%3A%22F%22%7D%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_nonpolymer_instance_annotation.type%22%2C%22operator%22%3A%22exact_match%22%2C%22value%22%3A%22HAS_NO_COVALENT_LINKAGE%22%2C%22negation%22%3Afalse%7D%7D%5D%2C%22label%22%3A%22nested-attribute%22%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_chem_comp_container_identifiers.comp_id%22%2C%22operator%22%3A%22in%22%2C%22negation%22%3Afalse%2C%22value%22%3A%5B%22FE%22%2C%22CO%22%2C%22NI%22%2C%22CU%22%2C%22MN%22%5D%7D%7D%5D%7D%5D%2C%22label%22%3A%22text%22%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%2C%22return_type%22%3A%22entry%22%2C%22request_options%22%3A%7B%22paginate%22%3A%7B%22start%22%3A0%2C%22rows%22%3A25%7D%2C%22results_content_type%22%3A%5B%22experimental%22%5D%2C%22sort%22%3A%5B%7B%22sort_by%22%3A%22score%22%2C%22direction%22%3A%22desc%22%7D%5D%2C%22scoring_strategy%22%3A%22combined%22%7D%2C%22request_info%22%3A%7B%22query_id%22%3A%22aef02f4fa52770114ba30e25b310e775%22%7D%7D"


rc_path = Path("/home/qin/.pymolrc.py")
exec(compile(rc_path.read_text(), str(rc_path), "exec"), globals(), globals())

build_metal_f_from_rcsb(
    RCSB_SEARCH_URL,
    output="/home/qin/Metal-F_project/Metal-F_rebuilt_scenes.pse",
    cif_dir="/home/qin/Metal-F_project/cif",
    f_cutoff=4.0,
    coord_cutoff=2.8,
    around=4.0,
    only_close=0,
    png_dir="",
    csv_path="/home/qin/Metal-F_project/Metal-F_scene_summary.csv",
    fasta_path="/home/qin/Metal-F_project/Metal-F_protein_entities.fasta",
    entity_csv_path="/home/qin/Metal-F_project/Metal-F_protein_entities.csv",
    literature_csv_path="/home/qin/Metal-F_project/Metal-F_literature_seed.csv",
    overwrite=0,
    ray=0,
)
