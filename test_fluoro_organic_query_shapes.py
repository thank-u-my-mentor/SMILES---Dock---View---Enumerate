import json
import urllib.error
import urllib.request


ENDPOINT = "https://search.rcsb.org/rcsbsearch/v2/query"


def run(name, chemical_params):
    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "chemical", "parameters": chemical_params},
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
                        "operator": "in",
                        "value": ["FE", "CO", "NI", "CU", "MN"],
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 10},
            "results_content_type": ["experimental"],
        },
    }
    req = urllib.request.Request(
        ENDPOINT,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as response:
            data = response.read().decode()
            print("OK", name, data[:600])
    except urllib.error.HTTPError as exc:
        print("ERR", name, exc.code, exc.read().decode("utf-8", "replace")[:800])


shapes = {
    "type_formula_match_subset": {"type": "formula", "value": "C F", "match_subset": True},
    "query_type_formula_match_subset": {"query_type": "formula", "value": "C F", "match_subset": True},
    "formula_upper": {"type": "formula", "value": "C1 F1", "match_subset": True},
    "descriptor_cf3": {
        "type": "descriptor",
        "descriptor_type": "SMILES",
        "match_type": "sub-struct-graph-relaxed",
        "value": "C(F)(F)F",
    },
    "query_type_descriptor_cf3": {
        "query_type": "descriptor",
        "descriptor_type": "SMILES",
        "match_type": "sub-struct-graph-relaxed",
        "value": "C(F)(F)F",
    },
}

for name, params in shapes.items():
    run(name, params)
