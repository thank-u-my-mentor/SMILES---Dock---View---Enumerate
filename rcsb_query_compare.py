import json
import urllib.error
import urllib.request


ENDPOINT = "https://search.rcsb.org/rcsbsearch/v2/query"


def post(payload):
    payload = dict(payload)
    payload["request_options"] = dict(payload.get("request_options", {}))
    payload["request_options"]["paginate"] = {"start": 0, "rows": 1000}
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        ENDPOINT,
        data=data,
        headers={"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            body = r.read()
            if not body.strip():
                return 0, []
            out = json.loads(body.decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", "replace")
        print("HTTP_ERROR", exc.code, body[:1000])
        raise
    ids = [row["identifier"] for row in out.get("result_set", [])]
    return out.get("total_count", len(ids)), ids


def nested_instance(comp_node, type_value):
    return {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            comp_node,
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_nonpolymer_instance_annotation.type",
                    "operator": "exact_match",
                    "value": type_value,
                },
            },
        ],
        "label": "nested-attribute",
    }


def comp_exact(value):
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
            "operator": "exact_match",
            "value": value,
        },
    }


def comp_in(values):
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
            "operator": "in",
            "value": values,
        },
    }


def chem_comp_in(values):
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
            "operator": "in",
            "value": values,
        },
    }


metals = ["FE", "CO", "NI", "CU", "MN"]
f_no_cov = nested_instance(comp_exact("F"), "HAS_NO_COVALENT_LINKAGE")

base_opts = {
    "results_content_type": ["experimental"],
    "sort": [{"sort_by": "score", "direction": "desc"}],
    "scoring_strategy": "combined",
}

queries = {
    "current_chem_comp_metal": {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [f_no_cov, chem_comp_in(metals)],
        },
        "return_type": "entry",
        "request_options": base_opts,
    },
    "metal_instance_any_type": {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [f_no_cov, comp_in(metals)],
        },
        "return_type": "entry",
        "request_options": base_opts,
    },
    "metal_instance_has_metal_coordination": {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [f_no_cov, nested_instance(comp_in(metals), "HAS_METAL_COORDINATION")],
        },
        "return_type": "entry",
        "request_options": base_opts,
    },
    "metal_instance_no_covalent": {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [f_no_cov, nested_instance(comp_in(metals), "HAS_NO_COVALENT_LINKAGE")],
        },
        "return_type": "entry",
        "request_options": base_opts,
    },
}


for name, payload in queries.items():
    total, ids = post(payload)
    print(name, total)
    print(" ".join(ids))
