import json
import urllib.request


data = json.load(urllib.request.urlopen("https://data.rcsb.org/rest/v1/core/polymer_entity/1BS3/1", timeout=30))

def walk(obj, path=""):
    if isinstance(obj, dict):
        for key, value in obj.items():
            if any(token in key.lower() for token in ["uniprot", "embl", "accession", "database", "identifier"]):
                print(path + "/" + key, json.dumps(value)[:1000])
            walk(value, path + "/" + key)
    elif isinstance(obj, list):
        for i, value in enumerate(obj[:5]):
            walk(value, path + f"[{i}]")

walk(data)
