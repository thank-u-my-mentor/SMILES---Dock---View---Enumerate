import ast
import json
import re
import urllib.parse
import urllib.request


text = open("/mnt/e/Codex/build_metal_f_from_rcsb.py", encoding="utf-8").read()
match = re.search(r"RCSB_SEARCH_URL = (\".*\")", text)
url = ast.literal_eval(match.group(1))

payload = json.loads(urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["request"][0])
payload.pop("request_info", None)
payload.setdefault("request_options", {})
payload["request_options"]["paginate"] = {"start": 0, "rows": 1000}
payload["return_type"] = "entry"

request = urllib.request.Request(
    "https://search.rcsb.org/rcsbsearch/v2/query",
    data=json.dumps(payload).encode("utf-8"),
    headers={"Content-Type": "application/json"},
)
with urllib.request.urlopen(request, timeout=60) as response:
    result = json.load(response)

ids = [row["identifier"] for row in result.get("result_set", [])]
print(len(ids))
print(" ".join(ids))
