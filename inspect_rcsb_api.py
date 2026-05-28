import json
import urllib.request


for url in [
    "https://data.rcsb.org/rest/v1/core/entry/1BS3",
    "https://data.rcsb.org/rest/v1/core/polymer_entity/1BS3/1",
]:
    print("\nURL", url)
    data = json.load(urllib.request.urlopen(url, timeout=30))
    print("top keys:", ", ".join(list(data)[:30]))
    print(json.dumps(data, indent=2)[:3500])
