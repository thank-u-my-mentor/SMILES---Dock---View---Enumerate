import importlib.util
import sys

print(sys.version)
for name in ["pandas", "xlrd", "openpyxl", "pyxlsb", "olefile", "Bio"]:
    print(name, "OK" if importlib.util.find_spec(name) else "MISSING")
