from pathlib import Path
import sys


NAME = sys.argv[1] if len(sys.argv) > 1 else "CF3"

rc_path = Path("/home/qin/.pymolrc.py")
exec(compile(rc_path.read_text(), str(rc_path), "exec"), globals(), globals())

build_metal_name_scenes(NAME)
