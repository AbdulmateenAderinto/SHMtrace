"""Write a figure's numbers (DATA) and the shared code-plot helpers (HELPERS) into its SciCanvas builder.

The builders run sandboxed and read no files, and their code plots must re-run on their own inside SciCanvas, so
both go into the builder itself, between /*DATA*/.../*END*/ and /*HELPERS*/.../*END*/ markers.
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def inject(builder: Path, data) -> int:
    src = builder.read_text()
    new = re.sub(r"/\*DATA\*/.*?/\*END\*/", lambda _m: "/*DATA*/" + json.dumps(data, separators=(",", ":")) + "/*END*/",
                 src, flags=re.S)
    helpers = (HERE / "plots_common.py").read_text()
    new = re.sub(r"/\*HELPERS\*/.*?/\*END\*/", lambda _m: "/*HELPERS*/" + json.dumps(helpers) + "/*END*/", new, flags=re.S)
    builder.write_text(new)
    return len(new) - len(src)
