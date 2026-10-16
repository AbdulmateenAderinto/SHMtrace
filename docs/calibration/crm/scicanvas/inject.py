"""Write a figure's numbers (DATA) into its SciCanvas builder, between /*DATA*/.../*END*/ markers.

The builders run sandboxed and read no files, so the numbers go into the builder itself. The shared plot helpers
(plots_common.py) are not injected: build.js saves them as the figure's plot prelude.
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def inject(builder: Path, data) -> int:
    src = builder.read_text()
    new = re.sub(r"/\*DATA\*/.*?/\*END\*/", lambda _m: "/*DATA*/" + json.dumps(data, separators=(",", ":")) + "/*END*/",
                 src, flags=re.S)
    builder.write_text(new)
    return len(new) - len(src)
