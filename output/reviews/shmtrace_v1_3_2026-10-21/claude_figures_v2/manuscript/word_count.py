"""Journal word count from recorded document operations, using the project's own word_count().

Mirrors build_crm.journal_word_count: everything in the article before STAR Methods counts (title through the
supplemental-information index) plus main-figure titles and legends; references, STAR Methods and Document S1 do
not. Calibrated against the Codex 1 count (6,823)."""
import json
import sys
from pathlib import Path

REVIEW = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REVIEW / "codex_version/sources"))
from docbuild import word_count  # noqa: E402


def count(ops):
    total, parts, sec = 0, {}, "front"
    for name, args, _kw in ops[3:]:                      # ops 0-2: version title page
        if name == "heading" and args[1] == 1:
            sec = args[0]
            if sec.startswith(("STAR", "Document S1", "References")):
                break
            continue
        if name == "heading" and sec in ("Results", "Resource availability"):
            n = word_count([("p", args[0])])
        elif name == "para":
            n = word_count([("p", args[0])])
        elif name == "caption":
            n = word_count([("p", args[0]), ("p", args[1])])
            sec_key = "legends"
            parts[sec_key] = parts.get(sec_key, 0) + n
            total += n
            continue
        else:
            continue
        parts[sec] = parts.get(sec, 0) + n
        total += n
    return total, parts


if __name__ == "__main__":
    for p in sys.argv[1:]:
        t, parts = count(json.loads(Path(p).read_text()))
        print(p, t, parts)
