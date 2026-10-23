"""Read-only manuscript audit; writes evidence into this review folder only."""

# ruff: noqa: E402
import difflib
import json
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
CRM = ROOT / "docs/calibration/crm"
HERE = Path(__file__).parent
sys.path.insert(0, str(CRM / "manuscript"))
import build_crm as B
from docbuild import MANUAL, REFS, Manuscript, word_count

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def accepted(el):
    if el.tag == W + "del":
        return ""
    if el.tag == W + "t":
        return el.text or ""
    if el.tag == W + "tab":
        return "\t"
    if el.tag == W + "br":
        return "\n"
    return "".join(accepted(ch) for ch in el)


def docx_text(path):
    with zipfile.ZipFile(path) as z:
        xml = ET.fromstring(z.read("word/document.xml"))
        paragraphs = [accepted(p) for p in xml.find(W + "body").findall(W + "p")]
        comments = (
            z.read("word/comments.xml").decode() if "word/comments.xml" in z.namelist() else ""
        )
        return paragraphs, {
            "insertions": len(xml.findall(".//" + W + "ins")),
            "deletions": len(xml.findall(".//" + W + "del")),
            "comments_xml": comments,
        }


n = {**B.crm_numbers.compute(), **B.public_numbers(), **B.dataset_counts()}
raw = {
    "INTRO": B.T.INTRO,
    "RESULTS": B.T.RESULTS,
    "DISCUSSION": B.T.DISCUSSION,
    "LIMITATIONS": B.T.LIMITATIONS,
    "RESOURCE": B.T.RESOURCE,
    "LEGENDS": B.M.LEGENDS,
}
(HERE / "source_blocks.json").write_text(json.dumps(raw, indent=2, ensure_ascii=False))
B.fill_modules(n)
count, parts = B.journal_word_count()
m = Manuscript()
B.write_article(m, n)
m.references()
fresh = [p.text for p in m.doc.paragraphs]
on_disk, changes = docx_text(CRM / "build/SHMtrace_CellReportsMethods.docx")
combined, combined_changes = docx_text(CRM / "build/SHMtrace_CellReportsMethods_complete_v1.3.docx")


def norm(ps):
    return [re.sub(r"\s+", " ", p).strip() for p in ps if p.strip()]


diff = "\n".join(
    difflib.unified_diff(
        norm(fresh),
        norm(on_disk),
        fromfile="current source",
        tofile="current article DOCX",
        lineterm="",
    )
)
(HERE / "source_vs_article.diff").write_text(diff)
(HERE / "accepted_complete.txt").write_text("\n\n".join(combined))
article_before_refs = norm(fresh[: next(i for i, p in enumerate(fresh) if p == "References")])
normalized_combined = norm(combined)
article_start = normalized_combined.index(article_before_refs[0])
combined_prefix = normalized_combined[article_start : article_start + len(article_before_refs)]
prefix_diff = "\n".join(
    difflib.unified_diff(
        article_before_refs,
        combined_prefix,
        fromfile="source article before refs",
        tofile="complete DOCX article prefix",
        lineterm="",
    )
)
(HERE / "source_vs_complete.diff").write_text(prefix_diff)
blocks = {
    "INTRO": B.T.INTRO,
    "RESULTS": B.T.RESULTS,
    "DISCUSSION": B.T.DISCUSSION,
    "LIMITATIONS": B.T.LIMITATIONS,
    "RESOURCE": B.T.RESOURCE,
    "LEGENDS": B.M.LEGENDS,
}
(HERE / "filled_blocks.json").write_text(json.dumps(blocks, indent=2, ensure_ascii=False))
(HERE / "numbers.json").write_text(json.dumps(n, indent=2, ensure_ascii=False, default=str))
lines = []
for name, value in blocks.items():
    lines.append("# " + name)
    for i, item in enumerate(value):
        txt = item if isinstance(item, str) else item[1]
        lines.append(f"[{i}; {word_count([('p', txt)])} words] {txt}")
(HERE / "filled_blocks.txt").write_text("\n\n".join(lines))
refs = [
    {
        "number": i,
        "key": key,
        "citation": MANUAL.get(key, REFS.get(key, {}).get("cell")),
        "metadata": REFS.get(key, {}),
    }
    for i, key in enumerate(m.order, 1)
]
(HERE / "references.json").write_text(json.dumps(refs, indent=2, ensure_ascii=False))
hardcoded = []
for name, value in vars(B.T).items():
    if name.isupper():
        hardcoded += [
            (name, x) for x in re.findall(r".{0,65}\^\{\d+(?:[–,]\d+)*\}.{0,55}", str(value))
        ]
summary = {
    "journal_words": count,
    "parts": dict(parts),
    "references_article": len(m.order),
    "article_diff_lines": len(diff.splitlines()),
    "complete_prefix_diff_lines": len(prefix_diff.splitlines()),
    "article_changes": changes,
    "complete_changes": combined_changes,
    "hardcoded_superscripts": hardcoded,
    "dataset_counts": B.dataset_counts(),
}
(HERE / "audit_snapshot.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
print(
    json.dumps(
        {k: v for k, v in summary.items() if k not in ("article_changes", "complete_changes")},
        indent=2,
    )
)
print(
    "Tracked changes article/complete:",
    changes["insertions"],
    changes["deletions"],
    combined_changes["insertions"],
    combined_changes["deletions"],
)
print("Reference order:", [(r["number"], r["key"]) for r in refs])
