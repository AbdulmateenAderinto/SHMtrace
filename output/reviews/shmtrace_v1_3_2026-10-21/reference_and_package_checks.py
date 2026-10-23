"""Focused structural checks; no manuscript or package mutations."""

# ruff: noqa: E402
import collections
import json
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

HERE = Path(__file__).parent
ROOT = HERE.resolve().parents[2]
CRM = ROOT / "docs/calibration/crm"
sys.path.insert(0, str(CRM / "manuscript"))
sys.path.insert(0, str(CRM))
import build_crm as B
import lab_details
import make_review_package as PKG
from docbuild import MANUAL, REFS

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
with zipfile.ZipFile(CRM / "build/SHMtrace_CellReportsMethods_complete_v1.3.docx") as z:
    root = ET.fromstring(z.read("word/document.xml"))
    names = [b.get(W + "name") for b in root.iter(W + "bookmarkStart")]
    links = [h.get(W + "anchor") for h in root.iter(W + "hyperlink") if h.get(W + "anchor")]
    broken = sorted(set(links) - set(names))
    duplicate_names = [k for k, n in collections.Counter(names).items() if n > 1]
sources = "\n".join(
    (CRM / "manuscript" / f).read_text() for f in ["crm_text.py", "crm_methods.py", "build_crm.py"]
)
keys = {
    k.strip().lstrip("@")
    for match in re.findall(r"\[@([^\]]+)\]", sources)
    for k in match.split(",")
}
unknown = sorted(keys - set(MANUAL) - set(REFS))
refs = json.loads((HERE / "references.json").read_text())
doi_groups = collections.defaultdict(list)
for r in refs:
    doi = re.search(r"10\.\d{4,9}/\S+", r["citation"] or "")
    if doi:
        doi_groups[doi.group().rstrip(".").lower()].append(r["key"])
duplicates = {k: v for k, v in doi_groups.items() if len(v) > 1}
missing_package_images = [
    x + ".png" for x in PKG.FIGURES if not (CRM / "figures" / f"{x}.png").exists()
]
supp_imgs = {
    Path(p).name
    for _, paths, _ in B.SUPP_FIGS
    for p in (paths if isinstance(paths, list) else [paths])
}
packaged = {x + ".png" for x in PKG.FIGURES}
facts = {
    "cited_source_keys": len(keys),
    "unknown_reference_keys": unknown,
    "article_reference_entries": len(refs),
    "duplicate_article_dois": duplicates,
    "complete_bookmarks": len(names),
    "complete_internal_links": len(links),
    "broken_internal_links": broken,
    "duplicate_bookmark_names": duplicate_names,
    "missing_package_images": missing_package_images,
    "supplement_images_not_listed_in_package": sorted(supp_imgs - packaged),
    "laboratory_blockers": lab_details.missing(),
}
(HERE / "reference_and_package_checks.json").write_text(json.dumps(facts, indent=2))
print(json.dumps(facts, indent=2))
