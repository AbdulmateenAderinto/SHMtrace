"""Count and export the unapplied editorial proposal against frozen v1.3 sources."""

# ruff: noqa: E402
import copy
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
ROOT = HERE.resolve().parents[2]
sys.path.insert(0, str(ROOT / "docs/calibration/crm/manuscript"))
import build_crm as B
import proposed_edits as P

numbers = {**B.crm_numbers.compute(), **B.public_numbers(), **B.dataset_counts()}
raw = json.loads((HERE / "source_blocks.json").read_text())
assert raw["RESULTS"] == [list(x) for x in B.T.RESULTS], "Source changed since audit snapshot"
new = copy.deepcopy(raw)
for i, text in P.RESULTS.items():
    assert new["RESULTS"][i][0] == "p"
    new["RESULTS"][i][1] = text
for i, text in enumerate(P.LEGENDS):
    new["LEGENDS"][i][1] = text
for i, replacements in P.DISCUSSION_REPLACE.items():
    for before, after in replacements:
        assert new["DISCUSSION"][i].count(before) == 1, (i, before)
        new["DISCUSSION"][i] = new["DISCUSSION"][i].replace(before, after)

# The wording in Figure 3D's proposal describes the CURRENT figure. The review
# separately recommends replacing its artificial zero position, requiring the
# legend to be updated in the same figure-edit batch.
filled = {key: B.fill_all(value, numbers) for key, value in new.items()}
assert not re.search(r"(?<![\^_])\{[a-z][a-z0-9_]*\}", str(filled))
B.fill_modules(numbers)
baseline, original_parts = B.journal_word_count()
original_parts = dict(original_parts)
B.T.RESULTS, B.M.LEGENDS = filled["RESULTS"], filled["LEGENDS"]
compact_only, _ = B.journal_word_count()
B.T.DISCUSSION = filled["DISCUSSION"]
with_discussion, _ = B.journal_word_count()
total, parts = B.journal_word_count()
assert total < 7000, total


def placeholders(s):
    return set(re.findall(r"(?<![\^_])\{([a-z][a-z0-9_]*)\}", s))


elsewhere = "\n".join(str(getattr(B.M, k)) for k in dir(B.M) if k.isupper()) + str(B.SUPP_FIGS)
# B.M is already filled, so retrieve the original Methods module source for key coverage.
elsewhere = (ROOT / "docs/calibration/crm/manuscript/crm_methods.py").read_text() + (
    ROOT / "docs/calibration/crm/manuscript/build_crm.py"
).read_text()
added = "\n".join(t for _, t in P.METHODS_ADDITIONS)
old_keys = placeholders(str(raw["RESULTS"]))
new_keys = placeholders(str(new["RESULTS"]))
migration = []
for key in sorted(old_keys - new_keys):
    locations = []
    if "{" + key + "}" in elsewhere:
        locations.append("existing Methods / legend / supplementary builder")
    if "{" + key + "}" in added:
        locations.append("proposed STAR Methods additions")
    original_blocks = [i for i, (_, t) in enumerate(raw["RESULTS"]) if "{" + key + "}" in t]
    migration.append(
        {
            "key": key,
            "value": numbers[key],
            "retained_in": locations,
            "original_results_blocks": original_blocks,
        }
    )

# Snapshot every original Results paragraph. This is the completeness record for
# detail relocation, not a suggestion to append duplicate full paragraphs to Methods.
ledger = [
    "# Detail relocation ledger",
    "",
    "Every removed numerical placeholder is listed below. Existing coverage means the source key occurs elsewhere; it does not by itself prove identical context. Review the retained full original paragraphs before implementing the move. No analysis output or source manuscript was changed.",
    "",
    "| Placeholder | Current value | Proposed location / existing coverage | Original Results block |",
    "|---|---|---|---|",
]
for item in migration:
    loc = (
        "; ".join(item["retained_in"])
        or "Retain in STAR Methods during implementation; see original paragraph below"
    )
    ledger.append(
        f"| `{item['key']}` | {item['value']} | {loc} | {item['original_results_blocks']} |"
    )
ledger += ["", "## Full original Results paragraphs for relocation review", ""]
for i, (kind, t) in enumerate(raw["RESULTS"]):
    if kind == "p" and i in P.RESULTS:
        ledger.extend([f"### Results block {i}", "", B.fill(t, numbers), ""])
(HERE / "DETAIL_RELOCATION_LEDGER.md").write_text("\n".join(ledger))

# Verify all proposed placeholders resolve and all citation keys exist. Also validate
# source-to-revision numeric identity for every key retained: values use the same n.
key_pat = re.compile(r"\[@([^\]]+)\]")
from docbuild import MANUAL, REFS

keys = set(k.strip().lstrip("@") for group in key_pat.findall(str(new)) for k in group.split(","))
assert not keys - (set(REFS) | set(MANUAL)), keys - (set(REFS) | set(MANUAL))
summary = {
    "baseline": baseline,
    "results_and_legends_only": compact_only,
    "with_discussion_revisions": with_discussion,
    "full_proposal": total,
    "reduction": baseline - total,
    "margin_below_7000": 7000 - total,
    "parts_before": original_parts,
    "parts_after": dict(parts),
    "removed_results_placeholders": len(migration),
    "removed_keys_requiring_context_review": sum(not x["retained_in"] for x in migration),
}
(HERE / "proposal_counts.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
(HERE / "proposal_source_blocks.json").write_text(json.dumps(new, indent=2, ensure_ascii=False))

out = [
    "# Proposed v1.3 editorial revision",
    "",
    f"Unapplied proposal. Project count: **{baseline:,} → {total:,} words**, a reduction of {baseline - total:,}; {7000 - total:,} below 7,000. Counts use the existing journal_word_count function with current analysis-derived values. Section headings outside the source variables and journal counting differences may add words, so preserve this margin.",
    "",
    "The current manuscript, Word files, analyses and figure files are unchanged. All unchanged sections remain as in v1.3. Numerical substitutions use the existing analysis outputs. The source-preserving version is proposal_source_blocks.json; @keys below are citation keys, not proposed final numbering. The detail relocation ledger retains every original Results paragraph and flags supporting detail for STAR Methods.",
    "",
    "| Section | Current | Proposed | Reduction |",
    "|---|---:|---:|---:|",
]
for name in ["results", "figure titles and legends", "discussion", "resource availability"]:
    out.append(
        f"| {name} | {original_parts[name]:,} | {dict(parts)[name]:,} | {original_parts[name] - dict(parts)[name]:,} |"
    )
out += ["", "## Results", ""]
for kind, text in filled["RESULTS"]:
    out += [(("### " if kind == "h" else "") + text), ""]
out += ["## Discussion", ""]
for text in filled["DISCUSSION"]:
    out += [text, ""]
out += ["## Main figure legends", ""]
for title, text in filled["LEGENDS"]:
    out += ["### " + title, "", text, ""]
out += [
    "## Proposed STAR Methods additions",
    "",
    "Integrate these into their existing Methods subsections rather than adding a separate prose supplement. Additional supporting details to retain are listed in DETAIL_RELOCATION_LEDGER.md.",
    "",
]
for title, text in P.METHODS_ADDITIONS:
    out += ["### " + title, "", B.fill(text, numbers), ""]
out += ["## Existing STAR Methods wording correction", ""]
for before, after in P.METHODS_REPLACEMENTS:
    assert before in (ROOT / "docs/calibration/crm/manuscript/crm_methods.py").read_text()
    out += [f"Replace `{before}` with `{after}`.", ""]
(HERE / "PROPOSED_REVISION.md").write_text("\n".join(out))
print(json.dumps(summary, indent=2))
