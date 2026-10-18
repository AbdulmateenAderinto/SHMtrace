"""Package the revision for the co-author who reviewed it on 17 Oct 2026.

    python docs/calibration/crm/make_coauthor_package.py [their_review.zip]
        -> crm/build/SHMtrace_coauthor_package_<date>.zip

Unlike make_review_package.py, which is a cold package for an outside reviewer, this one answers a specific
review. Its cover note (coauthor_response.md) has placeholders that are filled from evidence computed here, so
the note cannot claim a number the package does not show:

  - the documents are rebuilt and the build log parsed for word counts and blockers;
  - the test suite is run and its result recorded;
  - the regenerated outputs are compared field by field with the reviewer's own, read straight from their zip
    (JSON and the Table S2 workbook only; nothing in it is executed);
  - source diffs are taken from git against their snapshot and against their patch.

It refuses to run with tracked changes uncommitted, because the note names the commit it describes.
"""
import datetime
import json
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

CRM = Path(__file__).resolve().parent
CAL = CRM.parent
REPO = CAL.parents[1]
BUILD, FIGS, OUT = CRM / "build", CRM / "figures", CRM / "out"
PY = sys.executable
SOFFICE = "/Applications/LibreOffice.app/Contents/MacOS/soffice"
THEIR_ZIP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/nanoshm/output/reviews/"
                           "SHMtrace_revision_2026-10-17/SHMtrace_review_and_edits.zip")
THEIR_SNAPSHOT, THEIR_PATCH = "83a908c", "3317fd9"
THEIR_OUT = "revised_data/docs/calibration/crm/out/"
DATA = ["fdr_validation", "fdr_null_reps", "fdr_estimator", "benchmark_real_jh4", "benchmark_real_n114",
        "figS15_numbers", "figS15_estimator", "fig2_numbers"]
MAIN_FIGS = ["GraphicalAbstract_preview", "Figure1_preview", "Figure2", "Figure3", "Figure4", "Figure5",
             "Figure6", "Figure7"]
NOT_SOURCE = [f":(exclude)*.{x}" for x in ("png", "pdf", "tif", "svg", "scifig", "docx", "xlsx", "zip")]


def sh(*cmd, cwd=REPO, check=True):
    return subprocess.run(cmd, cwd=cwd, check=check, capture_output=True, text=True)


def build_documents(evidence):
    """Rebuild, so the package cannot ship documents older than the source, and parse what the build reports."""
    log = sh(PY, "-W", "ignore", "crm/manuscript/build_crm.py", "crm/build", cwd=CAL).stdout
    (evidence / "manuscript_build.log").write_text(log)
    m = re.search(r"summary (\d+) words; main text \(introduction to limitations\) (\d+) words", log)
    h = re.search(r"highlights: \[([\d, ]+)\] eTOC words: (\d+)", log)
    b = re.search(r"submission blockers: (\d+)", log)
    if not (m and h):
        sys.exit("could not read word counts from the build log")
    return {"summary_words": m.group(1), "main_words": f"{int(m.group(2)):,}", "etoc_words": h.group(2),
            "highlights": h.group(1).replace(" ", ""), "blockers": b.group(1) if b else "0"}


def run_tests(evidence):
    """The project config adds -q, so pytest prints progress dots and no summary line: count them."""
    r = sh(PY, "-m", "pytest", "-p", "no:cacheprovider", "-W", "ignore", check=False)
    progress = "".join(line.split("[")[0] for line in r.stdout.splitlines() if re.match(r"^[.sFEx]+\s+\[", line))
    res = {"exit_code": r.returncode, "passed": progress.count("."), "failed": progress.count("F"),
           "errors": progress.count("E"), "skipped": progress.count("s")}
    (evidence / "tests.json").write_text(json.dumps(res, indent=1))
    if r.returncode != 0:
        (evidence / "tests_output.txt").write_text(r.stdout + r.stderr)
    return res


def compare_with_reviewer(evidence):
    """Field-by-field comparison of the regenerated outputs with the reviewer's, read from their zip."""
    if not THEIR_ZIP.exists():
        note = {"available": False, "reason": f"{THEIR_ZIP} not found"}
        (evidence / "comparison_with_reviewer.json").write_text(json.dumps(note, indent=1))
        return note
    z = zipfile.ZipFile(THEIR_ZIP)
    key = lambda r: (r["signal"], r["controls"], r["noise"], r["exclude_noisy"], r["depth"], r["seed"],  # noqa: E731
                     r.get("estimate", "mean"))
    fields = ("discoveries", "true_calls", "false_calls", "fdp", "any_discovery", "planted", "decoys_set_aside")
    fdr = {"total": 0, "identical": 0, "order_only": 0, "differing": []}
    for name in ("fdr_validation", "fdr_null_reps", "fdr_estimator"):
        mine = json.loads((OUT / f"{name}.json").read_text())["runs"]
        theirs = {key(r): r for r in json.loads(z.read(f"{THEIR_OUT}{name}.json"))["runs"]}
        for r in mine:
            fdr["total"] += 1
            o = theirs.get(key(r))
            if o is not None and all(r[f] == o[f] for f in fields):
                fdr["identical"] += 1
                if r["control_names"] != o["control_names"] and sorted(r["control_names"]) == sorted(o["control_names"]):
                    fdr["order_only"] += 1
            else:
                fdr["differing"].append([name, list(key(r))])
    bench = {"total": 0, "identical": 0, "differing": []}
    for name in ("benchmark_real_jh4", "benchmark_real_n114"):
        rows = lambda d: {(r["sample"], r["mode"]): (r["calls"], r["bases"]) for r in d["rows"]}  # noqa: E731
        mine, theirs = rows(json.loads((OUT / f"{name}.json").read_text())), rows(json.loads(z.read(f"{THEIR_OUT}{name}.json")))
        for k, v in mine.items():
            bench["total"] += 1
            if theirs.get(k) == v:
                bench["identical"] += 1
            else:
                bench["differing"].append([name, list(k)])
    s2 = compare_table_s2(z)
    res = {"available": True, "source": str(THEIR_ZIP), "fdr": fdr, "benchmark": bench, "table_s2": s2}
    (evidence / "comparison_with_reviewer.json").write_text(json.dumps(res, indent=1))
    return res


def compare_table_s2(z):
    import io

    import openpyxl
    mine = openpyxl.load_workbook(OUT / "Table_S2_benchmark.xlsx", data_only=True)
    theirs = openpyxl.load_workbook(io.BytesIO(z.read("deliverables/Table_S2_benchmark.xlsx")), data_only=True)
    out = {}
    for sheet in mine.sheetnames:
        a = [list(r) for r in mine[sheet].iter_rows(values_only=True)]
        b = [list(r) for r in theirs[sheet].iter_rows(values_only=True)] if sheet in theirs.sheetnames else None
        out[sheet] = None if b is None else sum(x != y for x, y in zip(a, b, strict=False)) + abs(len(a) - len(b))
    return out


def fill_note(values, cmp_):
    note = (CRM / "coauthor_response.md").read_text()
    if cmp_.get("available"):
        f, b, s2 = cmp_["fdr"], cmp_["benchmark"], cmp_["table_s2"]
        values.update(fdr_same=f["identical"], fdr_total=f["total"], bench_same=b["identical"], bench_total=b["total"])
        values["s2_summary"] = ("identical, sheet for sheet" if all(v == 0 for v in s2.values())
                                else "differs: " + ", ".join(f"{k} {v} rows" for k, v in s2.items() if v))
        values["fdr_order_note"] = (
            f"{f['order_only']} of the identical validation runs list the same three controls in a different order. "
            "Your package kept the unmodified-background cells from before the fix, generated with the old sorted "
            "selection; here every cell was rerun with the balanced selection, which lists them differently. The "
            "results are the same." if f["order_only"] else "")
    else:
        values.update(fdr_same="n/a", fdr_total="n/a", bench_same="n/a", bench_total="n/a",
                      s2_summary="not compared (reviewer package not found)", fdr_order_note="")
    missing = sorted(set(re.findall(r"\{([a-z_0-9]+)\}", note)) - set(values))
    if missing:
        sys.exit(f"cover note placeholders without a value: {missing}")
    return note.format(**values)


def main():
    if sh("git", "status", "--porcelain", "--untracked-files=no").stdout.strip():
        sys.exit("commit tracked changes first: the cover note names the commit it describes")
    head = sh("git", "rev-parse", "--short", "HEAD").stdout.strip()
    stamp = datetime.date.today().isoformat()
    pkg = BUILD / f"SHMtrace_coauthor_package_{stamp}"
    if pkg.exists():
        shutil.rmtree(pkg)
    for d in ("figures/supplement", "patches", "data", "evidence"):
        (pkg / d).mkdir(parents=True)
    evidence = pkg / "evidence"

    print("rebuilding documents ...", flush=True)
    values = {"head": head, **build_documents(evidence)}
    print("running tests ...", flush=True)
    tests = run_tests(evidence)
    values["tests_passed"] = tests["passed"]
    if tests["exit_code"] != 0:
        print(f"WARNING: tests exited {tests['exit_code']}; see evidence/tests_output.txt")
    print("comparing with the reviewer's outputs ...", flush=True)
    cmp_ = compare_with_reviewer(evidence)

    # documents
    docs = {"SHMtrace_CellReportsMethods": "01_Manuscript", "Document_S1": "02_Supplement_Document_S1",
            "Highlights_eTOC_Motivation": "03_Highlights_eTOC_Motivation"}
    sh(SOFFICE, "--headless", "--convert-to", "pdf", *[f"{d}.docx" for d in docs], "--outdir", str(BUILD), cwd=BUILD)
    for src, name in docs.items():
        for ext in ("docx", "pdf"):
            shutil.copy(BUILD / f"{src}.{ext}", pkg / f"{name}.{ext}")
        if name != "03_Highlights_eTOC_Motivation":
            sh("pdftotext", "-layout", str(pkg / f"{name}.pdf"), str(pkg / f"{name}.txt"))
    shutil.copy(OUT / "Table_S2_benchmark.xlsx", pkg / "04_Table_S2_benchmark.xlsx")
    shutil.copy(CRM / "REVIEW_RESPONSE.md", pkg / "05_REVIEW_RESPONSE.md")
    shutil.copy(CRM / "HANDOVER.md", pkg / "06_HANDOVER.md")
    shutil.copy(CRM / "REFEREE_RESPONSE.md", pkg / "07_REFEREE_RESPONSE.md")

    # figures: the main figures, then every image Document S1 embeds, taken from the build's own list so the set
    # cannot drift from the document
    for f in MAIN_FIGS:
        shutil.copy(FIGS / f"{f}.png", pkg / "figures" / f"{f.replace('_preview', '')}.png")
    for f in ("Figure1", "GraphicalAbstract"):
        shutil.copy(FIGS / f"{f}.scifig", pkg / "figures" / f"{f}.scifig")
    sys.path.insert(0, str(CRM / "manuscript"))
    import build_crm
    for title, img, _legend in build_crm.SUPP_FIGS:
        label = re.match(r"Figure (S\d+)", title).group(1)
        part = re.search(r"page (\d) of", title)
        for k, p in enumerate(img if isinstance(img, list) else [img]):
            suffix = (f"_p{part.group(1)}" if part else "") + (f"_{k + 1}" if isinstance(img, list) and len(img) > 1 else "")
            shutil.copy(p, pkg / "figures/supplement" / f"Figure{label}{suffix}_{Path(p).stem}.png")

    # source diffs, text only
    for name, base in (("changes_since_your_snapshot", THEIR_SNAPSHOT), ("changes_beyond_your_patch", THEIR_PATCH)):
        diff = sh("git", "diff", base, "HEAD", "--", ".", *NOT_SOURCE).stdout
        (pkg / "patches" / f"{name}.patch").write_text(diff)
    (pkg / "patches" / "commits.txt").write_text(sh("git", "log", "--oneline", f"{THEIR_SNAPSHOT}..HEAD").stdout)

    for d in DATA:
        shutil.copy(OUT / f"{d}.json", pkg / "data" / f"{d}.json")

    (pkg / "00_README.md").write_text(fill_note(values, cmp_))
    (evidence / "package.json").write_text(json.dumps({"head": head, "date": stamp, **values}, indent=1))

    zp = BUILD / f"{pkg.name}.zip"
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(pkg.rglob("*")):
            if f.is_file():
                z.write(f, f.relative_to(BUILD))
    print(f"{zp} ({zp.stat().st_size / 1e6:.1f} MB) at {head}")
    print(f"tests {tests['passed']} passed, exit {tests['exit_code']}; "
          f"FDR {values['fdr_same']}/{values['fdr_total']}, benchmark {values['bench_same']}/{values['bench_total']}, "
          f"Table S2 {values['s2_summary']}")


if __name__ == "__main__":
    main()
