"""SHMtrace against the Immcantation route (IgBLAST, then SHazaM) on the same N114 reads.

    python crm/immcantation_compare.py IGDATA       -> crm/out/immcantation_compare.json

IGDATA is IgBLAST's data folder (internal_data/mouse and optional_file/mouse_gl.aux); the IMGT mouse IGHV/IGHD/IGHJ
germlines are read from ~/immcantation_dbs/imgt_mouse/vdj. MiXCR was not run (its licence was not available).

The reads: the 500 reads SHMtrace analysed from each N114 library (`project.run_project`, three pMX-Aicda and
three empty-vector subclones, the empty-vector libraries as controls). Each is passed to IgBLAST (AIRR output) and
SHazaM's observedMutations counts substitutions over the V segment (v_sequence_alignment against
v_germline_alignment, frequency per read; the library value is the mean over reads), three ways:
  imgt              the IMGT germline IgBLAST assigns: what a default run reports
  parental          the parental N114 V segment as the germline, as a careful user of a cell line would supply
  parental_q20      the same, with every base below Q20 masked to N first (pRESTO FilterSeq maskqual's rule)
SHMtrace's frequency over the same V segment: supporting reads of its calls there, over analysed reads x V length.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from datasets import DATASETS  # noqa: E402

from shmtrace import project  # noqa: E402
from shmtrace.io_utils import read_fasta, read_fastq  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
IGB = Path("/usr/local/ncbi/igblast/bin")
IMGT = Path.home() / "immcantation_dbs/imgt_mouse/vdj"
COMP = str.maketrans("ACGTN", "TGCAN")
R_SCRIPT = r"""
suppressMessages(library(shazam))
args <- commandArgs(trailingOnly = TRUE)
db <- read.delim(args[1], stringsAsFactors = FALSE)
db <- db[!is.na(db$v_sequence_alignment) & nchar(db$v_sequence_alignment) > 0, ]
db <- observedMutations(db, sequenceColumn = "v_sequence_alignment", germlineColumn = "v_germline_alignment",
                        regionDefinition = NULL, frequency = TRUE, combine = TRUE, nproc = 1)
write.table(db[, c("sequence_id", "mu_freq")], args[2], sep = "\t", row.names = FALSE, quote = FALSE)
"""


def run(cmd, **kw):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, **kw)


def main(igdata):
    ds = DATASETS["n114"]
    work = Path(tempfile.mkdtemp(prefix="shmtrace_immc_"))
    env = {**os.environ, "IGDATA": igdata}
    aux = str(Path(igdata) / "optional_file/mouse_gl.aux")
    ref = read_fasta(ds["reference"])[0][1]

    # SHMtrace on the six libraries, controls the empty-vector subclones
    proj = {"name": "n114", "region": "n114_vregion", "reference": ds["reference"], "mode": "recommended",
            "depth": ds["depth"],
            "samples": [{"name": k, "path": s["fastq"], "group": "empty" if s["control"] else "AID",
                         "control": s["control"]} for k, s in ds["samples"].items()]}
    out = project.run_project(proj)

    # germline databases: IMGT, and IMGT with the parental V segment in place of the V genes
    def db(name, fasta_text):
        fa = work / f"{name}.fa"
        fa.write_text(fasta_text)
        run([str(IGB / "makeblastdb"), "-parse_seqids", "-dbtype", "nucl", "-in", str(fa), "-out", str(work / name)])
    for g in "VDJ":
        clean = subprocess.run(["perl", str(IGB / "edit_imgt_file.pl"), str(IMGT / f"IMGT_mouse_IGH{g}.fasta")],
                               capture_output=True, text=True, check=True).stdout
        db(f"IGH{g}", clean)
    # where IgBLAST puts the V segment on the parental reference
    (work / "ref.fa").write_text(f">ref\n{ref}\n")
    run([str(IGB / "igblastn"), "-germline_db_V", str(work / "IGHV"), "-germline_db_D", str(work / "IGHD"),
         "-germline_db_J", str(work / "IGHJ"), "-organism", "mouse", "-ig_seqtype", "Ig", "-auxiliary_data", aux,
         "-query", str(work / "ref.fa"), "-outfmt", "19", "-out", str(work / "ref.tsv")], env=env)
    r0 = next(csv.DictReader(open(work / "ref.tsv"), delimiter="\t"))
    sense = ref.translate(COMP)[::-1] if r0["rev_comp"] == "T" else ref
    vs, ve = int(r0["v_sequence_start"]), int(r0["v_sequence_end"])
    db("PARV", f">{r0['v_call'].split(',')[0]}_N114parental\n{sense[vs - 1:ve]}\n")
    # the V segment in forward reference coordinates (1-based, inclusive)
    v_lo, v_hi = (len(ref) - ve + 1, len(ref) - vs + 1) if r0["rev_comp"] == "T" else (vs, ve)

    rows = []
    for name, s in ds["samples"].items():
        res = out["_details"][name]
        keep = {a.qname for a in res["_aligned"]}
        reads = [r for r in read_fastq(s["fastq"]) if r[0].split()[0] in keep]
        sh_mut = sum(v.supporting_reads for v in res["variants"] if v_lo <= v.pos <= v_hi)
        row = {"library": name, "group": "empty vector" if s["control"] else "pMX-Aicda", "reads": len(reads),
               "shmtrace": sh_mut / (len(reads) * (v_hi - v_lo + 1)),
               "shmtrace_calls_in_v": sum(1 for v in res["variants"] if v_lo <= v.pos <= v_hi)}
        for variant, vdb, mask in (("imgt", "IGHV", False), ("parental", "PARV", False),
                                   ("parental_q20", "PARV", True)):
            fa = work / f"{name}_{variant}.fa"
            with open(fa, "w") as fh:
                for n, seq, q in reads:
                    if mask:
                        seq = "".join(b if ord(c) - 33 >= 20 else "N" for b, c in zip(seq, q, strict=True))
                    fh.write(f">{n.split()[0]}\n{seq}\n")
            tsv = work / f"{name}_{variant}.tsv"
            run([str(IGB / "igblastn"), "-germline_db_V", str(work / vdb), "-germline_db_D", str(work / "IGHD"),
                 "-germline_db_J", str(work / "IGHJ"), "-organism", "mouse", "-ig_seqtype", "Ig",
                 "-auxiliary_data", aux, "-query", str(fa), "-outfmt", "19", "-out", str(tsv),
                 "-num_threads", "4"], env=env)
            mu = work / f"{name}_{variant}_mu.tsv"
            run(["Rscript", "-e", R_SCRIPT, str(tsv), str(mu)])
            vals = [float(r["mu_freq"]) for r in csv.DictReader(open(mu), delimiter="\t") if r["mu_freq"] != "NA"]
            row[variant] = sum(vals) / len(vals)
            row[f"{variant}_reads"] = len(vals)
        rows.append(row)
        print(name, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    OUT.mkdir(exist_ok=True)
    json.dump({"v_call": r0["v_call"], "v_identity_parental_to_imgt": float(r0["v_identity"]),
               "v_region_reference": [v_lo, v_hi], "rows": rows},
              open(OUT / "immcantation_compare.json", "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
