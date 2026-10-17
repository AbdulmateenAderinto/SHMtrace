"""Build the mouse Sμ and Sγ1 switch references used for the PacBio CSR analysis.

    python crm/build_switch_reference.py      -> crm/out/switch_ref/mouse_Smu_Sg1.fasta

The repository bundles only a 651-bp Sμ fragment (`mouse_igh_smu_ch12.fasta`), which is too short to carry a
class-switch junction: a junction read needs the donor region, the acceptor region and the primer sites at both
ends. This fetches the IgH locus from Ensembl (GRCm39, chromosome 12, minus strand, so the sequence runs in the
locus's own 5'->3' direction) and cuts the two windows out of it.

Two checks run before anything is written, because a fetch that silently returns the wrong window would be
invisible later: the project's own J_H4 intron reference must align inside the fetched locus, and so must the
bundled Sμ fragment. The window edges come from where the reads themselves start and end (crm/switch_pacbio.py
records this), not from a guess: the amplicon's forward primer sits just upstream of Sμ and its reverse primer
inside Cγ1, so the acceptor window has to reach past Sγ1 into Cγ1.
"""
import json
import sys
import urllib.request
from pathlib import Path

CRM = Path(__file__).resolve().parent
REPO = CRM.parents[2]
sys.path.insert(0, str(REPO / "src"))
import mappy as mp  # noqa: E402

REF = REPO / "src/shmtrace/references"
OUT = CRM / "out/switch_ref"
# the fetched locus, in chromosome coordinates; index i of the fetched string is coordinate TOP - i
LO, TOP = 113_280_000, 113_396_000
URL = f"https://rest.ensembl.org/sequence/region/mus_musculus/12:{LO}..{TOP}:-1?content-type=text/x-fasta"
# donor and acceptor windows, chosen to span the amplicon the reads show (see the module docstring)
WINDOWS = [("Smu", 113_386_000, 113_392_500, "donor: Sμ core plus the forward primer site upstream"),
           ("Sg1", 113_290_000, 113_304_000, "acceptor: Sγ1 through Cγ1, to the reverse primer site")]


def fasta(path):
    return "".join(x for x in Path(path).read_text().split("\n") if not x.startswith(">"))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    raw = OUT / "igh_locus_GRCm39.fa"
    if not raw.exists():
        with urllib.request.urlopen(URL, timeout=180) as fh:          # noqa: S310 (fixed https URL)
            raw.write_bytes(fh.read())
    header = raw.read_text().split("\n", 1)[0]
    if f"{LO}:{TOP}" not in header:
        sys.exit(f"Ensembl returned a different window: {header}")
    locus = fasta(raw)
    if len(locus) != TOP - LO + 1:
        sys.exit(f"locus length {len(locus)}, expected {TOP - LO + 1}")
    # the two checks: both must land, or the window is not what it claims to be
    al = mp.Aligner(seq=locus, preset="asm5")
    checks = {}
    for name in ("mouse_jh4_intron.fasta", "mouse_igh_smu_ch12.fasta"):
        q = fasta(REF / name)
        hit = max(al.map(q), key=lambda h: h.mlen, default=None)
        if hit is None or hit.mlen < 0.8 * len(q):
            sys.exit(f"{name} does not align inside the fetched locus; the window is wrong")
        checks[name] = {"locus_start": hit.r_st, "locus_end": hit.r_en, "strand": hit.strand,
                        "identity_pct": round(100 * hit.mlen / hit.blen, 2),
                        "chr12": [TOP - hit.r_en, TOP - hit.r_st]}
    records = []
    for name, lo, hi, note in WINDOWS:
        seq = locus[TOP - hi:TOP - lo]
        records.append((name, seq, f"chr12:{lo}-{hi} GRCm39 ({note})"))
    out = OUT / "mouse_Smu_Sg1.fasta"
    with open(out, "w") as fh:
        for name, seq, note in records:
            fh.write(f">{name} {note}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
    meta = {"source": URL, "assembly": "GRCm39", "header": header,
            "windows": {n: {"chr12": [lo, hi], "bp": hi - lo, "note": note} for n, lo, hi, note in WINDOWS},
            "checks": checks}
    json.dump(meta, open(OUT / "reference_provenance.json", "w"), indent=1)
    for name, seq, note in records:
        print(f"{name}: {len(seq)} bp  {note}")
    for k, v in checks.items():
        print(f"check {k}: {v['identity_pct']}% identity at chr12:{v['chr12'][0]}-{v['chr12'][1]}")


if __name__ == "__main__":
    main()
