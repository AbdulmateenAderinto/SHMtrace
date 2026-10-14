from shmtrace import recal
from shmtrace.align import AlignedRead


def _read(name, quals, ref, errs=()):
    pairs = [(i, ("A" if ref[i] != "A" else "C") if i in errs else ref[i], q) for i, q in enumerate(quals)]
    return AlignedRead(name, "ref", 0, len(ref), pairs, "".join(b for _r, b, _q in pairs))


def test_bases_next_to_poor_ones_are_demoted_and_clean_ones_promoted():
    ref = "ACGT" * 25
    reads = []
    # Q20 bases: wrong 1 in 20 when a neighbour is Q5, 1 in 400 when every neighbour is Q30
    for k in range(400):
        noisy = [30] * 100
        for i in range(10, 90, 10):
            noisy[i - 1] = 5
            noisy[i] = 20
        reads.append(_read(f"n{k}", noisy, ref, errs={i for i in range(10, 90, 10)} if k % 20 == 0 else set()))
        clean = [30] * 100
        for i in range(15, 95, 10):
            clean[i] = 20
        reads.append(_read(f"c{k}", clean, ref, errs={i for i in range(15, 95, 10)} if k % 400 == 0 else set()))
    r = recal.fit([reads], ref, 0, 99, 0.5)
    t = r.table
    assert t[20][0] < 20 < t[20][2]             # demoted beside a Q5 base, promoted in a clean stretch
    assert all(t[40][nb] == 40 for nb in range(recal.N_NB))   # no evidence at Q35+: unchanged
    q = r.quals(reads[0].pairs)
    assert q[10] == t[20][0] and q[55] == t[30][2]          # a Q30 base among Q30s
