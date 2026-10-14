# Base-quality recalibration from the controls (tested, off by default)

*14 October 2026. `shmtrace/recal.py`, `Params.recalibrate`.*

## Why it was tried

Variant callers such as Clair3 and ClairS learn that a base's error rate depends on its context, not only on its
reported quality. On the AID-null controls this is true. Error rate (× 10⁻⁴) of bases at the same reported quality,
split by the lowest quality among the two bases either side:

| Reported quality | Neighbours < Q10 | 10–19 | ≥ Q20 |
|---|---:|---:|---:|
| J_H4, Q15–24 | 90 | 73 | 22 |
| J_H4, Q25–34 | 22 | 14 | 5.9 |
| N114, Q15–24 | 112 | 63 | 14 |
| N114, Q25–34 | 7.4 | 3.7 | 1.5 |

Distance from the read end and homopolymer length had weak, inconsistent effects (homopolymers are a property
of the reference position, which the per-site control rates already capture).

## What it does

For each reported-quality band (Q0–14, 15–24, 25–34, 35+) and neighbourhood, the controls give a shift in dB
relative to the band as a whole; each base's quality is shifted by it. Shrinkage is 20 pseudo-errors at the band
rate. The table is fitted on the same controls as the null it belongs to, so no library is ever recalibrated
with its own errors. Fitted shifts (rows are bands, columns neighbourhoods < 10, 10–19, ≥ 20):

    J_H4  [[-1, 3, 3], [-1, -1, 3], [-1, -1, 1], [0, 0, 0]]
    N114  [[-1, 3, 2], [-3, -1, 4], [0, -1, 0], [0, 0, 0]]

A first version remapped each reported quality to its empirical rate; it collapsed the high qualities, which the
controls are too few to tell apart, and lowered one-read recall (J_H4 0.30 → 0.23, N114 0.58 → 0.46). It was
replaced by the shift.

## Result (planted mutations, `crm/spikein.py --tag _recal2`)

| | 1-read | 2-read | all | false calls / library | FDR |
|---|---:|---:|---:|---:|---:|
| J_H4 recommended | 0.30 | 0.76 | 0.809 | 2.44 | 0.048 |
| J_H4 recommended + recalibration | 0.30 | 0.78 | 0.819 | 2.56 | 0.049 |
| J_H4 sensitive | 0.72 | 0.89 | 0.919 | 11.89 | 0.177 |
| J_H4 sensitive + recalibration | 0.68 | 0.89 | 0.911 | 10.11 | 0.156 |
| N114 recommended | 0.58 | 0.91 | 0.900 | 2.56 | 0.045 |
| N114 recommended + recalibration | 0.58 | 0.89 | 0.896 | 2.33 | 0.042 |
| N114 sensitive | 0.87 | 0.96 | 0.963 | 11.00 | 0.160 |
| N114 sensitive + recalibration | 0.87 | 0.96 | 0.963 | 11.00 | 0.160 |

No measurable gain: differences are about 0.02 in recall and two calls per library, in both directions. The
shifts are 0 at Q35+, where most calls are made, and the per-site control rates already absorb the rest. It
stays available (`recalibrate: true`) for reads whose qualities are less well calibrated than R10.4.1 SUP,
where it has not been tested.
