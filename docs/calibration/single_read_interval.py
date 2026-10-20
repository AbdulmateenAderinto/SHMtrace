"""Coverage of the single-read interval at low counts (audit F04, October 2026).

    python docs/calibration/single_read_interval.py

Repeated experiments with a known number of real single-read mutations in a sample
and Poisson (or overdispersed) background in the sample and in three controls at
the same depth. For each, the share of 95% intervals that contain the truth: the
profile-likelihood interval singletons.estimate now reports, and the normal
(Wald) interval it replaced. Efficiency is 1 here; it rescales both alike.
Writes single_read_interval/coverage.json; SINGLE_READ_INTERVAL.md reports it.
"""
import json
import math
from pathlib import Path

import numpy as np

from shmtrace.singletons import profile_interval

OUT = Path(__file__).parent / "single_read_interval"
REPS = 4000
CONTROLS = 3


def phi_of(kcs, ncs):
    kc, nc = sum(kcs), sum(ncs)
    if len(kcs) < 3 or not kc:
        return 1.0
    rate = kc / nc
    x2 = sum((k - n * rate) ** 2 / (n * rate) for k, n in zip(kcs, ncs, strict=True))
    return max(1.0, x2 / (len(kcs) - 1))


def wald(n, kc, tau, phi, z=1.96):
    exp = kc / tau
    var = n + (phi - 1.0) * exp + phi * exp * exp / max(kc, 1)
    se = math.sqrt(max(var, 0.0))
    return n - exp - z * se, n - exp + z * se


def run(signal, background, cv, seed):
    rng = np.random.default_rng(seed)
    hit_p = hit_w = zero_w = 0
    for _ in range(REPS):
        # each library's background rate varies around `background` with this CV
        if cv:
            shape = 1 / cv ** 2
            rates = rng.gamma(shape, background / shape, CONTROLS + 1)
        else:
            rates = np.full(CONTROLS + 1, background)
        n = int(rng.poisson(signal + rates[0]))
        kcs = [int(x) for x in rng.poisson(rates[1:])]
        phi = phi_of(kcs, [1.0] * CONTROLS)
        _, lo, hi = profile_interval(n, sum(kcs), float(CONTROLS), phi)
        hit_p += lo <= signal <= hi
        wl, wh = wald(n, sum(kcs), float(CONTROLS), phi)
        hit_w += wl <= signal <= wh
        zero_w += wh - wl == 0
    return {"signal": signal, "background": background, "cv": cv,
            "coverage_profile": hit_p / REPS, "coverage_wald": hit_w / REPS,
            "wald_zero_width": zero_w / REPS}


def main():
    rows = []
    seed = 1
    for cv in (0.0, 0.5):
        for background in (0.2, 1.0, 3.0, 10.0):
            for signal in (0, 1, 3, 10, 30):
                rows.append(run(signal, background, cv, seed))
                seed += 1
                r = rows[-1]
                print(f"cv {cv} bg {background:>4} signal {signal:>3}: profile {r['coverage_profile']:.3f}"
                      f"  wald {r['coverage_wald']:.3f}  (wald zero width {r['wald_zero_width']:.3f})")
    OUT.mkdir(exist_ok=True)
    (OUT / "coverage.json").write_text(json.dumps({"reps": REPS, "controls": CONTROLS, "rows": rows}, indent=1))


if __name__ == "__main__":
    main()
