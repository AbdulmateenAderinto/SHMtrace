"""AID targeting: how well called mutations follow a 5-mer targeting model.

The hotspot figures count mutations in binary motif classes (WRC/GYW and so
on). A 5-mer targeting model replaces the classes with a continuous mutability
for every 5-mer, learned from hundreds of thousands of silent mutations in
other repertoires (the S5F models of Yaari et al., Front Immunol 2013, as
implemented in SHazaM). It gives each position of the amplicon an expected
share of the mutations, and each substitution an expected direction.

Two uses here, both reported rather than used to make calls:

* **Observed against expected.** Real AID mutations should land where the
  model says AID lands; sequencing errors should not care. The targeting index
  is the mean model mutability at the mutated positions divided by the mean
  over the window: about 1 for errors, above 1 for AID. A permutation test
  over the window's positions gives its p-value.
* **A per-mutation annotation.** Every call carries the model's mutability of
  its 5-mer and the model's probability of its substitution, relative to the
  window. A likelihood-ratio term built on these could one day join the error
  model; it is not used to call anything yet, because that needs its own
  validation.

The model is a file, not part of nanoshm. SHazaM is AGPL-3 and nanoshm is MIT,
so nanoshm does not ship SHazaM's models; `scripts/export_shazam_model.R`
writes one from your own SHazaM installation in the format read here:

    fivemer,mutability,to_A,to_C,to_G,to_T

A model learned from one amplicon would be no use: each 5-mer occurs about
once in a 529-bp window, so a "model" fitted there is the positional map in
disguise. The expectation has to come from elsewhere.
"""
from __future__ import annotations

import csv
import random

BASES = "ACGT"
HOT = ("WRCY", "RGYW")       # classical AID hotspots, for marking a model plot
COLD = ("SYC", "GRS")        # classical coldspots
_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "W": "AT", "R": "AG",
          "Y": "CT", "S": "GC", "N": "ACGT"}


class TargetingModel:
    """Mutability for each 5-mer and the substitution probabilities of its centre."""

    def __init__(self, mutability: dict, substitution: dict, name: str = "model"):
        self.mutability = mutability          # "ACGTA" -> relative mutability
        self.substitution = substitution      # "ACGTA" -> {"A": p, ...}, centre excluded
        self.name = name

    def at(self, fivemer: str) -> float | None:
        return self.mutability.get(fivemer.upper())

    def sub(self, fivemer: str, alt: str) -> float | None:
        d = self.substitution.get(fivemer.upper())
        return None if d is None else d.get(alt)


def load(path: str, name: str | None = None) -> TargetingModel:
    """Read a model written by scripts/export_shazam_model.R (or any CSV of that shape)."""
    mut, sub = {}, {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            k = (row.get("fivemer") or "").strip().strip('"').upper()
            if len(k) != 5 or any(c not in BASES for c in k):
                continue                      # 5-mers containing N are skipped
            try:
                m = float(row["mutability"])
            except (TypeError, ValueError):
                continue
            if m != m:                        # NaN
                continue
            mut[k] = m
            probs = {}
            for b in BASES:
                v = row.get(f"to_{b}")
                try:
                    f = float(v)
                    if f == f:
                        probs[b] = f
                except (TypeError, ValueError):
                    pass
            sub[k] = probs
    if not mut:
        raise ValueError(f"{path}: no 5-mer mutabilities found; expected columns "
                         "fivemer,mutability,to_A,to_C,to_G,to_T")
    return TargetingModel(mut, sub, name or path.rsplit("/", 1)[-1].rsplit(".", 1)[0])


def window_profile(model: TargetingModel, ref: str, lo: int, hi: int) -> dict:
    """Model mutability at each analysed position (0-based lo..hi), and its mean.

    Positions within two bases of the reference ends have no full 5-mer and
    are left out, as are 5-mers the model does not cover.
    """
    ref = ref.upper()
    prof = {}
    for i in range(max(lo, 2), min(hi, len(ref) - 3) + 1):
        m = model.at(ref[i - 2:i + 3])
        if m is not None:
            prof[i] = m
    mean = sum(prof.values()) / len(prof) if prof else 0.0
    return {"mutability": prof, "mean": mean}


MIN_MUTATIONS = 10   # below this an index is noise, and is not reported


def targeting_index(positions0: list[int], profile: dict, n_perm: int = 2000,
                    seed: int = 0, min_n: int = MIN_MUTATIONS) -> dict:
    """Mean model mutability at mutated positions over the window mean.

    About 1 means the mutations ignore the model, as sequencing error does;
    above 1 means they land where AID lands. The one-sided p-value is the
    share of random position sets of the same size whose index is as high.
    """
    prof, mean = profile["mutability"], profile["mean"]
    hits = [prof[p] for p in positions0 if p in prof]
    if not hits or mean <= 0:
        return {"n": len(hits), "index": None, "p": None}
    if len(hits) < min_n:
        # one or four mutations landing on a hotspot says nothing about AID
        return {"n": len(hits), "index": None, "p": None, "note": f"fewer than {min_n} mutations"}
    obs = (sum(hits) / len(hits)) / mean
    pool = list(prof.values())
    k = len(hits)
    rng = random.Random(seed)
    ge = 0
    for _ in range(n_perm):
        draw = rng.sample(pool, k) if k <= len(pool) else rng.choices(pool, k=k)
        if (sum(draw) / k) / mean >= obs:
            ge += 1
    return {"n": k, "index": round(obs, 4), "p": (ge + 1) / (n_perm + 1)}


def expected_spectrum(model: TargetingModel, ref: str, profile: dict) -> dict:
    """The twelve-class spectrum the model predicts for this window.

    Each position contributes its mutability, split across the three possible
    substitutions by the model's substitution probabilities.
    """
    ref = ref.upper()
    out = {f"{r}>{a}": 0.0 for r in BASES for a in BASES if a != r}
    for i, m in profile["mutability"].items():
        r = ref[i]
        probs = model.substitution.get(ref[i - 2:i + 3]) or {}
        tot = sum(v for b, v in probs.items() if b != r) or 0.0
        for a in BASES:
            if a != r and tot:
                out[f"{r}>{a}"] += m * probs.get(a, 0.0) / tot
    s = sum(out.values()) or 1.0
    return {k: round(v / s, 6) for k, v in out.items()}


def annotate(variants, model: TargetingModel, ref: str, profile: dict) -> list[dict]:
    """Per call: the model's mutability and substitution probability, relative to the window."""
    ref = ref.upper()
    mean = profile["mean"] or 1.0
    out = []
    for v in variants:
        i = v.pos - 1
        k = ref[i - 2:i + 3] if 2 <= i < len(ref) - 2 else ""
        m = model.at(k) if k else None
        s = model.sub(k, v.alt) if k else None
        out.append({"pos": v.pos, "alt": v.alt, "fivemer": k,
                    "s5f_mutability_rel": round(m / mean, 4) if m is not None else None,
                    "s5f_substitution": round(s, 4) if s is not None else None})
    return out


# (motif, index of the mutated base): the motif sits so that base is the
# 5-mer's centre
_MOTIF_CENTRE = {"WRCY": 2, "RGYW": 1, "SYC": 2, "GRS": 0}


def matches(fivemer: str, motif: str) -> bool:
    """Whether a 5-mer's centre base is the mutated base of a hot or cold motif."""
    start = 2 - _MOTIF_CENTRE[motif]
    seg = fivemer[start:start + len(motif)]
    if start < 0 or len(seg) < len(motif):
        return False
    return all(c in _IUPAC[m] for c, m in zip(seg, motif, strict=True))
