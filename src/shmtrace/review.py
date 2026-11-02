"""A suggested verdict for each called mutation, from the evidence the run already holds.

Every call SHMtrace lists has passed the calibrated test, so a suggestion here is not a
second test. It weighs what the test does not see directly, and says which way each item
points, so a person checking calls by eye starts from the same facts every time:

Decisive (suggest rejecting whatever else holds):
- the call is in a control library, where by design nothing is somatic;
- the controls' own error rate at this position predicts as many supporting reads at this
  depth as the call has;
- the same change is called in a control library (a germline difference, a PCR or
  sequencing artefact, or carry-over).

Against: supporting bases of low quality; most supporting reads ending within a few bases
of the call (where alignments are least reliable); a long homopolymer run; every
supporting read on one strand when the library has both; a one-read call in a sample
whose one-read events are mostly background by SHMtrace's single-read estimate.

For: three or more supporting reads; every supporting base at high quality; support well
above what the controls' error rate predicts; a P value far past the threshold.

Two things are deliberately not used: whether the call sits in an AID hotspot, and whether
other samples share it. Both are what an SHM analysis goes on to measure, so letting them
steer which calls survive would make the hotspot and recurrence results circular.

Suggest accepting at a score of +2 or more with nothing strongly against; rejecting on a
decisive item, a score of -2 or less, or something strongly against (-2) that the rest does
not outweigh; otherwise checking by eye. The weights are plain
integers so the reasons can be read off directly.
"""
from __future__ import annotations

from statistics import median

LOW_Q = 20            # Q20: a 1-in-100 chance the base is wrong, near the somatic rate itself
GOOD_Q = 30           # Q30: 1 in 1,000
END_BP = 10           # reads ending this close to a call support it less
FAR_ABOVE_BACKGROUND = 5
STRONG_P = 1e-6
BACKGROUND_SHARE = 0.5


def _key(v: dict) -> str:
    return f"{v['sample']}|{v['pos']}|{v['ref']}>{v['alt']}"


def _support(res: dict) -> dict:
    """(sample, pos, alt) -> [(quality, near an end, strand)] over the analysed reads."""
    out: dict = {}
    reads = {**(res.get("_call_reads") or {}), **(res.get("_alignments") or {})}
    for sample, data in reads.items():
        for r in data.get("reads", []):
            if r.get("fate") != "kept":
                continue
            for pos, base, q in r.get("diffs", []):
                near = pos - r["start"] < END_BP or r["end"] - pos < END_BP
                out.setdefault((sample, pos, base), []).append((q, near, r.get("strand", 1)))
    return out


def suggest(res: dict) -> dict[str, dict]:
    """{call key: {"suggestion": accept|reject|check, "score": int, "reasons": [...]}}.

    Each reason is {"weight": -2..2, "decisive": bool, "text": ...}.
    """
    samples = {s["sample"]: s for s in res.get("samples") or []}
    controls = {n for n, s in samples.items() if s.get("control")}
    # the controls' rate for this very change; the per-position figure is the highest of the
    # three alternatives, which would overstate the background for the other two
    bg_alt = {k: float(v) for k, v in (res.get("control_background_alt") or {}).items()}
    hp_flag = int((res.get("params") or {}).get("homopolymer_flag") or 4)
    variants = res.get("variants") or []
    in_control: dict = {}
    for v in variants:
        if v["sample"] in controls:
            in_control.setdefault((v["pos"], v["alt"]), []).append(v["sample"])
    both_strands = {n: any(v["minus_reads"] for v in variants if v["sample"] == n)
                    and any(v["plus_reads"] for v in variants if v["sample"] == n) for n in samples}
    support = _support(res)

    out = {}
    for v in variants:
        reasons = []

        def say(weight, text, decisive=False):
            reasons.append({"weight": weight, "decisive": decisive, "text": text})

        name, s = v["sample"], samples.get(v["sample"], {})
        if name in controls:
            say(-2, "Called in a control library, where nothing is expected to be somatic.", True)
        expected = bg_alt.get(f"{v['pos']}{v['alt']}", 0.0) * v["depth"]
        if expected and v["reads"] <= expected:
            say(-2, f"The controls' error rate here predicts {expected:.1f} supporting reads at this depth; "
                    f"the call has {v['reads']}.", True)
        elif expected and v["reads"] >= FAR_ABOVE_BACKGROUND * expected:
            say(1, f"{v['reads'] / expected:.0f} times the reads the controls' error rate here predicts "
                   f"({expected:.1f}).")
        shared = [c for c in in_control.get((v["pos"], v["alt"]), []) if c != name]
        if shared and name not in controls:
            say(-2, f"The same change is called in {'a control' if len(shared) == 1 else 'controls'} "
                    f"({', '.join(shared)}).", True)

        qs = [q for q, _n, _s in support.get((name, v["pos"], v["alt"]), [])]
        if qs:
            m = median(qs)
            if m < LOW_Q:
                say(-2, f"The supporting bases are of low quality (median {m:.0f}).")
            elif sum(q < LOW_Q for q in qs) > len(qs) / 2:
                say(-1, f"Most supporting bases are below quality {LOW_Q}.")
            elif min(qs) >= GOOD_Q:
                say(1, f"Every supporting base is quality {GOOD_Q} or better.")
            ends = sum(n for _q, n, _s in support[(name, v["pos"], v["alt"])])
            if ends and ends >= len(qs) / 2:
                say(-1, f"{ends} of {len(qs)} supporting reads end within {END_BP} bp of the call.")
        if (v.get("homopolymer") or 1) >= hp_flag:
            say(-1, f"Next to a run of {v['homopolymer']} identical bases, where reads miscall most often.")
        if both_strands.get(name) and v["reads"] >= 3 and (not v["plus_reads"] or not v["minus_reads"]):
            say(-1, f"Every supporting read is on one strand ({v['plus_reads']} forward, "
                    f"{v['minus_reads']} reverse) although the library has both.")
        if v["reads"] >= 3:
            say(1, f"Seen in {v['reads']} reads.")
        elif v["reads"] == 1:
            obs, back = s.get("single_read_observed"), s.get("single_read_background")
            if obs and back is not None:
                share = min(1.0, back / obs)
                say(-1 if share >= BACKGROUND_SHARE else 0,
                    f"Seen in one read. In this sample about {100 * share:.0f}% of one-read calls are "
                    f"expected to be background (single-read estimate).")
        p = v.get("p_value")
        if p is not None and p < STRONG_P:
            say(1, f"P = {p:.1g}, far past the calibrated threshold.")

        score = sum(r["weight"] for r in reasons)
        strong_against = any(r["weight"] <= -2 for r in reasons)
        if any(r["decisive"] for r in reasons) or score <= -2 or (strong_against and score <= 0):
            verdict = "reject"
        elif score >= 2 and not strong_against:
            verdict = "accept"
        else:
            verdict = "check"
        out[_key(v)] = {"suggestion": verdict, "score": score, "reasons": reasons}
    return out


def cached(res: dict) -> dict[str, dict]:
    """suggest(res), worked out once per run."""
    if "_suggestions" not in res:
        res["_suggestions"] = suggest(res)
    return res["_suggestions"]
