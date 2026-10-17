"""Numbers for the manuscript's conditional resampling validation, from saved runs."""
import json
import statistics
from collections import defaultdict

from scipy.stats import beta


def numbers(out):
    groups = defaultdict(list)
    for filename in ("fdr_validation", "fdr_null_reps", "fdr_estimator"):
        data = json.loads((out / f"{filename}.json").read_text())
        if data.get("revision_provenance", {}).get("complete") is False:
            raise RuntimeError(f"{filename} is still being regenerated; do not build stale validation claims")
        for row in data["runs"]:
            if row["controls"] != 3 or row["depth"] != 500:
                continue
            setting = "on" if row["exclude_noisy"] else "off"
            estimator = "max_" if filename == "fdr_estimator" else ""
            key = f"fdv_{estimator}{row['signal']}_{row['noise']}_{setting}"
            groups[key].append(row)
    result = {}
    for key, rows in groups.items():
        n = len(rows)
        k = sum(bool(r["any_discovery"]) for r in rows)
        result.update({
            f"{key}_n": str(n),
            f"{key}_any": str(k),
            f"{key}_excluded": str(sum(bool(r["decoys_set_aside"]) for r in rows)),
            f"{key}_true": f"{statistics.mean(r['true_calls'] for r in rows):.1f}",
            f"{key}_false": f"{statistics.mean(r['false_calls'] for r in rows):.1f}",
            f"{key}_fdp": f"{100 * statistics.mean(r['fdp'] for r in rows):.1f}",
            f"{key}_upper": f"{100 * beta.ppf(0.975, k + 1, n - k):.0f}" if k < n else "100",
        })
    return result
