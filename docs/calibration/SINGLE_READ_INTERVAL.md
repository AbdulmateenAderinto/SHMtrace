# Single-read interval at low counts

*19 October 2026. Audit finding F04. Script: [`single_read_interval.py`](single_read_interval.py); numbers: `single_read_interval/coverage.json`.*

## What changed

`singletons.estimate` reported `estimate ± 1.96 × SE` with a plug-in variance. With no single-read events in the sample and none in the controls, that variance is zero, and the interval was 0 [0, 0]. A finite number of reads cannot show that the underlying rate is zero.

The interval now inverts a profile likelihood (Rolke, López and Conrad 2005) for the count model

- sample events ~ Poisson(signal + background)
- control events ~ Poisson(τ × background), with τ = control bases / sample bases

The likelihood ratio is divided by the between-control overdispersion φ that the estimator already computed (quasi-likelihood). The signal may be negative, as the background-subtracted point estimate can be. The point estimate is unchanged. With 0 sample events, 0 control events and τ = 3, the interval is now −0.64 to 1.92 events, divided by the efficiency.

## Coverage

Repeated experiments: one sample, three controls at the same depth, known signal (0, 1, 3, 10, 30 events), 4,000 replicates per cell. Efficiency is 1. In the overdispersed scenario, each library's background rate is drawn from a gamma distribution with a coefficient of variation of 0.5.

| background between libraries | background events | profile likelihood (new) | Wald (old) | Wald zero-width |
|---|---|---|---|---|
| Poisson | 0.2–10 | 0.925–0.988 | 0.717–0.998 | up to 45.4% |
| CV 0.5 | 0.2–1 | 0.928–0.988 | 0.703–0.998 | up to 45.8% |
| CV 0.5 | 3–10 | 0.867–0.965 | 0.867–0.943 | none |

## Limits

- **Overdispersed controls with several background events:** coverage falls to 0.87–0.90. Three controls estimate φ from two degrees of freedom, and the interval does not carry that uncertainty. The old interval had the same shortfall. More controls per run is the remedy. A heavier correction (an F reference with two degrees of freedom) would make every interval several times wider.
- **Efficiency is treated as known.** It is a ratio of large base counts.
- **The combined interval does not cover the called mutations.** `combined_density_ci_low/high` is the single-read interval shifted by the called count. The called mutations are held fixed, and every sample row says so in `combined_density_ci_scope`.
- **These results are from simulation only.** They validate the interval under the stated count model. They do not validate that model against the single-read events in real libraries.
