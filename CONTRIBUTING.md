# Contributing

## Setup

```bash
git clone <repo-url> && cd nanoshm
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,web]"
pytest
```

## Ground rules

**The web layer never does science.** `src/nanoshm/webapp/` marshals uploads and
parameters into `nanoshm.pipeline.run()` and formats the result. Any analysis
logic belongs in the engine modules, where it can be tested without a server.

**Every parameter lives in `config.Params`.** If you add a knob, add it there,
give it a defensible default, and surface it in the presets rather than
hard-coding it at a call site.

**Determinism.** Given the same inputs, parameters and seed, output must not
change. Anything stochastic takes an explicit seed from `Params`.

**Changes to variant-calling rules need a test.** The filtering rules are the
scientific content of this tool. If you change one, add a test that pins the new
behaviour and say why in the PR.

## Before opening a PR

```bash
ruff check src tests
pytest
```

Update `CHANGELOG.md` under `[Unreleased]`.

## Calibration work

The open problem is the per-column error model (see `docs/ROADMAP.md` §6 and
Phase 1). If you are working on it, validate against the manuscript values on
the bundled example rather than only on unit tests.
