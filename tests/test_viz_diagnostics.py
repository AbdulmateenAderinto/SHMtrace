"""The diagnostic figures: they draw when the data is there and keep quiet when it is not."""
import pytest

from shmtrace import project, viz_diagnostics as vd

PNG = b"\x89PNG"


@pytest.fixture(scope="module")
def jh4(examples):
    """A small two-sample project with indels and the saturation curve on."""
    return project.run_project({
        "name": "t", "region": "jh4_intron",
        "reference": str(examples / "JH4_reference.fasta"),
        "contaminants": [str(examples / "JH4_contaminant.fasta")],
        "saturation": [50, 120], "saturation_draws": 2,
        "overrides": {"exclude_indels": False},
        "samples": [
            {"name": "WT", "path": str(examples / "WT3_JH4_sub500.fastq"), "group": "WT"},
            {"name": "KO", "path": str(examples / "KO3_JH4_sub500.fastq"),
             "group": "KO", "control": True},
        ]})


def test_every_drawable_figure_renders(jh4):
    names = vd.drawable(jh4)
    assert "quality_histogram" in names and "volcano" in names
    figs = vd.diagnostic_figures(jh4, "png")
    for name in names:
        assert figs[f"{name}.png"].startswith(PNG), name


def test_drawable_agrees_with_what_is_actually_drawn(jh4):
    """The cheap predicate must not promise a figure the renderer then skips."""
    drawn = {k.rsplit(".", 1)[0] for k in vd.diagnostic_figures(jh4, "png")}
    assert set(vd.drawable(jh4)) == drawn


def test_saturation_curve_is_measured_when_asked(jh4):
    rows = jh4["_details"]["WT"]["rarefaction"]
    assert [r["depth"] for r in rows] == [50, 120]
    # more reads cannot find fewer unique mutations
    assert rows[0]["unique_variants_mean"] <= rows[1]["unique_variants_mean"]


def test_nothing_is_drawn_from_an_empty_result():
    empty = {"samples": [], "group_order": [], "groups": []}
    assert vd.drawable(empty) == []
    assert vd.diagnostic_figures(empty, "png") == {}


def test_a_figure_that_raises_is_skipped_not_fatal(monkeypatch, jh4):
    def boom(res, fmt="png"):
        raise RuntimeError("no")
    monkeypatch.setattr(vd, "FIGURES",
                        [("boom", boom, "boom", lambda r: True)] + list(vd.FIGURES)[:1])
    out = vd.diagnostic_figures(jh4, "png")
    assert "boom.png" not in out and len(out) == 1


def test_public_lists_the_figures_without_drawing_them(jh4, monkeypatch):
    for name, fn, _label, _ok in vd.FIGURES:
        monkeypatch.setattr(vd, fn.__name__,
                            lambda *a, **k: pytest.fail("public() must not render"))
    pub = project.public(jh4)
    assert "quality_histogram" in pub["available_figures"]
    assert "density_by_group" in pub["available_figures"]
    assert not any(k.startswith("_") for k in pub)


def test_svg_and_pdf_come_out_too(jh4):
    assert vd.quality_histogram(jh4, "svg").lstrip().startswith(b"<?xml")
    assert vd.volcano(jh4, "pdf").startswith(b"%PDF")


def test_corrected_spectrum_uses_shmtool_denominator():
    from collections import Counter
    from shmtrace.project import corrected_spectrum
    # 6 C>T over 3 animals in a window with 20 C: 6 / (20 x 3)
    out = corrected_spectrum(Counter({"C>T": 6}), {"A": 10, "C": 20, "G": 5, "T": 0}, 3)
    assert out["C>T"] == 0.1
    assert out["A>G"] == 0.0
    assert out["T>C"] is None          # no T in the window: undefined, not zero
    assert len(out) == 12


def test_coverage_dips_name_the_low_stretch_only():
    from shmtrace.project import coverage_dips
    flat = [500] * 200
    assert coverage_dips(flat, 101) == ""
    dip = [500] * 80 + [100] * 30 + [500] * 90
    msg = coverage_dips(dip, 101)
    assert "181-210" in msg and "500" in msg
    # a dip shorter than the minimum run is a single bad column, not a stretch
    assert coverage_dips([500] * 100 + [10] * 3 + [500] * 100, 1) == ""


def test_project_reports_coverage_and_corrected_spectrum(jh4):
    assert set(jh4["coverage"]) == {"WT", "KO"}
    c = jh4["coverage"]["WT"]
    assert c["start"] == 101 and len(c["depth"]) == 529
    assert sum(jh4["base_composition"].values()) == 529
    assert "C>T" in jh4["spectrum_corrected_group"]["WT"]
