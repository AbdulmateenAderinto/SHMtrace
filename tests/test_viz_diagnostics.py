"""The diagnostic figures: they draw when the data is there and keep quiet when it is not."""
import pytest

from nanoshm import project, viz_diagnostics as vd

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
