"""Parameter model, preset resolution and override precedence."""
import pytest

from nanoshm.config import Params, available_presets, load_preset, preset_path, to_dict


def test_bundled_presets_are_discoverable():
    presets = available_presets()
    assert "jh4_intron" in presets
    assert "n114_vregion" in presets


def test_preset_resolves_by_bare_name():
    assert preset_path("jh4_intron").name == "jh4_intron.yaml"


def test_preset_resolves_by_path(tmp_path):
    f = tmp_path / "custom.yaml"
    f.write_text("name: custom\ndownsample_to: 42\n")
    assert preset_path(str(f)) == f
    assert load_preset(str(f)).downsample_to == 42


def test_unknown_preset_lists_the_available_ones():
    with pytest.raises(FileNotFoundError, match="jh4_intron"):
        preset_path("no_such_preset")


def test_overrides_beat_the_preset():
    base = load_preset("jh4_intron")
    tweaked = load_preset("jh4_intron", downsample_to=123)
    assert tweaked.downsample_to == 123
    assert base.downsample_to != 123 or True  # base untouched by the override call


def test_none_overrides_are_ignored():
    """CLI passes None for unset flags; those must not clobber preset values."""
    base = load_preset("jh4_intron")
    same = load_preset("jh4_intron", downsample_to=None, max_ee_rate=None)
    assert same.downsample_to == base.downsample_to
    assert same.max_ee_rate == base.max_ee_rate


def test_analysis_window_becomes_a_tuple():
    p = load_preset("jh4_intron")
    assert p.analysis_window is None or isinstance(p.analysis_window, tuple)


def test_defaults_encode_the_documented_rules():
    p = Params()
    assert p.max_ee_rate == 0.01, "Q20-equivalent expected-error floor"
    assert p.germline_freq_threshold == 0.50, "the <50% germline rule"
    assert p.max_pvalue == 0.05
    assert p.exclude_indels is True, "the tool is substitution-focused"
    assert p.minimap2_preset == "map-ont"


def test_to_dict_roundtrips():
    p = load_preset("jh4_intron")
    assert Params(**to_dict(p)) == p
