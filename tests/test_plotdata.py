"""Plot data: every figure's numbers as tidy tables, taken from the run without recalculation."""
import csv
import io
import zipfile

from shmtrace import plotdata


def _res():
    samples = [
        {"sample": "WT1", "group": "Wt", "control": False, "reads": 500, "unique_mutations": 40,
         "unique_density_x1e5": 15.2, "one_event_density_x1e5": 0.38, "density_bases": 263000,
         "mutation_frequency_x1e3": 2.1, "mutations_total": 550, "frequency_bases": 262000,
         "mutated_reads_pct": 20.0, "pct_GC_in_WRC_GYW": 61.0, "gc_mutations": 18},
        {"sample": "WT2", "group": "Wt", "control": False, "reads": 500, "unique_mutations": 55,
         "unique_density_x1e5": 20.9, "mutated_reads_pct": 25.0, "pct_GC_in_WRC_GYW": 58.0, "gc_mutations": 22},
        {"sample": "KO1", "group": "KO", "control": True, "reads": 500, "unique_mutations": 0,
         "unique_density_x1e5": 0, "one_event_density_x1e5": 0.38, "mutated_reads_pct": 0.4},
    ]
    return {
        "samples": samples, "group_order": ["Wt", "KO"],
        "groups": [{"group": "Wt", "n": 2, "control": False, "compared_with": "KO",
                    "unique_density_x1e5_mean": 18.05, "unique_density_x1e5_p_welch_log_vs_KO": 0.004,
                    "mutated_reads_pct_p_mwu_vs_KO": 0.2}],
        "spectrum_group": {"Wt": {"C>T": 30, "G>A": 10}, "KO": {}},
        "spectrum_sample": {"WT1": {"C>T": 15, "G>A": 5}, "WT2": {"C>T": 15, "G>A": 5}, "KO1": {}},
        "mutation_load": {"WT1": {"0": 400, "1": 60, "5": 40}, "KO1": {"0": 498, "1": 2}},
        "positional": [{"group": "Wt", "pos": 12, "ref": "C", "mutated_reads": 4, "frequency_pct": 0.8}],
        "positional_sample": [{"sample": "WT1", "group": "Wt", "pos": 12, "mutated_reads": 4, "frequency_pct": 0.8}],
        "hotspots": [12], "reference": "A" * 11 + "C" + "A" * 20,
    }


def test_each_figure_has_a_tidy_table_with_the_runs_own_numbers():
    t = plotdata.tables(_res())
    d = {r["sample"]: r for r in t["density_by_group.csv"]}
    assert d["WT1"]["unique_mutations_per_1e5_bp"] == 15.2 and d["KO1"]["zero_drawn_at"] == 0.38
    assert d["KO1"]["control"] == "yes"
    m = [r for r in t["substitution_matrix.csv"] if r["group"] == "Wt" and r["from_base"] == "C" and r["to_base"] == "T"]
    assert m[0]["count"] == 30 and m[0]["percent_of_group"] == 75.0
    assert len([r for r in t["substitution_matrix.csv"] if r["group"] == "Wt"]) == 12
    sp = [r for r in t["spectrum_by_group.csv"] if r["group"] == "Wt" and r["substitution"] == "C>T"][0]
    assert sp["mean_percent"] == 75.0 and sp["samples"] == 2 and sp["group_calls"] == 40
    load = {(r["group"], r["load_class"]): r["reads"] for r in t["mutation_load.csv"]}
    assert load[("Wt", "0")] == 400 and load[("Wt", "4-10")] == 40 and load[("KO", "1-3")] == 2
    pos = t["positional_by_sample.csv"][0]
    assert pos["reference_base"] == "C" and pos["wrc_gyw_hotspot"] == "yes"


def test_long_tables_stack_measures_and_statistics():
    t = plotdata.tables(_res())
    long = t["per_sample_long.csv"]
    assert {"sample", "group", "control", "measure", "value", "unit"} == set(long[0])
    assert {r["measure"] for r in long if r["sample"] == "WT1"} >= {"unique_mutations_per_1e5_bp", "reads_mutated_pct"}
    stats = {(r["measure"], r["statistic"]): r["value"] for r in t["group_statistics_long.csv"]}
    assert stats[("unique_mutations_per_1e5_bp", "welch_log")] == 0.004
    assert stats[("reads_mutated_pct", "mann_whitney")] == 0.2
    assert stats[("unique_mutations_per_1e5_bp", "mean")] == 18.05


def test_prism_files_have_one_column_per_group_and_the_zip_has_a_dictionary():
    files = plotdata.files(_res())
    rows = list(csv.reader(io.StringIO(files["prism/unique_mutations_per_1e5_bp.csv"])))
    assert rows[0] == ["Wt", "KO"] and rows[1] == ["15.2", "0"] and rows[2] == ["20.9", ""]
    dic = list(csv.DictReader(io.StringIO(files["data_dictionary.csv"])))
    cols = {(r["file"], r["column"]) for r in dic}
    assert ("density_by_group.csv", "zero_drawn_at") in cols
    assert all(r["meaning"] for r in dic if r["file"] == "density_by_group.csv")
    z = zipfile.ZipFile(io.BytesIO(plotdata.zip_bytes(_res())))
    assert "README.txt" in z.namelist() and "per_sample_long.csv" in z.namelist()
