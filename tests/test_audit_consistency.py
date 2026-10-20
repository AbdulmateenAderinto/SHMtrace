"""Scientific invariants from the October 2026 computational methods audit."""
import copy
import importlib.util
import itertools
import json
import pickle
import random
import tracemalloc
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from shmtrace import batch, clones, indels, null_model, pipeline, project, umi, variants
from shmtrace.align import AlignedRead
from shmtrace.config import Params
from shmtrace.io_utils import FASTA_DEFAULT_Q, read_fasta
from shmtrace.parallel import pmap
from shmtrace.recal import Recal


def library():
    reads = []
    for i in range(10):
        pairs = []
        for pos in range(14):
            if pos == 9 and i >= 5:  # unequal coverage is intentional
                continue
            base, q = 'A', 40
            if pos in (1, 2, 7, 8, 9) and i == 0:
                base, q = 'C', {1: 10, 2: 30, 7: 29, 8: 31, 9: 40}[pos]
            elif pos == 3 or (pos == 4 and i < 6):
                base = 'C'  # fixation and a dominant allele
            elif pos == 5:
                base = 'C' if i < 4 else 'G' if i < 8 else 'A'
            elif pos == 6 and i < 2:
                base, q = 'C', 10  # quality floor is singleton-only
            pairs.append((pos, base, q))
        reads.append(AlignedRead(str(i), 'target', 0, 14, pairs, 'A' * 14,
                                 1 if i % 2 else -1))
    return reads


@pytest.mark.parametrize('floor', [None, 30, 31])
@pytest.mark.parametrize('ancestral', [False, True])
@pytest.mark.parametrize('recalibrate', [False, True])
@pytest.mark.parametrize('window,minq,trim,min_reads', [
    (None, 0, 0, 1), ((2, 10), 20, 2, 1), (None, 0, 1, 2)])
def test_sample_and_decoy_candidates_have_identical_pvalues(
        tmp_path, floor, ancestral, recalibrate, window, minq, trim, min_reads):
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>target\n' + 'A' * 14 + '\n')
    p = Params(error_model='qcal', multiple_testing='none', max_pvalue=1,
               strand_filter='off', single_read_call_min_quality=floor,
               ancestral_reference=ancestral, analysis_window=window,
               min_base_quality=minq, end_trim_bp=trim, min_variant_reads=min_reads)
    null = null_model.NullModel(qtail={('A', a, s): [1e-4] * 61
                                       for a in 'CGT' for s in (-1, 1)})
    if recalibrate:
        null.recal = Recal([[max(0, q - 5)] * 3 for q in range(61)])
    reads = library()
    grid = variants.qcal_pvalues_by_qscale(reads, str(fasta), 'target', p, null, [1, 10])
    for scale in (1, 10):
        null.qscale = scale
        _, stats = variants.call_variants(reads, str(fasta), 'target', p, null)
        candidates = stats['_candidate_keys']
        assert len(candidates) == len(grid[scale])
        assert sorted(pv for _, _, pv in candidates) == pytest.approx(grid[scale])
        keys = {(pos, alt) for pos, alt, _ in candidates}
        if window is None and minq == 0 and min_reads == 1:
            # Both alternatives remain below 50%, even though their sum is 80%.
            assert {(6, 'C'), (6, 'G'), (7, 'C')} <= keys
            assert ((4, 'C') in keys) == ancestral
            assert ((5, 'C') in keys) == ancestral
            assert ((2, 'C') in keys) == (floor is None)


def test_panel_file_roundtrip_preserves_threshold_calls_and_recalibration(tmp_path):
    ref = tmp_path / 'ref.fa'
    ref.write_text('>target\nA\n')
    reads = [AlignedRead(str(i), 'target', 0, 1, [(0, b, 40)], b)
             for i, b in enumerate('CAAAAAAAAA')]
    null = null_model.NullModel(qtail={('A', 'C', 1): [1e-6] * 61}, qscale=5,
                               meta={'run': 'same-run'})
    null.decoys = {'K1': ([1e-12], 10), 'K2': ([], 20)}
    null.recal = Recal([[q] * 3 for q in range(61)])
    path = tmp_path / 'panel.json'
    null_model.save_panel(null, str(path))
    back = null_model.load_panel(str(path))
    assert json.loads(path.read_text())['version'] == 3
    assert back.decoys == null.decoys
    assert back.meta == null.meta and back.qscale == null.qscale
    assert back.recal.to_json() == null.recal.to_json()
    p = Params(error_model='qcal', multiple_testing='control', strand_filter='off')
    before = variants.call_variants(reads, str(ref), 'target', p, null)
    after = variants.call_variants(reads, str(ref), 'target', p, back)
    assert after == before
    assert after[1]['multiple_testing'] == 'control'


@pytest.mark.parametrize('version', [1, 2])
def test_legacy_rate_only_panels_remain_readable_but_warn_about_bh(tmp_path, version):
    ref = tmp_path / 'ref.fa'
    ref.write_text('>target\nA\n')
    null = null_model.NullModel(qtail={('A', 'C', 1): [1e-6] * 61})
    data = null.to_json()
    data['version'] = version
    del data['decoys']
    back = null_model.NullModel.from_json(data)
    assert back.decoys == {} and back.qtail == null.qtail
    reads = [AlignedRead(str(i), 'target', 0, 1, [(0, b, 40)], b)
             for i, b in enumerate('CAAAAAAAAA')]
    p = Params(error_model='qcal', multiple_testing='control', strand_filter='off')
    result = pipeline.analyse(reads, str(ref), p, null_rates=back)
    assert result['summary']['calling']['multiple_testing'] == 'bh'
    assert any('Benjamini-Hochberg' in x for x in result['summary']['warnings'])
    metrics = project.sample_metrics('case', result, 'A', 0, 0, 'case', False)
    assert metrics['warnings'].count('control-based thresholding requested') == 1


@pytest.mark.parametrize('empty_name', ['case', 'control'])
@pytest.mark.parametrize('depth', [None, 5])
def test_empty_case_or_control_stops_batch_before_all_samples_become_zero(
        tmp_path, monkeypatch, empty_name, depth):
    ref = tmp_path / 'ref.fa'
    ref.write_text('>target\nA\n')
    good = [AlignedRead(str(i), 'target', 0, 1, [(0, 'A', 40)], 'A') for i in range(10)]
    def prepare(fq, *_):
        return ([] if fq == empty_name else good), {}, {}, []
    monkeypatch.setattr(batch, '_prepare_one', prepare)
    with pytest.raises(ValueError, match=f'No usable reads.*{empty_name}'):
        batch.prepare_batch({'case': 'case', 'control': 'control'}, str(ref), [],
                            Params(phase_reads=False), controls={'control'}, depth=depth)
    assert len(good) == 10


@pytest.mark.parametrize('depth,configured', [(0, None), (-1, None), (None, 0), (None, -1)])
def test_nonpositive_batch_depth_is_rejected(tmp_path, monkeypatch, depth, configured):
    ref = tmp_path / 'ref.fa'
    ref.write_text('>target\nA\n')
    read = AlignedRead('r', 'target', 0, 1, [(0, 'A', 40)], 'A')
    monkeypatch.setattr(batch, '_prepare_one', lambda *_: ([read], {}, {}, []))
    with pytest.raises(ValueError, match='positive'):
        batch.prepare_batch({'case': 'case'}, str(ref), [],
                            Params(phase_reads=False, downsample_to=configured), depth=depth)


@pytest.mark.parametrize('umi_mode', ['header', 'flank'])
def test_umi_indels_fail_explicitly_in_consensus_and_direct_caller(umi_mode):
    p = Params(umi=umi_mode, exclude_indels=False)
    with pytest.raises(ValueError, match='Indel analysis after UMI consensus is not supported'):
        umi.consensus_reads([], 'A', p)
    with pytest.raises(ValueError, match='Indel analysis after UMI consensus is not supported'):
        indels.call(SimpleNamespace(ref='A'), None, p)


def test_spikein_harness_refits_controls_for_a_different_candidate_policy(tmp_path, monkeypatch):
    crm = Path(__file__).resolve().parents[1] / 'docs/calibration/crm'
    monkeypatch.syspath_prepend(str(crm))
    spec = importlib.util.spec_from_file_location('spikein_audit_test', crm / 'spikein.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    reference = tmp_path / 'ref.fa'
    reference.write_text('>target\n' + 'A' * 100 + '\n')
    monkeypatch.setattr(module, 'DATASETS', {'toy': {
        'preset': 'generic_amplicon', 'reference': str(reference), 'depth': 1,
        'samples': {name: {'fastq': name, 'control': True} for name in ('K1', 'K2')}}})
    monkeypatch.setattr(module, 'OUT', tmp_path)
    monkeypatch.setattr(module, 'PROFILES_RUN', ['recommended', 'recommended_q30'])
    for flag in ('WITH_LOFREQ', 'WITH_CLAIRS', 'WITH_CLAIR3', 'WITH_DV', 'WITH_DS'):
        monkeypatch.setattr(module, flag, False, raising=False)
    read = SimpleNamespace(qname='r')
    monkeypatch.setattr(module.pipeline, 'prepare_reads', lambda *_: ([read], {}))
    monkeypatch.setattr(module, 'read_fastq', lambda *_: [('r', 'A' * 100, 'I' * 100)])
    monkeypatch.setattr(module.simulate, 'choose_sites', lambda *_: [])
    monkeypatch.setattr(module.simulate, 'spike', lambda reads, *_: (reads, {}))
    fitted = []
    def fit(peers, reference, contaminants, p):
        fitted.append((tuple(peers), p.single_read_call_min_quality))
        return SimpleNamespace(floor=p.single_read_call_min_quality), {}
    monkeypatch.setattr(module.null_model, 'fit', fit)
    def run(fastq, reference, contaminants, p, null):
        assert null.floor == p.single_read_call_min_quality
        return {'variants': [], '_aligned': [read]}
    monkeypatch.setattr(module.pipeline, 'run', run)
    module.run('toy', reps=1, n_sites=0)
    assert fitted == [(('K2',), None), (('K2',), 30),
                      (('K1',), None), (('K1',), 30)]
    assert len(json.loads((tmp_path / 'spikein_toy.json').read_text())['rows']) == 4


# ---- round 2: zero counts (F05), sparse intervals (F04), legacy compare (F06),
# ---- exposures (F08) ----------------------------------------------------------

def _row(name, group, k, control, n=1000, reads=10, bp=100):
    return {"sample": name, "group": group, "control": control, "reads": reads,
            "analysed_bp": bp, "unique_mutations": k, "mutations_total": k,
            "density_bases": n, "frequency_bases": n,
            "unique_density_x1e5": k * 1e5 / n, "one_event_density_x1e5": 1e5 / n,
            "mutation_frequency_x1e3": k * 1e3 / n}


def _rows(ctrl=(0, 0, 1), case=(2, 3, 4), n=1000):
    return ([_row(f'K{i}', 'KO', k, True, n) for i, k in enumerate(ctrl)]
            + [_row(f'W{i}', 'WT', k, False, n) for i, k in enumerate(case)])


@pytest.mark.filterwarnings('ignore::RuntimeWarning')
def test_one_zero_count_rule_gives_one_p_across_project_rows_and_animal_tests():
    from shmtrace import replicates
    rows = _rows()
    g = next(r for r in project.group_table(rows, ['unique_density_x1e5', 'mutation_frequency_x1e3'])
             if r['group'] == 'WT')
    tests = {t['metric']: t for t in replicates.group_tests(rows, boot=19)}
    primary = [g['unique_density_x1e5_p_welch_log_vs_KO'], g['mutation_frequency_x1e3_p_welch_log_vs_KO'],
               tests['unique_density']['p_welch_log'], tests['mutation_frequency']['p_welch_log']]
    # the audit's example: 0.034 / 0.022 / 0.008 before; one primary value now
    assert primary == pytest.approx([0.03418] * 4, abs=5e-5)
    # the other conventions survive only as labelled sensitivity analyses
    assert g['unique_density_x1e5_sensitivity_p_welch_log_half_event_vs_KO'] == pytest.approx(0.00805, abs=5e-5)
    assert g['unique_density_x1e5_sensitivity_p_welch_log_third_of_smallest_vs_KO'] == pytest.approx(0.02173, abs=5e-5)
    assert tests['unique_density']['p_welch_log_half_event'] == pytest.approx(0.00805, abs=5e-5)
    assert tests['unique_density']['p_welch_log_third_of_smallest'] == pytest.approx(0.02173, abs=5e-5)
    assert 'fold_is_lower_bound' not in tests['unique_density']


@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.parametrize('scale', [1.0, 1e3, 1e5])
def test_log_welch_is_invariant_to_the_unit_of_the_rate(scale):
    from shmtrace import replicates
    ks_c, ks_x, n = [0, 0, 1], [2, 3, 4], [900, 1000, 1100]
    def logs(ks):
        return replicates.censored_log10([k * scale / m for k, m in zip(ks, n, strict=True)], [scale / m for m in n])
    p = replicates.welch_p(logs(ks_x), logs(ks_c))
    assert p == pytest.approx(replicates.welch_p(
        replicates._log_rates(ks_x, n), replicates._log_rates(ks_c, n)))


def test_mixed_and_all_zero_groups_do_not_crash():
    rows = _rows(ctrl=(0, 0, 0), case=(0, 5, 2))
    g = next(r for r in project.group_table(rows, ['unique_density_x1e5']) if r['group'] == 'WT')
    assert g['unique_density_x1e5_p_welch_log_vs_KO'] is not None
    rows = _rows(ctrl=(0, 0, 0), case=(0, 0, 0))
    g = next(r for r in project.group_table(rows, ['unique_density_x1e5']) if r['group'] == 'WT')
    assert g['unique_density_x1e5_p_welch_log_vs_KO'] is None   # no variance: nothing to test


def test_zero_single_read_events_do_not_give_a_zero_width_interval():
    from shmtrace.singletons import profile_interval
    est, lo, hi = profile_interval(0, 0, 3.0)
    assert est == 0 and lo < 0 < hi
    assert hi == pytest.approx(1.92, abs=0.01)         # chi2(1) / 2 events
    _est, _lo, hi_od = profile_interval(0, 0, 3.0, phi=3.0)
    assert hi_od > hi                                  # overdispersion widens it
    for n, m in [(0, 6), (5, 0), (10, 6), (200, 30)]:
        est, lo, hi = profile_interval(n, m, 3.0)
        assert lo < est < hi


def test_profile_interval_covers_at_low_counts():
    import numpy as np

    from shmtrace.singletons import profile_interval
    rng = np.random.default_rng(7)
    for signal, bg in [(0, 0.2), (1, 1.0), (3, 3.0)]:
        hits = 0
        for _ in range(600):
            n, m = rng.poisson(signal + bg), rng.poisson(3 * bg)
            _, lo, hi = profile_interval(int(n), int(m), 3.0)
            hits += lo <= signal <= hi
        assert hits / 600 >= 0.91, (signal, bg, hits / 600)


def test_single_read_estimate_reports_its_interval_method():
    from shmtrace import singletons
    reads = [AlignedRead(str(i), 'target', 0, 20, [(j, 'A', 40) for j in range(20)], 'A' * 20)
             for i in range(10)]
    s = singletons.tally(reads, 'A' * 20, 0, 19)
    c = [singletons.tally(reads, 'A' * 20, 0, 19) for _ in range(3)]
    est = singletons.estimate(s, c, 30, 0.5, 'A' * 20)
    assert est['observed'] == 0 and est['control_events'] == 0
    assert est['ci_low'] < 0 < est['ci_high']
    assert 'profile likelihood' in est['ci_method']


def _fake_result(model, notes=()):
    return {"variants": [], "summary": {"reads": {"aligned_to_target": 5},
                                        "calling": {"error_model": model, "multiple_testing": "bh",
                                                    "notes": list(notes)}}}


def test_legacy_compare_says_it_is_pooled_and_names_the_model_each_file_used(tmp_path, monkeypatch):
    from shmtrace import compare
    ref = tmp_path / 'ref.fa'
    ref.write_text('>t\n' + 'ACGT' * 25 + '\n')
    fallback = ("error_model 'qcal' needs control reads (--null or a panel); "
                "used the read-quality model instead")
    def groups(control_files, case_files, *_):
        out = {"control_0": _fake_result('read_quality', [fallback]),
               "case_0": _fake_result('read_quality', [fallback])}
        return ([out["control_0"]], [out["case_0"]],
                {"depth": 5, "known_blocks": [], "calling": compare._calling_provenance(out)})
    monkeypatch.setattr(compare, '_run_groups', groups)
    res = compare.compare(['c'], ['x'], str(ref), [], Params())
    assert res['analysis_route'].startswith('legacy pooled')
    assert res['calling']['error_models_used'] == {'read_quality': 2}
    assert any('not evidence at the level of biological replicates' in w for w in res['warnings'])
    # a missing empirical null cannot pass silently as the calibrated analysis
    assert any('read-quality model' in w for w in res['warnings'])
    md = compare.report_markdown(res)
    assert md.splitlines()[0].endswith('(legacy, pooled, descriptive)')
    assert 'pooled p' in md and 'read_quality × 2' in md


def test_density_exposure_excludes_end_trimmed_positions(tmp_path):
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>target\n' + 'A' * 100 + '\n')
    reads = [AlignedRead(str(i), 'target', 0, 100, [(j, 'A', 40) for j in range(100)], 'A' * 100)
             for i in range(10)]
    p = Params(error_model='auto', multiple_testing='none', end_trim_bp=10, strand_filter='off')
    _, stats = variants.call_variants(reads, str(fasta), 'target', p)
    assert stats['analysed_positions'] == 80
    assert stats['analysed_bases'] == 800
    res = pipeline.analyse(reads, str(fasta), p)
    per_read = res['per_read']
    assert sum(r['bases'] for r in per_read) == stats['analysed_bases']


def test_observed_and_fixed_window_frequencies_are_both_reported_for_partial_coverage(tmp_path):
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>target\n' + 'A' * 100 + '\n')
    # ten reads cover the first 50 of 100 bases; three carry one mutation (30%, not germline)
    reads = []
    for i in range(10):
        pairs = [(j, 'C' if j == 5 and i < 3 else 'A', 40) for j in range(50)]
        reads.append(AlignedRead(str(i), 'target', 0, 50, pairs, 'A' * 50))
    p = Params(error_model='auto', multiple_testing='none', strand_filter='off')
    res = pipeline.analyse(reads, str(fasta), p)
    assert len(res['variants']) == 1
    m = project.sample_metrics('s', res, 'A' * 100, 0, 99, 's', False)
    assert m['mutation_frequency_x1e3'] == pytest.approx(1e3 * 3 / 1000)          # fixed window
    assert m['mutation_frequency_observed_x1e3'] == pytest.approx(1e3 * 3 / 500)  # observed bases
    assert m['density_bases'] == 500


def test_indels_in_control_mode_say_they_used_bh():
    p = Params(multiple_testing='control')
    _called, stats = indels.call(SimpleNamespace(ref='A', events={}, span=lambda *a: 0), None, p)
    assert stats['indel_multiple_testing'] == 'bh'
    assert any('Benjamini-Hochberg' in n for n in stats['notes'])


# ---- third pass: F07, F10, F11, F12, F13, F14, F15 --------------------------

def _random_seq(n, seed):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


def _mutate(seq, sites, seed):
    rng = random.Random(seed)
    s = list(seq)
    for i in sites:
        s[i] = rng.choice([b for b in 'ACGT' if b != s[i]])
    return ''.join(s)


@pytest.fixture
def control_library(tmp_path):
    """A reference, a shorter paralogue differing at 15 sites, and two control
    FASTQs: clean reads with scattered errors, paralogue-allele reads that the
    full-length reference wins, and UMI tags in the read names."""
    ref = _random_seq(1000, 1)
    para_sites = list(range(260, 740, 32))
    paralogue = _mutate(ref, para_sites, 2)
    (tmp_path / 'ref.fa').write_text('>target\n' + ref + '\n')
    (tmp_path / 'para.fa').write_text('>paralogue\n' + paralogue[200:800] + '\n')
    files = []
    for f in range(2):
        lines = []
        for i in range(24):
            umi_tag = _random_seq(12, 100 + f * 10 + i // 4)    # six molecules of four reads
            if i % 6 == 5:
                seq = _mutate(paralogue, [400 + 7 * f], 50 + i)
            else:
                seq = _mutate(ref, [(37 * i + 11 * f) % 900 + 50, 500], 10 + i)
            lines += [f'@r{i}_{umi_tag}', seq, '+', '?' * len(seq)]   # Q30
        path = tmp_path / f'ctl{f}.fastq'
        path.write_text('\n'.join(lines) + '\n')
        files.append(str(path))
    return str(tmp_path / 'ref.fa'), [str(tmp_path / 'para.fa')], files


@pytest.mark.parametrize('extra', [
    {},                                          # paralogue-allele filter
    {'ignore_qualities': True},
    {'umi': 'header', 'umi_min_reads': 3, 'umi_max_edits': 0},
])
def test_standalone_null_fit_prepares_controls_as_samples_are(control_library, extra):
    """F07: null_model.fit sees the reads prepare_reads would give, phased as controls."""
    reference, contaminants, files = control_library
    p = Params(error_model='qcal', phase_reads=False, contaminant_min_alleles=10, **extra)
    null, stats = null_model.fit(files, reference, contaminants, p)
    prepared = [pipeline.prepare_reads(f, reference, contaminants, null_model.control_params(p))[0]
                for f in files]
    ref = dict(read_fasta(reference))['target']
    expect, _ = null_model.fit_aligned(prepared, ref, p, 0, len(ref) - 1, ref_name='target')
    assert stats['control_reads_used'] == sum(len(a) for a in prepared)
    assert dict(null) == dict(expect)
    assert null.qtail == expect.qtail and null.site_quals == expect.site_quals
    if not extra:
        assert stats['control_reads_used'] == 40           # the 8 paralogue reads are gone
    elif 'umi' in extra:
        assert stats['control_reads_used'] == 12           # six molecules per file
    else:
        assert {q for qs in null.site_quals.values() for q in qs} == {FASTA_DEFAULT_Q}


def _gated_library(n=20, length=50):
    """Reads on an all-A reference: every 10th base is Q30, the rest Q10. Q10
    bases err to C at 20%, Q30 bases at a fixed 2%."""
    rng = random.Random(3)
    reads = []
    for i in range(n):
        pairs = []
        for j in range(length):
            q = 30 if (i * length + j) % 10 == 0 else 10
            err = rng.random() < (0.02 if q == 30 else 0.2)
            pairs.append((j, 'C' if err else 'A', q))
        reads.append(AlignedRead(str(i), 'target', 0, length, pairs, 'A' * length,
                                 1 if i % 2 else -1))
    return reads


def test_null_denominator_counts_only_bases_passing_the_quality_gate(tmp_path):
    """F10: gating the pileup at Q20 is the same as fitting on the Q>=20 bases alone,
    and the decoys' base total matches the sample's."""
    lib = [_gated_library(), _gated_library(15)]
    eligible = [[AlignedRead(r.qname, r.ctg, r.r_st, r.r_en, [c for c in r.pairs if c[2] >= 20],
                             r.seq, r.strand) for r in f] for f in lib]
    ref = 'A' * 50
    gated = Params(error_model='qcal', min_base_quality=20, germline_freq_threshold=0.9)
    plain = replace(gated, min_base_quality=0)
    a, _ = null_model.fit_aligned(lib, ref, gated, 0, 49, calibrate=False)
    b, _ = null_model.fit_aligned(eligible, ref, plain, 0, 49, calibrate=False)
    assert dict(a) == dict(b) and a.strand_rates == b.strand_rates
    assert a.qtail == b.qtail and a.site_depth == b.site_depth
    # the quality tail above the gate is conditional on passing it: about 2%, not
    # the 2% x 10% share an ungated denominator gives
    k = sum(c[1] == 'C' for f in eligible for r in f for c in r.pairs if r.strand == 1)
    n = sum(1 for f in eligible for r in f for c in r.pairs if r.strand == 1)
    assert a.qtail[('A', 'C', 1)][20] == pytest.approx((k + 0.5) / (n + 1.0))
    ungated, _ = null_model.fit_aligned(lib, ref, plain, 0, 49, calibrate=False)
    assert ungated.qtail[('A', 'C', 1)][20] < a.qtail[('A', 'C', 1)][20] / 5
    # decoy exposure: the same bases the sample's control threshold counts
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>target\n' + ref + '\n')
    p = replace(gated, end_trim_bp=0, multiple_testing='none')
    _, stats = variants.call_variants(lib[0], str(fasta), 'target', p)
    assert null_model.window_bases(lib[0], 0, 49, variants.gate_quality(p)) == stats['analysed_bases']


INDEL_REF = ('GATCCATGCAAATTTGGGCCCATATCGCGATCGATCGTAGCTAGCATGCAGTCAGTCAGGATCCTAGC' * 4)


def _indel_reads(n, events, offset=0, start=0, end=None, q=40):
    end = len(INDEL_REF) if end is None else end
    return [AlignedRead(f'r{i}', 'ref', start + offset, end + offset, [], '', 1 if i % 2 else -1,
                        [(pos + offset, op, x, q) for pos, op, x in events(i)])
            for i in range(n)]


def _indel_noise(i):
    return [(60, 'D', 1)] if i % 20 == 0 else [(100 + i % 50, 'D', 1)] if i % 7 == 0 else []


def test_indel_class_rates_use_spanning_coverage_inside_the_window():
    """F11: padding the reference outside an unchanged window changes nothing, and
    class opportunities are the reads that actually span each event."""
    lo, hi = 40, 180
    pad = 'T' * 3 + 'GACGTCAGCTAGCATCGATCGTACG' * 20 + 'C' * 3
    ref2 = pad + INDEL_REF + pad
    off = len(pad)
    piles1 = [indels.Pileup(_indel_reads(200, _indel_noise, start=0 if k else 90), INDEL_REF, lo, hi)
              for k in range(2)]
    piles2 = [indels.Pileup(_indel_reads(200, _indel_noise, off, start=0 if k else 90), ref2,
                            lo + off, hi + off) for k in range(2)]
    n1, n2 = indels.IndelNull(piles1, INDEL_REF), indels.IndelNull(piles2, ref2)
    assert n1.cls_rate == n2.cls_rate
    sample = lambda i: _indel_noise(i) + ([(120, 'D', 12)] if i < 10 else [])   # noqa: E731
    p = Params(end_trim_bp=0, indel_min_reads=3)
    c1 = indels.call(indels.Pileup(_indel_reads(300, sample), INDEL_REF, lo, hi), n1, p)[0]
    c2 = indels.call(indels.Pileup(_indel_reads(300, sample, off), ref2, lo + off, hi + off), n2, p)[0]
    assert [(v.pos, v.p_value) for v in c1] == [(v.pos - off, v.p_value) for v in c2]
    # brute force: for each class, reads covering the base before and after the event
    opp = {}
    for pl in piles1:
        for pos0 in range(lo, hi + 1):
            for kind in ('del', 'ins'):
                for ln in (1, 2, 3, 4):
                    last = pos0 + (ln if kind == 'del' else 0)
                    spans = sum(1 for st, en in zip(pl.starts, pl.ends, strict=True)
                                if st <= pos0 - 1 and en > last)
                    c = indels.context(INDEL_REF, pos0, kind, 'N' * ln)
                    opp[c] = opp.get(c, 0) + spans
    k = {}
    for pl in piles1:
        for key, ev in pl.events.items():
            c = indels.context(INDEL_REF, *key)
            k[c] = k.get(c, 0) + len(ev)
    assert n1.cls_rate == {c: (k.get(c, 0) + 0.5) / (n + 1.0) for c, n in opp.items() if n}


@pytest.mark.parametrize('ancestral', [False, True])
def test_indel_frequency_filter_honours_ancestral_reference(ancestral):
    """F11: with an ancestral reference a fixed insertion is tested, not removed:
    called when the controls lack it, not called when they carry it too. As in a
    project, the controls' rates are inflated by the leave-one-out multiplier."""
    lo, hi = 0, len(INDEL_REF) - 1
    p = Params(end_trim_bp=0, indel_min_reads=3, ancestral_reference=ancestral)
    fixed = lambda i: [(150, 'I', 'GG')] + _indel_noise(i)                       # noqa: E731
    clean = [indels.Pileup(_indel_reads(200, _indel_noise), INDEL_REF, lo, hi) for _ in range(2)]
    carrying = [indels.Pileup(_indel_reads(200, fixed), INDEL_REF, lo, hi) for _ in range(2)]
    case = indels.Pileup(_indel_reads(200, fixed), INDEL_REF, lo, hi)
    key = next(k for k in case.events if k[1] == 'ins')
    for ctrls, expect_call in ((clean, ancestral), (carrying, False)):
        null = indels.IndelNull(ctrls, INDEL_REF, ancestral)
        null.qscale = indels.calibrate(ctrls, INDEL_REF, p)[0]
        called, stats = indels.call(case, null, p)
        assert any(v.kind == 'ins' and v.pos == key[0] + 1 for v in called) == expect_call
        assert stats['indel_germline_removed'] == (0 if ancestral else 1)
    # the fixed event is the controls' own site rate, never class error
    null = indels.IndelNull(carrying, INDEL_REF, True)
    assert null.site_k[key] == 400
    assert indels.IndelNull(carrying, INDEL_REF, False).cls_rate == null.cls_rate


def _run_samples(spec):
    return [{'name': n, 'run': r, 'path': None} for n, r in spec]


def test_each_case_is_checked_against_the_control_runs():
    """F12: controls from A, cases from A and B: the B case is flagged on its own."""
    samples = _run_samples([('K1', 'A'), ('K2', 'A'), ('W1', 'A'), ('W2', 'B'), ('W3', ['A', 'B'])])
    rc = project.run_compatibility(samples, {'K1', 'K2'})
    status = {n: c['status'] for n, c in rc['samples'].items()}
    assert status == {'W1': 'same run', 'W2': 'different run', 'W3': 'different run'}
    per_sample = [{'sample': n, 'warnings': ''} for n in ('W1', 'W2', 'W3')]
    notes = project.run_notes(rc, per_sample)
    assert any(n.startswith('CALIBRATION NOT SUPPORTED for W2 (run B), W3 (run B)') for n in notes)
    assert [bool(m['warnings']) for m in per_sample] == [False, True, True]
    # mixed control runs are named as a pooled, unmatched model
    rc = project.run_compatibility(_run_samples([('K1', 'A'), ('K2', 'B'), ('W1', 'A')]), {'K1', 'K2'})
    assert rc['samples']['W1']['status'] == 'same run'
    assert any('span sequencing runs A, B' in n for n in project.run_notes(rc, []))


def test_unknown_runs_are_not_reported_as_same_run(tmp_path):
    """F12: run IDs from read headers, across every file of a sample; an
    unidentified file leaves the case unverified rather than matched."""
    def fq(name, header):
        f = tmp_path / f'{name}.fastq'
        f.write_text(f'@{header}\nACGT\n+\nIIII\n')
        return str(f)
    samples = [
        {'name': 'K1', 'path': fq('k1', '0_RUNA_1')},
        {'name': 'W1', 'path': [fq('w1a', '0_RUNA_1'), fq('w1b', '0_RUNB_2')]},
        {'name': 'W2', 'path': [fq('w2a', '0_RUNA_1'), fq('w2b', 'read7')]},
        {'name': 'W3', 'path': fq('w3', 'read8')},
        {'name': 'W4', 'path': fq('w4', '0_RUNA_9')},
    ]
    rc = project.run_compatibility(samples, {'K1'})
    status = {n: c['status'] for n, c in rc['samples'].items()}
    assert status == {'W1': 'different run', 'W2': 'unverified', 'W3': 'unverified', 'W4': 'same run'}
    assert any('W2, W3' in n and 'not verified' in n for n in project.run_notes(rc, []))
    # nothing identifiable anywhere: no claim either way, and no note
    rc = project.run_compatibility([{'name': 'K', 'path': fq('k', 'a')}, {'name': 'W', 'path': fq('w', 'b')}],
                                   {'K'})
    assert rc['samples']['W']['status'] == 'unverified' and project.run_notes(rc, []) == []


def test_clone_summary_uses_the_callers_quality_gate():
    """F13: C is accepted elsewhere, but in this two-read lineage only a Q10 base
    supports it; at a Q30 gate it neither carries nor covers there."""
    j = 'ACGTACGTACGTACGTACGT'
    reads = [AlignedRead(n, 'target', 0, 20, [(5, b, q)] + [(i, 'A', 40) for i in range(6, 20)],
                         j + 'A' * 20, s, [], 20, 40)
             for n, b, q, s in (('lo', 'C', 10, 1), ('hi', 'A', 40, -1))]
    v = SimpleNamespace(pos=6, alt='C')
    base = dict(junction_ref_start=1, junction_length=20, junction_max_edits=0)
    rows = lambda minq: clones.call_clones(reads, [v], Params(min_base_quality=minq, **base),  # noqa: E731
                                           0, 19)['rows']
    assert rows(0)[0]['mutations'] == 1         # the historical, ungated count
    assert rows(30)[0]['mutations'] == 0
    # one eligible read at a tie: half of the covering reads carry it
    reads[1].pairs[0] = (5, 'C', 40)
    assert rows(30)[0]['mutations'] == 1


def test_clone_pair_sampling_is_bounded_and_exact_on_small_sets():
    """F14: pair indices are drawn without building every pair."""
    called = {(0, 'C'), (1, 'C')}
    def lib(n):
        rs = {f'r{i}': SimpleNamespace(pairs=[(i % 2, 'C', 40), (2, 'A', 40)]) for i in range(n)}
        return rs, {f'r{i}': i % 3 for i in range(n)}
    reads, clone_of = lib(60)
    full = clones.sharing(reads, clone_of, called, max_pairs=10 ** 6)
    within, between = [], []
    for a, b in itertools.combinations(reads, 2):
        A = {(p, x) for p, x, _ in reads[a].pairs if (p, x) in called}
        B = {(p, x) for p, x, _ in reads[b].pairs if (p, x) in called}
        (within if clone_of[a] == clone_of[b] else between).append(len(A & B) / len(A | B))
    assert full['pairs_within'] + full['pairs_between'] == 60 * 59 // 2
    assert full['within'] == round(sum(within) / len(within), 4)
    assert full['between'] == round(sum(between) / len(between), 4)
    capped = clones.sharing(reads, clone_of, called, max_pairs=100)
    assert capped == clones.sharing(reads, clone_of, called, max_pairs=100)   # fixed seed
    assert capped['pairs_within'] + capped['pairs_between'] == 100
    big, big_clones = lib(3000)                      # 4.5 million pairs
    tracemalloc.start()
    out = clones.sharing(big, big_clones, called, max_pairs=1000)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert out['pairs_within'] + out['pairs_between'] == 1000
    assert peak < 8 * 2 ** 20                        # all pairs would take hundreds of MB


def test_umi_members_survive_pickling_and_worker_processes():
    """F15: consensus reads keep the reads they came from, serial or parallel."""
    plain = AlignedRead('c', 'target', 0, 3, [(0, 'A', 40), (1, 'C', 30), (2, 'G', 20)], 'ACG')
    odd = AlignedRead('d', 'target', 0, 1, [(0, 'A', 40.5)], 'A')     # the unpacked path
    for r in (plain, odd):
        r.members = ['a', 'b', 'c']
        back = pickle.loads(pickle.dumps(r))
        assert back == r and back.members == ['a', 'b', 'c']
    assert not hasattr(pickle.loads(pickle.dumps(AlignedRead('e', 't', 0, 0, [], ''))), 'members')
    serial = pmap(copy.copy, [(plain,), (odd,)], 1)
    parallel = pmap(copy.copy, [(plain,), (odd,)], 2)
    assert [r.members for r in serial] == [r.members for r in parallel] == [['a', 'b', 'c']] * 2
