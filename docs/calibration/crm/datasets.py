"""Datasets for the Cell Reports Methods benchmarks: delivered reads and public deposits."""
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
REF = REPO / "src/shmtrace/references"
DEP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/"
                     "ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data")
PUB = Path.home() / "Nanoshm/public_data"
IRCM = PUB / "PRJNA1481264"

# N114 reads: ~/Nanoshm/n114_fastq links the delivered libraries; empty-vector subclone C3 is
# its second library (C3_Empty_March_2024_2), the first being noisier than its peers (CRM paper)
N114_DIR = Path.home() / "Nanoshm/n114_fastq"
JH4 = {"WT2": "DHJGKB_2_Wildtype_2", "WT3": "DHJGKB_3_Wildtype_3", "WT4": "DHJGKB_4_Wildtype_4",
       "KO1": "DHJGKB_5_AID_KO_1", "KO2": "DHJGKB_6_AID_KO_2", "KO3": "DHJGKB_7_AID_KO_3"}
N114 = {"C1": "75K3B3_6_C1_AID_10.06.24", "D3": "75K3B3_7_D3_AID_10.01.24", "D5": "75K3B3_10_D5_AID_6.28.24_",
        "C3": "75K3B3_2_C3_Empty_March_2024_2", "C4": "75K3B3_3_C4_Empty_March_2024",
        "D1": "75K3B3_4_D1_Empty_March_2024"}

DATASETS = {
    "jh4": {"preset": "jh4_intron", "reference": str(REF / "mouse_jh4_intron.fasta"),
            "contaminants": [str(REF / "mouse_jh4_contaminant.fasta")], "depth": 500,
            "samples": {k: {"fastq": str(DEP / "02_JH4_ONT/Original_FASTQ" / f"{v}.fastq"),
                            "control": k.startswith("KO")} for k, v in JH4.items()}},
    "n114": {"preset": "n114_vregion", "reference": str(REF / "n114_vregion.fasta"), "depth": 500,
             "samples": {k: {"fastq": str(N114_DIR / f"{v}.fastq"),
                             "control": k in ("C3", "C4", "D1")} for k, v in N114.items()}},
}

# N114 sensitivity sets: the second C3 library alone (now the default, so n114_c3b = n114), and all
# four empty-vector libraries as controls
HYB = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/Hybridoma manuscript/"
                     "75K3B3_fastq")
for _key, _ctl in (("n114_c3b", {"C3_2": "75K3B3_2_C3_Empty_March_2024_2", "C4": N114["C4"], "D1": N114["D1"]}),
                   ("n114_4ctl", {"C3": "75K3B3_1_C3_Empty_March_2024_1", "C3_2": "75K3B3_2_C3_Empty_March_2024_2", "C4": N114["C4"],
                                  "D1": N114["D1"]})):
    _s = {k: {"fastq": str(N114_DIR / f"{N114[k]}.fastq"), "control": False}
          for k in ("C1", "D3", "D5")}
    _s.update({k: {"fastq": str(HYB / f"{v}.fastq"), "control": True} for k, v in _ctl.items()})
    DATASETS[_key] = {**DATASETS["n114"], "samples": _s}


def ircm(locus, groups, reference, preset, depth=2000, subdir="sub"):
    """IRCM public ONT runs (PRJNA1481264): `groups` maps a sample to (library name, control).

    `subdir` picks the subsample: "sub" is 6,000 reads per library, "sub_deep" is 60,000, for the weak
    off-targets where 6,000 leaves nothing above the control background.
    """
    return {"preset": preset, "reference": str(reference), "depth": depth,
            "samples": {k: {"fastq": str(IRCM / subdir / f"{lib}.fastq"), "control": c}
                        for k, (lib, c) in groups.items()}}


# IRCM PRJNA1481264: AID-competent (DMSO-treated or WT) against Aicda-/- libraries, 2,000 reads each
_SU = {"min_align_len": 450, "max_total_len": 1000}
_PIM = {"min_align_len": 600, "max_total_len": 1400}
DATASETS["su_lps"] = {**ircm("Igh_Su", {f"WT{i}": (f"Igh_Su_LPS_DMSO_r{i}", False) for i in (1, 2, 3)}
                              | {"KO1": ("Igh_Su_LPS_AID_KO", True), "KO2": ("Igh_Su_LPS_AID_KO_r2", True),
                                 "KO3": ("Igh_Su_LPS_AID_KO_r3", True)},
                              IRCM / "references/Igh_Su.fasta", "generic_amplicon"), "overrides": _SU}
DATASETS["su_ch12"] = {**ircm("Igh_Su", {f"WT{i}": (f"Igh_Su_CH12_WT_r{i}", False) for i in (1, 2, 3, 4)}
                               | {f"KO{i}": (f"Igh_Su_CH12_AID_KO_r{i}", True) for i in (1, 2, 3)},
                               IRCM / "references/Igh_Su.fasta", "generic_amplicon"), "overrides": _SU}
DATASETS["pim1"] = {**ircm("Pim1", {f"WT{i}": (f"Pim1_LPS_DMSO_r{i}", False) for i in (1, 2, 3)}
                           | {"KO1": ("Pim1_LPS_AID_KO", True), "KO2": ("Pim1_LPS_AID_KO_r2", True),
                              "KO3": ("Pim1_LPS_AID_KO_r3", True)},
                           IRCM / "references/Pim1.fasta", "generic_amplicon"), "overrides": _PIM}
def _ref_gates(path, floor=0.6, ceil=2.0):
    """Length gates from the consensus reference itself: an alignment must cover `floor` of it, and a read
    longer than `ceil` times it is a concatemer. The second-pass loci have no published amplicon size."""
    try:
        n = len("".join(x.strip() for x in open(path) if not x.startswith(">")))
    except OSError:
        return {}
    return {"min_align_len": int(floor * n), "max_total_len": int(ceil * n)}


# Second pass over PRJNA1481264: the loci and cell contexts not used above, each with its own Aicda-/-
# controls. Cd83 and Ly6e are AID off-targets in LPS blasts and in iGB cells; OT1-OT3 are the CH12
# off-targets, which add an MLLT1/MLLT3 double knockout expected to lower mutation; the iMLLT arms are the
# MLLT inhibitor. Length gates come from each consensus reference (_ref_gates). Where a cell context has
# only one Aicda-/- library of its own, the other contexts' knockout libraries of the SAME locus join it as
# controls: their per-position error profiles agree as closely as two controls of one context do (Pearson
# r 0.76-0.86 against 0.76-0.77), so they are interchangeable as decoys. OT3 has only one knockout library
# in the whole project and cannot be calibrated.
for _key, _locus, _case, _ctl in (
        ("cd83_lps", "Cd83", {f"WT{i}": f"Cd83_LPS_DMSO_r{i}" for i in (1, 2, 3)}
         | {f"iMLLT{i}": f"Cd83_LPS_iMLLT_r{i}" for i in (1, 2, 3)},
         {"KO1": "Cd83_LPS_AID_KO", "KO2": "Cd83_LPS_AID_KO_r2", "KO3": "Cd83_LPS_AID_KO_r3"}),
        ("cd83_igb", "Cd83", {f"WT{i}": f"Cd83_iGB_DMSO_r{i}" for i in (1, 2)}
         | {f"iMLLT{i}": f"Cd83_iGB_iMLLT_r{i}" for i in (1, 2)},
         {"KO1": "Cd83_iGB_AID_KO", "KO2": "Cd83_LPS_AID_KO", "KO3": "Cd83_LPS_AID_KO_r2",
          "KO4": "Cd83_LPS_AID_KO_r3"}),
        ("ly6e_lps", "Ly6e", {f"WT{i}": f"Ly6e_LPS_DMSO_r{i}" for i in (1, 2, 3)}
         | {f"iMLLT{i}": f"Ly6e_LPS_iMLLT_r{i}" for i in (1, 2, 3)},
         {"KO1": "Ly6e_LPS_AID_KO", "KO2": "Ly6e_LPS_AID_KO_r2", "KO3": "Ly6e_LPS_AID_KO_r3"}),
        ("ly6e_igb", "Ly6e", {f"WT{i}": f"Ly6e_iGB_DMSO_r{i}" for i in (1, 2)}
         | {f"iMLLT{i}": f"Ly6e_iGB_iMLLT_r{i}" for i in (1, 2)},
         {"KO1": "Ly6e_iGB_AID_KO", "KO2": "Ly6e_LPS_AID_KO", "KO3": "Ly6e_LPS_AID_KO_r2",
          "KO4": "Ly6e_LPS_AID_KO_r3"}),
        ("pim1_igb", "Pim1", {f"WT{i}": f"Pim1_iGB_DMSO_r{i}" for i in (1, 2)}
         | {f"iMLLT{i}": f"Pim1_iGB_iMLLT_r{i}" for i in (1, 2)},
         {"KO1": "Pim1_iGB_AID_KO", "KO2": "Pim1_LPS_AID_KO", "KO3": "Pim1_LPS_AID_KO_r2",
          "KO4": "Pim1_LPS_AID_KO_r3"}),
        ("ot1", "OT1", {f"WT{i}": f"OT1_CH12_WT_r{i}" for i in (1, 2, 3)}
         | {f"dKO{i}": f"OT1_CH12_dKO_r{i}" for i in (1, 2, 3)},
         {f"KO{i}": f"OT1_CH12_AID_KO_r{i}" for i in (1, 2, 3)}),
        ("ot2", "OT2", {f"WT{i}": f"OT2_CH12_WT_r{i}" for i in (1, 2, 3)}
         | {f"dKO{i}": f"OT2_CH12_dKO_r{i}" for i in (1, 2, 3)},
         {f"KO{i}": f"OT2_CH12_AID_KO_r{i}" for i in (1, 2, 3)}),
        ("ot3", "OT3", {f"WT{i}": f"OT3_CH12_WT_r{i}" for i in (1, 2, 3)}
         | {f"dKO{i}": f"OT3_CH12_dKO_r{i}" for i in (1, 2, 3)}, {"KO1": "OT3_CH12_AID_KO"})):
    DATASETS[_key] = {**ircm(_locus, {k: (v, False) for k, v in _case.items()}
                             | {k: (v, True) for k, v in _ctl.items()},
                             IRCM / f"references/{_locus}.fasta", "generic_amplicon"),
                      "overrides": _ref_gates(IRCM / f"references/{_locus}.fasta")}

# The off-targets that gave nothing at 6,000 reads, re-run on 60,000-read subsamples (sub_deep)
for _k in ("cd83_lps", "cd83_igb", "ly6e_lps"):
    DATASETS[f"{_k}_deep"] = {**DATASETS[_k], "depth": 30000,
                              "samples": {n: {**s, "fastq": s["fastq"].replace("/sub/", "/sub_deep/")}
                                          for n, s in DATASETS[_k]["samples"].items()}}

# Pim1 again at 4,000 reads: an AID off-target near the detection floor, where depth should matter
DATASETS["pim1_deep"] = {**DATASETS["pim1"], "depth": 4000}
