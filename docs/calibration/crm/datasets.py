"""Datasets for the Cell Reports Methods benchmarks: delivered reads and public deposits."""
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
REF = REPO / "src/shmtrace/references"
DEP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/"
                     "ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data")
PUB = Path.home() / "Nanoshm/public_data"
IRCM = PUB / "PRJNA1481264"

JH4 = {"WT2": "DHJGKB_2_Wildtype_2", "WT3": "DHJGKB_3_Wildtype_3", "WT4": "DHJGKB_4_Wildtype_4",
       "KO1": "DHJGKB_5_AID_KO_1", "KO2": "DHJGKB_6_AID_KO_2", "KO3": "DHJGKB_7_AID_KO_3"}
N114 = {"C1": "75K3B3_6_C1_AID_10.06.24", "D3": "75K3B3_7_D3_AID_10.01.24", "D5": "75K3B3_10_D5_AID_6.28.24_",
        "C3": "75K3B3_1_C3_Empty_March_2024_1", "C4": "75K3B3_3_C4_Empty_March_2024",
        "D1": "75K3B3_4_D1_Empty_March_2024"}

DATASETS = {
    "jh4": {"preset": "jh4_intron", "reference": str(REF / "mouse_jh4_intron.fasta"),
            "contaminants": [str(REF / "mouse_jh4_contaminant.fasta")], "depth": 500,
            "samples": {k: {"fastq": str(DEP / "02_JH4_ONT/Original_FASTQ" / f"{v}.fastq"),
                            "control": k.startswith("KO")} for k, v in JH4.items()}},
    "n114": {"preset": "n114_vregion", "reference": str(REF / "n114_vregion.fasta"), "depth": 500,
             "samples": {k: {"fastq": str(DEP / "03_N114_ONT/Original_FASTQ" / f"{v}.fastq"),
                             "control": k in ("C3", "C4", "D1")} for k, v in N114.items()}},
}


def ircm(locus, groups, reference, preset, depth=2000):
    """IRCM public ONT runs (PRJNA1481264): `groups` maps a sample to (library name, control)."""
    return {"preset": preset, "reference": str(reference), "depth": depth,
            "samples": {k: {"fastq": str(IRCM / "sub" / f"{lib}.fastq"), "control": c}
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
# Pim1 again at 4,000 reads: an AID off-target near the detection floor, where depth should matter
DATASETS["pim1_deep"] = {**DATASETS["pim1"], "depth": 4000}
