# SHM comparative analysis — CONTROL vs CASE

Preset **jh4_intron**, window (101, 629). Control = 3 file(s), 1500 reads; Case = 3 file(s), 1498 reads.

Frequencies are base-composition corrected (unique mutations). Stars: * p<0.05, ** p<0.01, *** p<0.001.

## Mutation profile

| Metric | Control | Case | Δ | p |
|---|---:|---:|:--:|---|
| Overall mutation frequency | 9.14e-02 | 1.26e-02 | ↓ | *** |
| A:T mutations | 1.03e-01 | 9.23e-03 | ↓ | *** |
| G:C mutations | 7.55e-02 | 1.71e-02 | ↓ | *** |
| G:C transitions (C→T, G→A) | 1.18e-01 | 2.63e-02 | ↓ | *** |

## Motif targeting (hot/cold-spots)

| Motif class | Control | Case | Δ | p |
|---|---:|---:|:--:|---|
| AID_hotspot (WRC/GYW/RGYW/WRCY) | 4.23e-01 | 4.50e-02 | ↓ | *** |
| polymerase_eta_hotspot (WA/TW) | 4.12e-01 | 3.29e-02 | ↓ | *** |
| coldspot (SYC/GRS) | 1.08e-01 | 2.45e-02 | ↓ | ** |

## Strand bias (AID hotspot top vs bottom strand)

Control WRC/GYW = 28/37; Case WRC/GYW = 3/3. similar WRC and GYW mutation counts imply no strand bias.