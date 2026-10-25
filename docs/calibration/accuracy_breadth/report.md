# Accuracy beyond J꜀4 and N114: report

Estimate / truth: mean over held-out controls and draws, with the range in brackets (1.00 = exact). Null rows give the false calls and the unique density estimated where the truth is zero. Recall is of mutations carried by two or more analysed reads. Truth / error: the planted mutation frequency and the mean estimate minus it, per 10³ bp; at the loads of these studies a library holds a few dozen planted mutations, and the ratios swing with them.

## Breadth: every public set with three or more mutation-free controls

| Set | Scenario | Depth | Runs | Truth / error | Mutation frequency | Unique density | Combined density | False calls | Recall, 2+ reads |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Illumina J꜀4 (PRJNA973112) | low | 1000 | 8 | 0.52 / -0.09 | 0.77 (0.00–1.06) | 0.32 (0.00–0.67) | 0.67 (0.00–1.52) | 0.4 (max 1) | 0.50 |
| Illumina J꜀4 (PRJNA973112) | low | 5000 | 8 | 0.80 / -0.15 | 0.34 (0.00–0.99) | 0.04 (0.00–0.17) | 0.56 (0.00–1.63) | 0.0 (max 0) | 0.21 |
| Illumina J꜀4 (PRJNA973112) | null | 1000 | 8 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina J꜀4 (PRJNA973112) | null | 5000 | 8 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina J꜀4 (PRJNA973112) | polyclonal | 1000 | 8 | 3.29 / -0.49 | 0.68 (0.18–0.97) | 0.35 (0.11–0.58) | 0.50 (0.11–0.87) | 0.2 (max 1) | 0.56 |
| Illumina J꜀4 (PRJNA973112) | polyclonal | 5000 | 8 | 3.09 / -0.62 | 0.78 (0.47–0.94) | 0.16 (0.05–0.30) | 0.58 (0.18–0.99) | 0.0 (max 0) | 0.43 |
| Illumina CH12 VDJ (PRJNA918596) | low | 1000 | 8 | 0.02 / -0.01 | 0.25 (0.00–1.00) | 0.21 (0.00–1.00) | 20.71 (7.42–29.66) | 0.0 (max 0) | 0.75 |
| Illumina CH12 VDJ (PRJNA918596) | low | 5000 | 8 | 0.03 / -0.00 | 1.87 (0.00–8.74) | 0.16 (0.00–0.50) | 5.49 (3.43–9.09) | 0.5 (max 2) | 0.83 |
| Illumina CH12 VDJ (PRJNA918596) | null | 1000 | 8 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina CH12 VDJ (PRJNA918596) | null | 5000 | 8 | – | – | 0.04 (est) | – | 0.5 (max 2) | – |
| Illumina CH12 VDJ (PRJNA918596) | subclonal | 1000 | 8 | 0.12 / -0.03 | 0.72 (0.00–1.18) | 0.22 (0.00–0.50) | 4.44 (3.27–7.22) | 0.0 (max 0) | 0.71 |
| Illumina CH12 VDJ (PRJNA918596) | subclonal | 5000 | 8 | 0.21 / +0.00 | 0.93 (0.47–1.40) | 0.16 (0.05–0.41) | 1.27 (0.99–1.70) | 1.1 (max 6) | 0.81 |
| Illumina Ramos VDJ (PRJNA820976) | low | 1000 | 6 | 0.11 / -0.11 | 0.00 (0.00–0.00) | 0.00 (0.00–0.00) | 2.59 (0.00–11.01) | 0.0 (max 0) | 0.00 |
| Illumina Ramos VDJ (PRJNA820976) | low | 5000 | 6 | 0.09 / -0.09 | 0.00 (0.00–0.00) | 0.00 (0.00–0.00) | 4.17 (3.11–8.02) | 0.0 (max 0) | 0.00 |
| Illumina Ramos VDJ (PRJNA820976) | null | 1000 | 6 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina Ramos VDJ (PRJNA820976) | null | 5000 | 6 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina Ramos VDJ (PRJNA820976) | subclonal | 1000 | 6 | 0.25 / -0.19 | 0.14 (0.00–0.87) | 0.05 (0.00–0.31) | 0.63 (0.07–1.44) | 0.0 (max 0) | 0.14 |
| Illumina Ramos VDJ (PRJNA820976) | subclonal | 5000 | 6 | 0.32 / -0.29 | 0.15 (0.00–0.88) | 0.02 (0.00–0.11) | 0.73 (0.57–0.87) | 0.0 (max 0) | 0.13 |
| ONT Sμ, LPS blasts | low | 500 | 6 | 0.40 / -0.02 | 0.93 (0.84–1.03) | 0.52 (0.23–1.00) | 1.02 (0.40–1.53) | 1.2 (max 3) | 0.89 |
| ONT Sμ, LPS blasts | low | 2000 | 6 | 0.19 / -0.02 | 0.88 (0.78–0.93) | 0.19 (0.11–0.31) | 1.02 (0.43–2.01) | 0.3 (max 1) | 0.73 |
| ONT Sμ, LPS blasts | null | 500 | 6 | – | – | 0.33 (est) | – | 1.0 (max 3) | – |
| ONT Sμ, LPS blasts | null | 2000 | 6 | – | – | 0.03 (est) | – | 0.3 (max 1) | – |
| ONT Sμ, LPS blasts | polyclonal | 500 | 6 | 1.36 / -0.05 | 0.96 (0.91–1.01) | 0.73 (0.52–0.93) | 0.95 (0.77–1.12) | 1.7 (max 3) | 0.89 |
| ONT Sμ, LPS blasts | polyclonal | 2000 | 6 | 1.08 / -0.05 | 0.95 (0.91–0.98) | 0.33 (0.26–0.39) | 0.98 (0.75–1.25) | 0.5 (max 1) | 0.81 |
| ONT Sμ, CH12 | low | 500 | 6 | 0.37 / +0.03 | 1.93 (0.00–11.61) | 0.13 (0.00–0.76) | 0.65 (0.41–0.87) | 0.3 (max 2) | 0.15 |
| ONT Sμ, CH12 | low | 2000 | 6 | 0.31 / -0.03 | 0.91 (0.81–0.97) | 0.23 (0.12–0.38) | 1.00 (0.89–1.10) | 0.2 (max 1) | 0.79 |
| ONT Sμ, CH12 | null | 500 | 6 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| ONT Sμ, CH12 | null | 2000 | 6 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| ONT Sμ, CH12 | subclonal | 500 | 6 | 1.42 / +0.80 | 1.49 (0.86–2.48) | 0.76 (0.54–0.96) | 0.92 (0.81–0.99) | 5.8 (max 13) | 0.91 |
| ONT Sμ, CH12 | subclonal | 2000 | 6 | 1.36 / -0.06 | 0.95 (0.91–1.02) | 0.46 (0.25–0.77) | 0.89 (0.86–0.93) | 5.8 (max 15) | 0.86 |
| ONT Pim1 | low | 500 | 6 | 0.30 / -0.01 | 0.95 (0.88–1.03) | 0.59 (0.22–1.00) | 0.99 (0.46–1.49) | 0.7 (max 2) | 0.90 |
| ONT Pim1 | low | 2000 | 6 | 0.12 / -0.01 | 0.86 (0.73–0.93) | 0.22 (0.14–0.26) | 0.97 (0.79–1.15) | 0.0 (max 0) | 0.87 |
| ONT Pim1 | null | 500 | 6 | – | – | 0.08 (est) | – | 0.3 (max 1) | – |
| ONT Pim1 | null | 2000 | 6 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| ONT OT1 (CH12) | low | 500 | 6 | 0.27 / -0.03 | 0.71 (0.00–0.97) | 0.38 (0.00–0.65) | 1.39 (0.00–3.49) | 0.5 (max 1) | 0.62 |
| ONT OT1 (CH12) | low | 2000 | 6 | 0.22 / -0.03 | 0.85 (0.73–0.93) | 0.20 (0.14–0.26) | 0.89 (0.16–2.14) | 0.0 (max 0) | 0.70 |
| ONT OT1 (CH12) | null | 500 | 6 | – | – | 0.07 (est) | – | 0.3 (max 1) | – |
| ONT OT1 (CH12) | null | 2000 | 6 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| ONT OT2 (CH12) | low | 500 | 6 | 0.27 / -0.02 | 0.87 (0.70–0.99) | 0.48 (0.33–0.76) | 0.99 (0.33–1.71) | 0.3 (max 1) | 0.85 |
| ONT OT2 (CH12) | low | 2000 | 6 | 0.16 / -0.01 | 0.94 (0.82–1.12) | 0.28 (0.14–0.42) | 1.11 (0.68–2.00) | 0.5 (max 3) | 0.90 |
| ONT OT2 (CH12) | null | 500 | 6 | – | – | 0.07 (est) | – | 0.3 (max 1) | – |
| ONT OT2 (CH12) | null | 2000 | 6 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |

## Sμ with controls from the other cell context

| Set | Scenario | Depth | Runs | Truth / error | Mutation frequency | Unique density | Combined density | False calls | Recall, 2+ reads |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ONT Sμ, LPS blasts | null | 2000 | 6 | – | – | 0.14 (est) | – | 1.7 (max 3) | – |
| ONT Sμ, LPS blasts | polyclonal | 2000 | 6 | 1.08 / -0.01 | 0.99 (0.97–1.01) | 0.47 (0.41–0.56) | 0.94 (0.76–1.13) | 2.8 (max 5) | 0.92 |
| ONT Sμ, CH12 | null | 2000 | 6 | – | – | 0.17 (est) | – | 2.0 (max 2) | – |
| ONT Sμ, CH12 | subclonal | 2000 | 6 | 1.36 / +1.45 | 2.08 (1.89–2.29) | 0.66 (0.55–0.79) | 0.87 (0.81–0.94) | 11.8 (max 19) | 0.87 |

## Realism: the planted side made less ideal, one step at a time

| Set | Scenario | Planted | Runs | Truth / error | Mutation frequency | Unique density | Combined density | False calls | Recall, 2+ reads |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| Illumina J꜀4 (PRJNA973112) | null | base | 4 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina J꜀4 (PRJNA973112) | null | pcr3 | 4 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina J꜀4 (PRJNA973112) | null | pcr10 | 4 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina J꜀4 (PRJNA973112) | null | all | 4 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina J꜀4 (PRJNA973112) | polyclonal | base | 4 | 3.83 / -0.72 | 0.77 (0.47–0.92) | 0.18 (0.05–0.30) | 0.64 (0.39–0.99) | 0.0 (max 0) | 0.48 |
| Illumina J꜀4 (PRJNA973112) | polyclonal | qmap | 4 | 3.83 / -0.72 | 0.77 (0.47–0.92) | 0.18 (0.05–0.30) | 0.64 (0.39–0.99) | 0.0 (max 0) | 0.48 |
| Illumina J꜀4 (PRJNA973112) | polyclonal | pcr3 | 4 | 1.67 / -0.79 | 0.43 (0.00–0.90) | 0.15 (0.00–0.37) | 0.50 (0.07–0.93) | 0.2 (max 1) | 0.18 |
| Illumina J꜀4 (PRJNA973112) | polyclonal | pcr10 | 4 | 2.01 / -0.61 | 0.72 (0.55–0.83) | 0.32 (0.18–0.54) | 0.51 (0.18–0.80) | 0.2 (max 1) | 0.34 |
| Illumina J꜀4 (PRJNA973112) | polyclonal | deep | 4 | 3.71 / -0.87 | 0.67 (0.43–0.88) | 0.13 (0.05–0.23) | 0.35 (0.11–0.55) | 0.2 (max 1) | 0.19 |
| Illumina J꜀4 (PRJNA973112) | polyclonal | all | 4 | 2.12 / -0.63 | 0.69 (0.34–0.90) | 0.19 (0.05–0.31) | 0.50 (0.09–0.78) | 0.2 (max 1) | 0.23 |
| Illumina CH12 VDJ (PRJNA918596) | null | base | 4 | – | – | 0.04 (est) | – | 0.5 (max 2) | – |
| Illumina CH12 VDJ (PRJNA918596) | null | pcr3 | 4 | – | – | 0.15 (est) | – | 1.8 (max 3) | – |
| Illumina CH12 VDJ (PRJNA918596) | null | pcr10 | 4 | – | – | 0.31 (est) | – | 3.5 (max 8) | – |
| Illumina CH12 VDJ (PRJNA918596) | null | all | 4 | – | – | 0.15 (est) | – | 1.8 (max 3) | – |
| Illumina CH12 VDJ (PRJNA918596) | subclonal | base | 4 | 0.13 / +0.01 | 0.82 (0.47–1.40) | 0.16 (0.05–0.41) | 1.28 (0.99–1.70) | 1.5 (max 6) | 0.76 |
| Illumina CH12 VDJ (PRJNA918596) | subclonal | qmap | 4 | 0.13 / +0.01 | 0.82 (0.47–1.40) | 0.16 (0.05–0.41) | 1.28 (0.99–1.70) | 1.5 (max 6) | 0.76 |
| Illumina CH12 VDJ (PRJNA918596) | subclonal | pcr3 | 4 | 0.19 / +0.02 | 1.23 (0.85–1.67) | 0.58 (0.30–0.71) | 3.03 (2.24–4.54) | 4.8 (max 9) | 0.36 |
| Illumina CH12 VDJ (PRJNA918596) | subclonal | pcr10 | 4 | 0.13 / +0.10 | 2.24 (0.61–6.00) | 1.03 (0.17–2.33) | 4.72 (3.15–5.94) | 6.2 (max 19) | 0.35 |
| Illumina CH12 VDJ (PRJNA918596) | subclonal | deep | 4 | 0.23 / -0.04 | 0.73 (0.40–1.11) | 0.14 (0.03–0.27) | 1.14 (1.01–1.28) | 1.0 (max 3) | 0.26 |
| Illumina CH12 VDJ (PRJNA918596) | subclonal | all | 4 | 0.17 / +0.09 | 1.80 (0.61–4.69) | 0.74 (0.21–1.63) | 2.30 (1.60–3.88) | 8.5 (max 21) | 0.38 |
| Illumina Ramos VDJ (PRJNA820976) | null | base | 3 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina Ramos VDJ (PRJNA820976) | null | pcr3 | 3 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina Ramos VDJ (PRJNA820976) | null | pcr10 | 3 | – | – | 0.22 (est) | – | 4.7 (max 14) | – |
| Illumina Ramos VDJ (PRJNA820976) | null | all | 3 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina Ramos VDJ (PRJNA820976) | subclonal | base | 3 | 0.28 / -0.21 | 0.29 (0.00–0.88) | 0.04 (0.00–0.11) | 0.77 (0.57–0.87) | 0.0 (max 0) | 0.26 |
| Illumina Ramos VDJ (PRJNA820976) | subclonal | qmap | 3 | 0.28 / -0.21 | 0.29 (0.00–0.88) | 0.04 (0.00–0.11) | 0.79 (0.67–0.86) | 0.0 (max 0) | 0.26 |
| Illumina Ramos VDJ (PRJNA820976) | subclonal | pcr3 | 3 | 0.32 / -0.26 | 0.46 (0.00–1.37) | 0.18 (0.00–0.55) | 1.42 (0.92–1.91) | 6.3 (max 19) | 0.09 |
| Illumina Ramos VDJ (PRJNA820976) | subclonal | pcr10 | 3 | 0.28 / -0.15 | 0.38 (0.00–1.14) | 0.33 (0.00–1.00) | 1.97 (1.13–3.30) | 4.0 (max 12) | 0.20 |
| Illumina Ramos VDJ (PRJNA820976) | subclonal | deep | 3 | 0.37 / -0.21 | 0.32 (0.00–0.95) | 0.08 (0.00–0.23) | 0.73 (0.59–0.96) | 0.0 (max 0) | 0.17 |
| Illumina Ramos VDJ (PRJNA820976) | subclonal | all | 3 | 0.42 / -0.31 | 0.37 (0.00–1.10) | 0.22 (0.00–0.67) | 1.47 (1.14–2.06) | 4.7 (max 14) | 0.17 |
| ONT Sμ, LPS blasts | null | base | 3 | – | – | 0.03 (est) | – | 0.3 (max 1) | – |
| ONT Sμ, LPS blasts | null | pcr3 | 3 | – | – | 0.36 (est) | – | 4.3 (max 6) | – |
| ONT Sμ, LPS blasts | null | pcr10 | 3 | – | – | 0.91 (est) | – | 11.0 (max 17) | – |
| ONT Sμ, LPS blasts | null | all | 3 | – | – | 0.33 (est) | – | 4.0 (max 7) | – |
| ONT Sμ, LPS blasts | polyclonal | base | 3 | 0.98 / -0.04 | 0.95 (0.91–0.97) | 0.35 (0.30–0.39) | 0.96 (0.75–1.16) | 0.7 (max 1) | 0.84 |
| ONT Sμ, LPS blasts | polyclonal | qmap | 3 | 0.97 / -0.04 | 0.95 (0.91–0.97) | 0.35 (0.30–0.39) | 0.94 (0.75–1.11) | 0.7 (max 1) | 0.83 |
| ONT Sμ, LPS blasts | polyclonal | pcr3 | 3 | 1.21 / +0.09 | 1.07 (1.05–1.10) | 0.92 (0.65–1.28) | 1.58 (1.42–1.85) | 19.3 (max 37) | 0.77 |
| ONT Sμ, LPS blasts | polyclonal | pcr10 | 3 | 1.31 / +0.11 | 1.09 (1.06–1.15) | 1.07 (0.85–1.25) | 1.53 (1.11–2.13) | 13.0 (max 19) | 0.89 |
| ONT Sμ, LPS blasts | polyclonal | deep | 3 | 1.11 / +0.00 | 1.00 (0.99–1.01) | 0.57 (0.54–0.65) | 0.92 (0.77–1.02) | 3.3 (max 7) | 0.81 |
| ONT Sμ, LPS blasts | polyclonal | all | 3 | 1.08 / +0.07 | 1.07 (1.03–1.12) | 0.97 (0.81–1.14) | 1.46 (1.35–1.57) | 22.0 (max 29) | 0.83 |
| ONT Sμ, CH12 | null | base | 3 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| ONT Sμ, CH12 | null | pcr3 | 3 | – | – | 0.72 (est) | – | 8.7 (max 16) | – |
| ONT Sμ, CH12 | null | pcr10 | 3 | – | – | 1.02 (est) | – | 12.3 (max 16) | – |
| ONT Sμ, CH12 | null | all | 3 | – | – | 0.77 (est) | – | 9.3 (max 16) | – |
| ONT Sμ, CH12 | subclonal | base | 3 | 1.26 / -0.08 | 0.94 (0.91–0.97) | 0.41 (0.25–0.62) | 0.88 (0.86–0.92) | 4.3 (max 8) | 0.84 |
| ONT Sμ, CH12 | subclonal | qmap | 3 | 1.26 / -0.06 | 0.95 (0.92–0.98) | 0.44 (0.25–0.65) | 0.93 (0.92–0.94) | 4.7 (max 8) | 0.85 |
| ONT Sμ, CH12 | subclonal | pcr3 | 3 | 1.35 / +0.41 | 1.30 (1.10–1.66) | 0.95 (0.79–1.12) | 1.25 (1.17–1.32) | 36.0 (max 49) | 0.91 |
| ONT Sμ, CH12 | subclonal | pcr10 | 3 | 1.24 / +0.16 | 1.14 (1.12–1.16) | 1.03 (0.89–1.15) | 1.14 (1.05–1.22) | 22.3 (max 33) | 0.93 |
| ONT Sμ, CH12 | subclonal | deep | 3 | 1.77 / +0.35 | 1.19 (1.04–1.45) | 0.75 (0.69–0.82) | 0.92 (0.91–0.92) | 7.3 (max 10) | 0.92 |
| ONT Sμ, CH12 | subclonal | all | 3 | 1.41 / +0.13 | 1.09 (1.06–1.11) | 0.97 (0.87–1.16) | 1.15 (1.11–1.20) | 38.3 (max 55) | 0.91 |
| ONT Pim1 | null | base | 3 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| ONT Pim1 | null | pcr3 | 3 | – | – | 1.42 (est) | – | 24.3 (max 26) | – |
| ONT Pim1 | null | pcr10 | 3 | – | – | 1.44 (est) | – | 24.7 (max 27) | – |
| ONT Pim1 | null | all | 3 | – | – | 1.38 (est) | – | 23.7 (max 25) | – |
| ONT OT1 (CH12) | null | base | 3 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| ONT OT1 (CH12) | null | pcr3 | 3 | – | – | 0.16 (est) | – | 3.0 (max 6) | – |
| ONT OT1 (CH12) | null | pcr10 | 3 | – | – | 0.31 (est) | – | 6.0 (max 15) | – |
| ONT OT1 (CH12) | null | all | 3 | – | – | 0.16 (est) | – | 3.0 (max 6) | – |
| ONT OT2 (CH12) | null | base | 3 | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| ONT OT2 (CH12) | null | pcr3 | 3 | – | – | 0.84 (est) | – | 15.7 (max 31) | – |
| ONT OT2 (CH12) | null | pcr10 | 3 | – | – | 0.66 (est) | – | 12.3 (max 18) | – |
| ONT OT2 (CH12) | null | all | 3 | – | – | 0.88 (est) | – | 16.3 (max 31) | – |
| ONT J꜀4 (this study) | null | base | 3 | – | – | 0.25 (est) | – | 0.7 (max 2) | – |
| ONT J꜀4 (this study) | null | pcr3 | 3 | – | – | 0.38 (est) | – | 1.0 (max 3) | – |
| ONT J꜀4 (this study) | null | pcr10 | 3 | – | – | 0.25 (est) | – | 0.7 (max 2) | – |
| ONT J꜀4 (this study) | null | all | 3 | – | – | 0.38 (est) | – | 1.0 (max 3) | – |
| ONT J꜀4 (this study) | polyclonal | base | 3 | 13.80 / +0.01 | 1.00 (0.99–1.02) | 0.97 (0.94–1.01) | 0.99 (0.96–1.01) | 12.0 (max 15) | 0.98 |
| ONT J꜀4 (this study) | polyclonal | qmap | 3 | 13.67 / +0.01 | 1.00 (0.99–1.02) | 0.96 (0.94–1.00) | 0.98 (0.96–1.00) | 11.0 (max 14) | 0.98 |
| ONT J꜀4 (this study) | polyclonal | pcr3 | 3 | 14.19 / +0.20 | 1.02 (1.01–1.02) | 1.02 (0.99–1.05) | 1.04 (1.02–1.05) | 27.0 (max 33) | 0.98 |
| ONT J꜀4 (this study) | polyclonal | pcr10 | 3 | 13.45 / +0.24 | 1.02 (1.00–1.03) | 1.02 (0.99–1.06) | 1.03 (1.01–1.06) | 18.0 (max 25) | 0.99 |
| ONT J꜀4 (this study) | polyclonal | deep | 3 | 24.76 / +0.55 | 1.02 (1.02–1.03) | 1.00 (0.99–1.00) | 1.00 (0.99–1.00) | 20.7 (max 25) | 0.99 |
| ONT J꜀4 (this study) | polyclonal | all | 3 | 23.47 / +0.36 | 1.02 (1.01–1.02) | 1.01 (1.00–1.02) | 1.01 (1.00–1.02) | 29.3 (max 33) | 0.98 |
| ONT N114 (this study) | null | base | 3 | – | – | 0.11 (est) | – | 0.3 (max 1) | – |
| ONT N114 (this study) | null | pcr3 | 3 | – | – | 3.27 (est) | – | 10.0 (max 21) | – |
| ONT N114 (this study) | null | pcr10 | 3 | – | – | 2.40 (est) | – | 7.3 (max 9) | – |
| ONT N114 (this study) | null | all | 3 | – | – | 3.27 (est) | – | 10.0 (max 21) | – |
| ONT N114 (this study) | subclonal | base | 3 | 1.70 / -0.01 | 0.99 (0.97–1.01) | 0.87 (0.83–0.89) | 0.97 (0.95–0.99) | 6.7 (max 8) | 0.99 |
| ONT N114 (this study) | subclonal | qmap | 3 | 1.69 / -0.01 | 0.99 (0.97–1.02) | 0.87 (0.83–0.90) | 0.98 (0.96–1.00) | 6.7 (max 8) | 0.99 |
| ONT N114 (this study) | subclonal | pcr3 | 3 | 1.57 / +0.15 | 1.10 (1.05–1.13) | 1.12 (1.10–1.14) | 1.22 (1.15–1.29) | 24.7 (max 30) | 0.99 |
| ONT N114 (this study) | subclonal | pcr10 | 3 | 1.48 / +0.11 | 1.08 (1.07–1.08) | 1.12 (1.10–1.14) | 1.20 (1.14–1.25) | 17.7 (max 18) | 0.97 |
| ONT N114 (this study) | subclonal | deep | 3 | 3.64 / +0.12 | 1.03 (1.02–1.04) | 0.99 (0.99–0.99) | 1.00 (0.99–1.01) | 14.3 (max 18) | 0.99 |
| ONT N114 (this study) | subclonal | all | 3 | 2.81 / +0.15 | 1.06 (1.03–1.10) | 1.06 (1.04–1.11) | 1.08 (1.05–1.12) | 29.0 (max 33) | 0.98 |

## The opt-in caller fixes, against their own base

Each fix is off by default; `base` here is the same libraries called without it. A read-end margin also puts the planted mutations in the trimmed columns out of reach, so recall falls by about the share of the window it removes.

| Set | Scenario | Caller | Runs | Truth / error | Mutation frequency | Frequency, counted bases | Unique density | Combined density | False calls | Recall, 2+ reads |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Illumina J꜀4 (PRJNA973112) | null | base | 8 | – | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina J꜀4 (PRJNA973112) | null | margin4 | 8 | – | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina J꜀4 (PRJNA973112) | null | margin8 | 8 | – | – | – | 0.21 (est) | – | 0.6 (max 3) | – |
| Illumina J꜀4 (PRJNA973112) | null | margin8het | 8 | – | – | – | 0.21 (est) | – | 0.6 (max 3) | – |
| Illumina J꜀4 (PRJNA973112) | polyclonal | base | 8 | 3.09 / -0.62 | 0.78 (0.47–0.94) | 0.81 (0.50–0.98) | 0.16 (0.05–0.30) | 0.58 (0.18–0.99) | 0.0 (max 0) | 0.43 |
| Illumina J꜀4 (PRJNA973112) | polyclonal | margin4 | 8 | 3.09 / -0.84 | 0.71 (0.47–0.92) | 0.74 (0.50–0.93) | 0.12 (0.05–0.22) | 0.53 (0.18–0.99) | 0.0 (max 0) | 0.31 |
| Illumina J꜀4 (PRJNA973112) | polyclonal | margin8 | 8 | 3.09 / -1.00 | 0.65 (0.36–0.83) | 0.80 (0.41–1.03) | 0.16 (0.05–0.35) | 0.59 (0.25–1.04) | 0.5 (max 2) | 0.39 |
| Illumina J꜀4 (PRJNA973112) | polyclonal | margin8het | 8 | 3.09 / -1.00 | 0.65 (0.36–0.83) | 0.80 (0.41–1.03) | 0.16 (0.05–0.35) | 0.59 (0.25–1.04) | 0.5 (max 2) | 0.39 |
| Illumina CH12 VDJ (PRJNA918596) | null | base | 8 | – | – | – | 0.04 (est) | – | 0.5 (max 2) | – |
| Illumina CH12 VDJ (PRJNA918596) | null | het | 8 | – | – | – | 0.04 (est) | – | 0.5 (max 2) | – |
| Illumina CH12 VDJ (PRJNA918596) | subclonal | base | 8 | 0.21 / +0.00 | 0.93 (0.47–1.40) | – | 0.16 (0.05–0.41) | 1.27 (0.99–1.70) | 1.1 (max 6) | 0.81 |
| Illumina CH12 VDJ (PRJNA918596) | subclonal | het | 8 | 0.21 / +0.00 | 0.93 (0.47–1.40) | – | 0.16 (0.05–0.41) | 1.27 (0.99–1.70) | 1.1 (max 6) | 0.81 |
| Illumina Ramos VDJ (PRJNA820976) | null | base | 6 | – | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina Ramos VDJ (PRJNA820976) | null | het | 6 | – | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| Illumina Ramos VDJ (PRJNA820976) | subclonal | base | 6 | 0.32 / -0.29 | 0.15 (0.00–0.88) | – | 0.02 (0.00–0.11) | 0.73 (0.57–0.87) | 0.0 (max 0) | 0.13 |
| Illumina Ramos VDJ (PRJNA820976) | subclonal | het | 6 | 0.32 / -0.29 | 0.15 (0.00–0.88) | – | 0.02 (0.00–0.11) | 0.73 (0.57–0.87) | 0.0 (max 0) | 0.13 |
| ONT Sμ, CH12 | null | base | 6 | – | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| ONT Sμ, CH12 | null | het | 6 | – | – | – | 0.00 (est) | – | 0.0 (max 0) | – |
| ONT Sμ, CH12 | subclonal | base | 6 | 1.36 / -0.06 | 0.95 (0.91–1.02) | – | 0.46 (0.25–0.77) | 0.89 (0.86–0.93) | 5.8 (max 15) | 0.86 |
| ONT Sμ, CH12 | subclonal | het | 6 | 1.36 / -0.06 | 0.96 (0.94–1.02) | – | 0.47 (0.30–0.77) | 0.89 (0.84–0.93) | 6.2 (max 15) | 0.86 |
| ONT J꜀4 (this study) | null | base | 6 | – | – | – | 0.25 (est) | – | 0.7 (max 2) | – |
| ONT J꜀4 (this study) | null | margin4 | 6 | – | – | – | 0.25 (est) | – | 0.7 (max 2) | – |
| ONT J꜀4 (this study) | null | margin8 | 6 | – | – | – | 0.25 (est) | – | 0.7 (max 2) | – |
| ONT J꜀4 (this study) | polyclonal | base | 6 | 15.22 / +0.05 | 1.00 (0.99–1.02) | 1.01 (1.00–1.02) | 0.97 (0.94–1.02) | 0.99 (0.96–1.02) | 13.8 (max 22) | 0.98 |
| ONT J꜀4 (this study) | polyclonal | margin4 | 6 | 15.22 / +0.05 | 1.00 (0.99–1.02) | 1.01 (1.00–1.02) | 0.97 (0.94–1.02) | 0.99 (0.96–1.02) | 13.8 (max 22) | 0.98 |
| ONT J꜀4 (this study) | polyclonal | margin8 | 6 | 15.22 / +0.04 | 1.00 (0.99–1.02) | 1.01 (1.00–1.02) | 0.97 (0.94–1.02) | 0.99 (0.96–1.02) | 13.8 (max 22) | 0.98 |
