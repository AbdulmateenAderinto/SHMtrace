#!/usr/bin/env python3
"""Rebuild Supplementary Figures S5-S7 with explicit bounded cell geometry.

This is a ground-up renderer.  Every matrix count, cell percentage, and row
total is registered against its Rectangle patch and verified in display pixels
before the figure is saved.  Export fails if any registered text leaves its
cell boundary.
"""

from __future__ import annotations

import csv
from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle


from config import BASE as _BASE, DATA_SANGER as _DS, OUT_DIR as _OUT
ROOT = Path(_BASE)                       # original inputs (nanopore, hybridoma): read only
SANGER = Path(_DS)                       # corrected Sanger calls
OUT = Path(_OUT)
PDF_OUT = OUT / "pdf"
OUT.mkdir(exist_ok=True)
PDF_OUT.mkdir(parents=True, exist_ok=True)

BASES = ["A", "T", "C", "G"]
INDEX = {base: index for index, base in enumerate(BASES)}
TYPES = [f"{source}>{target}" for source in BASES for target in BASES if source != target]
SEQ = LinearSegmentedColormap.from_list(
    "bounded_seqblue", ["#F5F8FC", "#C4D8EE", "#6FA0CE", "#2F6BA6", "#173D63"]
)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "custom",
        "mathtext.rm": "DejaVu Sans",
        "mathtext.bf": "DejaVu Sans:bold",
        "mathtext.it": "DejaVu Sans:italic:bold",
    }
)


def as_float(value: str | None) -> float | None:
    try:
        return float(str(value).replace("%", "").strip())
    except (TypeError, ValueError):
        return None


def unique_variants(path: Path, window: tuple[int, int] | None) -> set[tuple[int, str, str]]:
    variants: set[tuple[int, str, str]] = set()
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle):
            if "SNP" not in row.get("Polymorphism Type", ""):
                continue
            frequency = as_float(row.get("Variant Frequency"))
            p_value = as_float(row.get("Variant P-Value (approximate)"))
            if frequency is None or frequency >= 50 or p_value is None or p_value > 0.05:
                continue
            try:
                position = int(row["Minimum"])
            except (KeyError, TypeError, ValueError):
                continue
            if window is not None and not (window[0] <= position <= window[1]):
                continue
            source = (row.get("Reference Nucleotide(s)") or "").strip().upper()
            target = (row.get("Variant Nucleotide(s)") or "").strip().upper()
            if source in INDEX and target in INDEX and source != target:
                variants.add((position, source, target))
    return variants


def mutation_matrix(variants: set[tuple[int, str, str]]) -> np.ndarray:
    matrix = np.zeros((4, 4), dtype=int)  # [to, from]
    for _, source, target in variants:
        matrix[INDEX[target], INDEX[source]] += 1
    return matrix


def display_label(label: str) -> str:
    if label.startswith("Aicda−/−"):
        return label.replace("Aicda−/−", r"$\mathbfit{Aicda}^{\mathbf{-/-}}$", 1)
    if label.startswith("pMX-Aicda"):
        return label.replace("pMX-Aicda", r"pMX-$\mathbfit{Aicda}$", 1)
    return label


def cell_percentage(value: float, total: int) -> str:
    if total == 0:
        return "—"
    if abs(value - round(value)) < 1e-9:
        return f"{int(round(value))}%"
    return f"{value:.1f}%"


def row_percentage(value: float, total: int) -> str:
    if total == 0:
        return "—"
    if abs(value - round(value)) < 1e-9:
        return f"{int(round(value))}"
    return f"{value:.1f}"


def bounded_text(
    ax,
    bounds: list[tuple[object, Rectangle, str]],
    rectangle: Rectangle,
    x: float,
    y: float,
    value: str,
    **kwargs,
):
    text = ax.text(x, y, value, clip_on=False, **kwargs)
    bounds.append((text, rectangle, value))
    return text


def draw_matrix(
    ax,
    matrix: np.ndarray,
    title: str,
    vmax: float,
    bounds: list[tuple[object, Rectangle, str]],
) -> None:
    shown = matrix.T  # rows are From; columns are To
    total = int(shown.sum())
    percentages = 100 * shown / total if total else np.zeros((4, 4), dtype=float)
    row_percentages = 100 * shown.sum(axis=1) / total if total else np.zeros(4)
    normalizer = Normalize(0, max(vmax, 1.0))

    cell_width = 1.12
    cell_height = 1.00
    step_x = 1.16
    step_y = 1.04
    matrix_x = 0.00
    matrix_y = 0.00
    summary_x = 5.18
    summary_width = 1.36

    for row in range(4):
        y = matrix_y + (3 - row) * step_y
        for column in range(4):
            x = matrix_x + column * step_x
            if row == column:
                patch = Rectangle((x, y), cell_width, cell_height, facecolor="white", edgecolor="white")
                ax.add_patch(patch)
                bounded_text(
                    ax,
                    bounds,
                    patch,
                    x + cell_width / 2,
                    y + cell_height / 2,
                    "—",
                    ha="center",
                    va="center",
                    fontsize=8.0,
                    color="#A9A9A9",
                )
                continue

            value = float(percentages[row, column])
            patch = Rectangle(
                (x, y),
                cell_width,
                cell_height,
                facecolor=SEQ(normalizer(value)),
                edgecolor="white",
                linewidth=0.9,
            )
            ax.add_patch(patch)
            colour = "white" if normalizer(value) > 0.56 else "#151515"
            bounded_text(
                ax,
                bounds,
                patch,
                x + cell_width / 2,
                y + 0.64,
                f"{int(shown[row, column])}",
                ha="center",
                va="center",
                fontsize=8.0,
                fontweight="bold",
                color=colour,
            )
            bounded_text(
                ax,
                bounds,
                patch,
                x + cell_width / 2,
                y + 0.27,
                cell_percentage(value, total),
                ha="center",
                va="center",
                fontsize=4.25,
                color=colour,
            )

        row_value = float(row_percentages[row])
        summary_patch = Rectangle(
            (summary_x, y),
            summary_width,
            cell_height,
            facecolor=SEQ(0.60 * row_value / max(float(row_percentages.max()), 1.0)),
            edgecolor="white",
            linewidth=0.9,
        )
        ax.add_patch(summary_patch)
        bounded_text(
            ax,
            bounds,
            summary_patch,
            summary_x + summary_width / 2,
            y + cell_height / 2,
            row_percentage(row_value, total),
            ha="center",
            va="center",
            fontsize=7.6,
            fontweight="bold",
            color="#101010",
        )

    column_centres = [matrix_x + column * step_x + cell_width / 2 for column in range(4)]
    row_centres = [matrix_y + (3 - row) * step_y + cell_height / 2 for row in range(4)]
    for centre, base in zip(column_centres, BASES):
        ax.text(centre, 4.45, base, ha="center", va="bottom", fontsize=8.6, fontweight="bold")
    ax.text(summary_x + summary_width / 2, 4.45, "%", ha="center", va="bottom", fontsize=8.6, fontweight="bold")
    ax.text(np.mean(column_centres), 4.95, "To", ha="center", va="bottom", fontsize=9.1, fontweight="bold")
    for centre, base in zip(row_centres, BASES):
        ax.text(-0.46, centre, base, ha="center", va="center", fontsize=8.6, fontweight="bold")
    ax.text(-1.05, 2.06, "From", ha="center", va="center", fontsize=8.6, fontweight="bold", rotation=90)
    ax.text(2.70, 5.58, title, ha="center", va="bottom", fontsize=9.2, fontweight="bold")

    ax.set_xlim(-1.38, 6.72)
    ax.set_ylim(-0.20, 5.82)
    ax.set_aspect("equal")
    ax.axis("off")


def verify_bounds(fig, bounds: list[tuple[object, Rectangle, str]], figure_name: str) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    failures = []
    for text, rectangle, value in bounds:
        text_box = text.get_window_extent(renderer=renderer)
        cell_box = rectangle.get_window_extent(renderer=renderer)
        # Require a small internal margin so touching the border also fails.
        margin = 0.35
        if not (
            text_box.x0 >= cell_box.x0 + margin
            and text_box.x1 <= cell_box.x1 - margin
            and text_box.y0 >= cell_box.y0 + margin
            and text_box.y1 <= cell_box.y1 - margin
        ):
            failures.append((value, tuple(round(v, 2) for v in text_box.bounds), tuple(round(v, 2) for v in cell_box.bounds)))
    if failures:
        raise RuntimeError(f"{figure_name}: {len(failures)} text-boundary failures: {failures[:8]}")


def matrix_vmax(matrices: list[np.ndarray]) -> float:
    values = []
    for matrix in matrices:
        total = int(matrix.sum())
        if total:
            displayed = 100 * matrix.T / total
            values.extend(displayed[~np.eye(4, dtype=bool)].tolist())
    return max(values, default=1.0)


def make_figure(
    filename: str,
    groups: list[tuple[str, OrderedDict[str, Path], tuple[int, int] | None]],
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    row_data: list[tuple[str, list[tuple[str, np.ndarray]]]] = []
    for group_name, samples, window in groups:
        matrices = []
        for sample_name, path in samples.items():
            matrix = mutation_matrix(unique_variants(path, window))
            matrices.append((sample_name, matrix))
            record: dict[str, object] = {
                "figure": filename,
                "group": group_name,
                "sample": sample_name,
                "source_file": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else path.name,
                "unique_variants": int(matrix.sum()),
            }
            for mutation in TYPES:
                source, target = mutation.split(">")
                record[mutation] = int(matrix[INDEX[target], INDEX[source]])
            records.append(record)
        row_data.append((group_name, matrices))

    # Reviewer request: tighten white space. Wider/shorter panels so the matrices
    # fill each cell, and the two groups sit closer together.
    fig = plt.figure(figsize=(8.4, 5.5))
    panel_positions = [
        (0.045, 0.590, 0.300, 0.300),
        (0.360, 0.590, 0.300, 0.300),
        (0.675, 0.590, 0.300, 0.300),
        (0.045, 0.130, 0.300, 0.300),
        (0.360, 0.130, 0.300, 0.300),
        (0.675, 0.130, 0.300, 0.300),
    ]
    all_bounds: list[tuple[object, Rectangle, str]] = []
    for row_index, (group_name, matrices) in enumerate(row_data):
        vmax = matrix_vmax([matrix for _, matrix in matrices])
        for column_index, (sample_name, matrix) in enumerate(matrices):
            ax = fig.add_axes(panel_positions[row_index * 3 + column_index])
            draw_matrix(
                ax,
                matrix,
                f"{display_label(sample_name)} ($n$ = {int(matrix.sum())})",
                vmax,
                all_bounds,
            )

        label_y = 0.945 if row_index == 0 else 0.485
        fig.text(0.012, label_y, chr(ord("A") + row_index), fontsize=17, fontweight="bold", va="top")
        fig.text(0.058, label_y, display_label(group_name), fontsize=16, fontweight="bold", va="top")
        cbar_y = 0.545 if row_index == 0 else 0.075
        colour_axis = fig.add_axes([0.355, cbar_y, 0.29, 0.012])
        colourbar = fig.colorbar(
            plt.cm.ScalarMappable(cmap=SEQ, norm=Normalize(0, vmax)),
            cax=colour_axis,
            orientation="horizontal",
        )
        colourbar.set_label("% of sample's unique variants", fontsize=8.0, labelpad=2)
        colourbar.ax.tick_params(labelsize=7.0, length=2)

    verify_bounds(fig, all_bounds, filename)
    for extension in ("png", "pdf"):
        path = OUT / f"{filename}.{extension}"
        fig.savefig(path, dpi=600 if extension == "png" else None, bbox_inches="tight", pad_inches=0.035)
        if extension == "pdf":
            (PDF_OUT / path.name).write_bytes(path.read_bytes())
    plt.close(fig)
    return records


def main() -> None:
    sanger_wt = OrderedDict(
        [
            ("Wt 2", SANGER / "WT2 corrected Sanger Annotations.csv"),
            ("Wt 3", SANGER / "WT3 corrected Sanger Annotations.csv"),
            ("Wt 4", SANGER / "WT4 corrected Sanger Annotations.csv"),
        ]
    )
    sanger_ko = OrderedDict(
        [
            ("Aicda−/− 1", SANGER / "Aicda_A1 corrected Sanger Annotations.csv"),
            ("Aicda−/− 2", SANGER / "Aicda_A2 corrected Sanger Annotations.csv"),
            ("Aicda−/− 3", SANGER / "Aicda_A3 corrected Sanger Annotations.csv"),
        ]
    )

    # JH4 ONT: SHMtrace calls (docs/calibration/jh4_shmtrace.py), 500 reads per mouse;
    # previously the Geneious sub500 exports in "July 13th JH4"
    jh4 = Path(__file__).resolve().parent.parent / "jh4_calls"
    ont_wt = OrderedDict([(f"Wt {i}", jh4 / f"WT{i} SHMtrace Annotations.csv") for i in (2, 3, 4)])
    ont_ko = OrderedDict([(f"Aicda−/− {i}", jh4 / f"KO{i} SHMtrace Annotations.csv") for i in (1, 2, 3)])

    # N114: the SHMtrace re-analysis of the deposited reads (docs/calibration/n114_shmtrace.py);
    # previously the Geneious sub500 exports in "July 10th /Hybridoma"
    n114 = Path(__file__).resolve().parent.parent / "n114_calls"
    aid = OrderedDict([(k, n114 / f"{k} SHMtrace Annotations.csv") for k in ("C1", "D3", "D5")])
    empty = OrderedDict([(k, n114 / f"{k} SHMtrace Annotations.csv") for k in ("C3", "C4", "D1")])

    records = []
    records.extend(
        make_figure(
            "FigureS5_individual_Sanger_matrices",
            [("Wt", sanger_wt, None), ("Aicda−/−", sanger_ko, None)],
        )
    )
    records.extend(
        make_figure(
            "FigureS6_individual_JH4_ONT_matrices",
            [("Wt", ont_wt, (101, 629)), ("Aicda−/−", ont_ko, (101, 629))],
        )
    )
    records.extend(
        make_figure(
            "FigureS7_individual_N114_matrices",
            [("pMX-Aicda", aid, (21, 635)), ("pMX-Empty vector", empty, (21, 635))],
        )
    )

    csv_path = OUT / "Individual_mutation_matrices_counts.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = ["figure", "group", "sample", "source_file", "unique_variants", *TYPES]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    for record in records:
        print(f"{record['figure']}: {record['group']} {record['sample']} n={record['unique_variants']}")
    print("All registered matrix text passed pixel-boundary validation.")


if __name__ == "__main__":
    main()
