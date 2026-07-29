#!/usr/bin/env python3
"""
H3K9me2 20-kb compartment analysis and Fig. 5a track preparation.

Analyses retained from the manuscript workflow:
  1. Per-bin H3K9me2 log2FC (spermatocyte / wing)
  2. Eu/PCH/Telo compartment summaries
  3. Kruskal-Wallis and pairwise Mann-Whitney U tests
  4. Late-primary-spermatocyte TE signal extraction
  5. Smoothed chromosome-track values

Required inputs:
  H3K9me2_sperm_wing_20kb_counts.tab
  SoloTE_late_primary_spermatocyte.bw

Required packages:
  numpy, pandas, pyBigWig, scipy
"""

from __future__ import annotations

import importlib.metadata
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import kruskal, mannwhitneyu


# ------------------------------------------------------------------------------
# Inputs, outputs, and parameters
# ------------------------------------------------------------------------------

COUNTS_TABLE = Path("H3K9me2_sperm_wing_20kb_counts.tab")
TE_BIGWIG = Path("SoloTE_late_primary_spermatocyte.bw")
OUTPUT_DIR = Path("H3K9me2_20kb_results")

PSEUDOCOUNT = 1e-3
INPUT_BIN_SIZE_BP = 20_000
DISPLAY_SMOOTHING_BINS = 5       # 100 kb
DISPLAY_LOG2FC_CLIP = 3.0
TE_PERCENTILE_FLOOR = 0.70

TARGET_CHROMOSOMES = ("2L", "2R", "3L", "3R", "X")
COMPARTMENT_ORDER = ("Eu", "PCH", "Telo")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------------------------
# Chromatin-compartment boundaries
# Coordinates are 1-based inclusive.
# ------------------------------------------------------------------------------

BOUNDARIES = pd.DataFrame(
    [
        ("2L", 1, 34_496, "Telo"),
        ("2L", 34_497, 21_806_874, "Eu"),
        ("2L", 21_806_875, 23_219_612, "PCH"),
        ("2R", 1, 4_651_836, "PCH"),
        ("2R", 4_651_837, 24_535_303, "Eu"),
        ("2R", 24_535_304, 24_550_778, "Telo"),
        ("3L", 1, 29_781, "Telo"),
        ("3L", 29_782, 22_811_581, "Eu"),
        ("3L", 22_811_582, 27_189_634, "PCH"),
        ("3R", 1, 3_415_214, "PCH"),
        ("3R", 3_415_215, 31_253_423, "Eu"),
        ("3R", 31_253_424, 31_281_498, "Telo"),
        ("X", 1, 138_429, "Telo"),
        ("X", 138_430, 22_110_469, "Eu"),
        ("X", 22_110_470, 22_604_182, "PCH"),
    ],
    columns=["chromosome", "start", "end", "compartment"],
)


# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required input file not found: {path}")


def bh_adjust(p_values: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg correction while preserving missing values."""
    p = np.asarray(list(p_values), dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)

    if not valid.any():
        return out

    valid_p = p[valid]
    order = np.argsort(valid_p)
    ranked = valid_p[order]
    n = len(ranked)

    adjusted_ranked = ranked * n / np.arange(1, n + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0, 1)

    adjusted_valid = np.empty(n, dtype=float)
    adjusted_valid[order] = adjusted_ranked
    out[valid] = adjusted_valid
    return out


def load_counts_table(path: Path) -> pd.DataFrame:
    """Read the first five columns of a multiBigwigSummary table."""
    raw = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        low_memory=False,
    )

    if raw.shape[1] < 5:
        raise ValueError(f"{path} has fewer than five columns.")

    counts = raw.iloc[:, :5].copy()
    counts.columns = [
        "chromosome",
        "start",
        "end",
        "spermatocyte",
        "wing",
    ]

    for column in ("start", "end", "spermatocyte", "wing"):
        counts[column] = pd.to_numeric(counts[column], errors="coerce")

    counts = counts.dropna(
        subset=["start", "end", "spermatocyte", "wing"]
    ).copy()

    counts["chromosome"] = (
        counts["chromosome"]
        .astype(str)
        .str.replace(r"^chr", "", regex=True)
    )
    counts["start"] = counts["start"].astype(int)
    counts["end"] = counts["end"].astype(int)
    counts["spermatocyte"] = counts["spermatocyte"].astype(float)
    counts["wing"] = counts["wing"].astype(float)

    if (counts["end"] <= counts["start"]).any():
        raise ValueError("At least one input interval has end <= start.")

    if (counts[["spermatocyte", "wing"]] < 0).any().any():
        raise ValueError("H3K9me2 signal contains negative values.")

    return counts.reset_index(drop=True)


def assign_compartments(counts: pd.DataFrame) -> pd.Series:
    """Assign Eu/PCH/Telo using each 20-kb bin midpoint."""
    labels = pd.Series(index=counts.index, dtype="object")

    for boundary in BOUNDARIES.itertuples(index=False):
        mask = (
            (counts["chromosome"] == boundary.chromosome)
            & (counts["bin_midpoint"] >= boundary.start)
            & (counts["bin_midpoint"] <= boundary.end)
        )

        if labels.loc[mask].notna().any():
            raise ValueError("A bin matched more than one compartment interval.")

        labels.loc[mask] = boundary.compartment

    return labels


def extract_bigwig_means(
    bigwig_path: Path,
    intervals: pd.DataFrame,
) -> np.ndarray:
    """Extract mean bigWig signal for each interval."""
    values = np.full(len(intervals), np.nan, dtype=float)

    with pyBigWig.open(str(bigwig_path)) as bigwig:
        chromosome_lengths = bigwig.chroms()

        for i, row in enumerate(intervals.itertuples(index=False)):
            if row.chromosome not in chromosome_lengths:
                continue

            start = max(0, int(row.start))
            end = min(int(row.end), int(chromosome_lengths[row.chromosome]))

            if end <= start:
                continue

            value = bigwig.stats(
                row.chromosome,
                start,
                end,
                type="mean",
            )[0]

            values[i] = 0.0 if value is None else float(value)

    return values


def rolling_mean(values: pd.Series, window: int) -> np.ndarray:
    return (
        values
        .rolling(window=window, center=True, min_periods=1)
        .mean()
        .to_numpy()
    )


def prepare_figure5a_tracks(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Smooth sperm and wing separately before calculating displayed H3K9me2
    log2FC. Smooth TE hotspot scores over the same 100-kb window.
    """
    output = []

    for chromosome in TARGET_CHROMOSOMES:
        chromosome_data = (
            counts[counts["chromosome"] == chromosome]
            .sort_values("bin_midpoint")
            .copy()
        )

        if chromosome_data.empty:
            continue

        sperm_smooth = rolling_mean(
            chromosome_data["spermatocyte"],
            DISPLAY_SMOOTHING_BINS,
        )
        wing_smooth = rolling_mean(
            chromosome_data["wing"],
            DISPLAY_SMOOTHING_BINS,
        )

        displayed_log2fc = np.log2(
            (sperm_smooth + PSEUDOCOUNT)
            / (wing_smooth + PSEUDOCOUNT)
        )

        chromosome_data["display_H3K9me2_log2FC"] = np.clip(
            displayed_log2fc,
            -DISPLAY_LOG2FC_CLIP,
            DISPLAY_LOG2FC_CLIP,
        )

        chromosome_data["display_TE_hotspot_score"] = np.clip(
            rolling_mean(
                chromosome_data["TE_hotspot_score"],
                DISPLAY_SMOOTHING_BINS,
            ),
            0,
            1,
        )

        output.append(chromosome_data)

    if not output:
        raise ValueError("No chromosome-track data were generated.")

    return pd.concat(output, ignore_index=True)


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


# ------------------------------------------------------------------------------
# Load H3K9me2 data and calculate log2FC
# ------------------------------------------------------------------------------

require_file(COUNTS_TABLE)
require_file(TE_BIGWIG)

counts = load_counts_table(COUNTS_TABLE)
counts = counts[
    counts["chromosome"].isin(TARGET_CHROMOSOMES)
].copy()

if counts.empty:
    raise ValueError("No rows remained for the five target chromosome arms.")

counts["bin_size_bp"] = counts["end"] - counts["start"]
counts["bin_midpoint"] = (counts["start"] + counts["end"]) // 2
counts["compartment"] = assign_compartments(counts)

counts[counts["compartment"].isna()].to_csv(
    OUTPUT_DIR / "H3K9me2_20kb_unassigned_bins.csv",
    index=False,
)

counts = counts[
    counts["compartment"].isin(COMPARTMENT_ORDER)
].copy()

if counts.empty:
    raise ValueError("No bins were assigned to Eu, PCH, or Telo.")

counts["compartment"] = pd.Categorical(
    counts["compartment"],
    categories=COMPARTMENT_ORDER,
    ordered=True,
)

counts["H3K9me2_log2FC"] = np.log2(
    (counts["spermatocyte"] + PSEUDOCOUNT)
    / (counts["wing"] + PSEUDOCOUNT)
)
counts["H3K9me2_log2FC_clipped"] = counts[
    "H3K9me2_log2FC"
].clip(-DISPLAY_LOG2FC_CLIP, DISPLAY_LOG2FC_CLIP)


# ------------------------------------------------------------------------------
# Compartment summaries
# ------------------------------------------------------------------------------

compartment_summary = (
    counts
    .groupby("compartment", observed=False)
    .agg(
        n_bins=("H3K9me2_log2FC", "size"),
        mean_spermatocyte=("spermatocyte", "mean"),
        mean_wing=("wing", "mean"),
        mean_per_bin_log2FC=("H3K9me2_log2FC", "mean"),
        median_per_bin_log2FC=("H3K9me2_log2FC", "median"),
    )
    .reset_index()
)

compartment_summary["compartment_level_log2FC"] = np.log2(
    (compartment_summary["mean_spermatocyte"] + PSEUDOCOUNT)
    / (compartment_summary["mean_wing"] + PSEUDOCOUNT)
)

compartment_summary.to_csv(
    OUTPUT_DIR / "H3K9me2_20kb_compartment_summary.csv",
    index=False,
)

chromosome_summary = (
    counts
    .groupby(["chromosome", "compartment"], observed=False)
    .agg(
        n_bins=("H3K9me2_log2FC", "size"),
        mean_spermatocyte=("spermatocyte", "mean"),
        mean_wing=("wing", "mean"),
        mean_per_bin_log2FC=("H3K9me2_log2FC", "mean"),
        median_per_bin_log2FC=("H3K9me2_log2FC", "median"),
    )
    .reset_index()
)

chromosome_summary = chromosome_summary[
    chromosome_summary["n_bins"] > 0
].copy()

chromosome_summary["compartment_level_log2FC"] = np.log2(
    (chromosome_summary["mean_spermatocyte"] + PSEUDOCOUNT)
    / (chromosome_summary["mean_wing"] + PSEUDOCOUNT)
)

chromosome_summary.to_csv(
    OUTPUT_DIR / "H3K9me2_20kb_chromosome_compartment_summary.csv",
    index=False,
)


# ------------------------------------------------------------------------------
# Kruskal-Wallis and pairwise Mann-Whitney U tests
# ------------------------------------------------------------------------------

values_by_compartment = {
    compartment: counts.loc[
        counts["compartment"] == compartment,
        "H3K9me2_log2FC",
    ].dropna().to_numpy()
    for compartment in COMPARTMENT_ORDER
}

for compartment, values in values_by_compartment.items():
    if len(values) == 0:
        raise ValueError(f"No values were available for {compartment}.")

kw = kruskal(
    *[values_by_compartment[x] for x in COMPARTMENT_ORDER],
    nan_policy="omit",
)

kruskal_table = pd.DataFrame(
    [{
        "comparison": "Eu vs PCH vs Telo",
        "n_Eu": len(values_by_compartment["Eu"]),
        "n_PCH": len(values_by_compartment["PCH"]),
        "n_Telo": len(values_by_compartment["Telo"]),
        "H_statistic": float(kw.statistic),
        "p_value": float(kw.pvalue),
    }]
)

kruskal_table.to_csv(
    OUTPUT_DIR / "H3K9me2_20kb_Kruskal_Wallis.csv",
    index=False,
)

pairwise_rows = []

for group_1, group_2 in (
    ("Eu", "PCH"),
    ("Eu", "Telo"),
    ("PCH", "Telo"),
):
    values_1 = values_by_compartment[group_1]
    values_2 = values_by_compartment[group_2]

    result = mannwhitneyu(
        values_1,
        values_2,
        alternative="two-sided",
        method="auto",
    )

    pairwise_rows.append(
        {
            "comparison": f"{group_1} vs {group_2}",
            "group_1": group_1,
            "group_2": group_2,
            "n_group_1": len(values_1),
            "n_group_2": len(values_2),
            "mean_group_1": np.mean(values_1),
            "mean_group_2": np.mean(values_2),
            "median_group_1": np.median(values_1),
            "median_group_2": np.median(values_2),
            "median_difference_group_1_minus_group_2": (
                np.median(values_1) - np.median(values_2)
            ),
            "Mann_Whitney_U": float(result.statistic),
            "p_value": float(result.pvalue),
        }
    )

pairwise_table = pd.DataFrame(pairwise_rows)
pairwise_table["p_adjusted_BH"] = bh_adjust(
    pairwise_table["p_value"]
)
pairwise_table["significant_FDR_0.05"] = (
    pairwise_table["p_adjusted_BH"] < 0.05
)

pairwise_table.to_csv(
    OUTPUT_DIR / "H3K9me2_20kb_pairwise_Mann_Whitney.csv",
    index=False,
)


# ------------------------------------------------------------------------------
# TE signal and Fig. 5a chromosome-track values
# ------------------------------------------------------------------------------

counts["late_primary_spermatocyte_TE_signal"] = extract_bigwig_means(
    TE_BIGWIG,
    counts,
)
counts["late_primary_spermatocyte_TE_log1p"] = np.log1p(
    counts["late_primary_spermatocyte_TE_signal"].fillna(0)
)
counts["TE_signal_percentile"] = (
    counts["late_primary_spermatocyte_TE_log1p"]
    .rank(method="average", pct=True)
)
counts["TE_hotspot_score"] = (
    (counts["TE_signal_percentile"] - TE_PERCENTILE_FLOOR)
    / (1 - TE_PERCENTILE_FLOOR)
).clip(0, 1)

counts.to_csv(
    OUTPUT_DIR / "H3K9me2_20kb_signal_with_TE.csv",
    index=False,
)

figure5a_tracks = prepare_figure5a_tracks(counts)

figure5a_tracks[
    [
        "chromosome",
        "start",
        "end",
        "bin_midpoint",
        "compartment",
        "spermatocyte",
        "wing",
        "H3K9me2_log2FC",
        "late_primary_spermatocyte_TE_signal",
        "late_primary_spermatocyte_TE_log1p",
        "TE_signal_percentile",
        "TE_hotspot_score",
        "display_H3K9me2_log2FC",
        "display_TE_hotspot_score",
    ]
].to_csv(
    OUTPUT_DIR / "Figure5a_H3K9me2_TE_tracks_20kb.csv",
    index=False,
)


# ------------------------------------------------------------------------------
# Parameters and software versions
# ------------------------------------------------------------------------------

pd.DataFrame(
    {
        "parameter": [
            "H3K9me2 counts table",
            "TE bigWig",
            "Input bin size",
            "Pseudocount",
            "Target chromosomes",
            "Compartments",
            "Omnibus test",
            "Pairwise test",
            "Pairwise correction",
            "Display smoothing",
            "Displayed log2FC clipping",
            "TE percentile floor",
        ],
        "value": [
            str(COUNTS_TABLE),
            str(TE_BIGWIG),
            INPUT_BIN_SIZE_BP,
            PSEUDOCOUNT,
            ", ".join(TARGET_CHROMOSOMES),
            ", ".join(COMPARTMENT_ORDER),
            "Kruskal-Wallis",
            "Two-sided Mann-Whitney U",
            "Benjamini-Hochberg across three comparisons",
            f"{DISPLAY_SMOOTHING_BINS} bins",
            f"-{DISPLAY_LOG2FC_CLIP} to +{DISPLAY_LOG2FC_CLIP}",
            TE_PERCENTILE_FLOOR,
        ],
    }
).to_csv(
    OUTPUT_DIR / "analysis_parameters.csv",
    index=False,
)

pd.DataFrame(
    {
        "software": ["Python", "numpy", "pandas", "pyBigWig", "scipy"],
        "version": [
            sys.version.replace("\n", " "),
            package_version("numpy"),
            package_version("pandas"),
            package_version("pyBigWig"),
            package_version("scipy"),
        ],
    }
).to_csv(
    OUTPUT_DIR / "software_versions.csv",
    index=False,
)


# ------------------------------------------------------------------------------
# Final summary
# ------------------------------------------------------------------------------

print("\nH3K9me2 20-kb analysis complete.")
print(f"Results were written to: {OUTPUT_DIR}")

print("\nCompartment summary:")
print(compartment_summary.to_string(index=False))

print("\nKruskal-Wallis test:")
print(kruskal_table.to_string(index=False))

print("\nPairwise Mann-Whitney tests:")
print(pairwise_table.to_string(index=False))
