#!/usr/bin/env python3

"""
H3K9me2 change and TE-expression correlation in 1-kb genomic bins.

This script:

1. Bins spermatocyte H3K9me2, wing H3K9me2, and late-primary-
   spermatocyte TE-expression bigWigs at approximately 1-kb resolution.
2. Assigns bins to euchromatin (Eu), pericentromeric heterochromatin
   (PCH), or telomeric heterochromatin (Telo).
3. Calculates:
       delta H3K9me2 = spermatocyte - wing
       H3K9me2 log2FC = log2((spermatocyte + pseudocount) /
                             (wing + pseudocount))
4. Tests the genome-wide and compartment-specific Spearman
   relationships between delta H3K9me2 and TE expression.
5. Divides Eu and PCH bins into H3K9me2-delta ventiles and summarizes
   TE expression within each ventile.

Required input files:
    H3K9me2_sperm.RPGC.bw
    H3K9me2_wing.RPGC.bw
    SoloTE_late_primary_spermatocyte.bw

Required Python packages:
    numpy
    pandas
    pyBigWig
    scipy
"""

from __future__ import annotations

import importlib.metadata
import math
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import spearmanr


# ==============================================================================
# 1. Inputs, outputs, and analysis parameters
# ==============================================================================

H3K9_SPERMATOCYTE_BIGWIG = Path("H3K9me2_sperm.RPGC.bw")
H3K9_WING_BIGWIG = Path("H3K9me2_wing.RPGC.bw")
TE_EXPRESSION_BIGWIG = Path("SoloTE_late_primary_spermatocyte.bw")

OUTPUT_DIR = Path("H3K9me2_1kb_TE_correlation_results")

BIN_SIZE_BP = 1_000
H3K9_PSEUDOCOUNT = 1e-3
N_VENTILES = 20

TARGET_CHROMOSOMES = ("2L", "2R", "3L", "3R", "X")
COMPARTMENT_ORDER = ("Eu", "PCH", "Telo")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ==============================================================================
# 2. Genomic-compartment boundaries
#
# Coordinates are 1-based and inclusive.
# ==============================================================================

COMPARTMENT_BOUNDARIES = pd.DataFrame(
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
    columns=("chromosome", "start", "end", "compartment"),
)


# ==============================================================================
# 3. Helper functions
# ==============================================================================

def require_file(path: Path) -> None:
    """Stop when a required input file is absent."""
    if not path.is_file():
        raise FileNotFoundError(f"Required input file not found: {path}")


def canonical_chromosome(name: str) -> str:
    """Remove an optional chr prefix."""
    name = str(name)
    return name[3:] if name.startswith("chr") else name


def chromosome_map(bigwig: pyBigWig.pyBigWig) -> dict[str, str]:
    """
    Map canonical chromosome names to the names used by a bigWig.

    Duplicate canonical names are rejected to avoid ambiguous lookups.
    """
    mapping: dict[str, str] = {}

    for bigwig_name in bigwig.chroms():
        canonical_name = canonical_chromosome(bigwig_name)

        if canonical_name in mapping:
            raise ValueError(
                "The bigWig contains multiple chromosome names that reduce to "
                f"the same canonical name: {canonical_name}"
            )

        mapping[canonical_name] = bigwig_name

    return mapping


def normalize_stats(values: Iterable[float | None]) -> np.ndarray:
    """Convert pyBigWig None values to zero."""
    return np.asarray(
        [
            0.0 if value is None else float(value)
            for value in values
        ],
        dtype=float,
    )


def assign_compartments(
    chromosome: str,
    midpoints: np.ndarray,
) -> np.ndarray:
    """Assign 1-based bin midpoints to Eu, PCH, or Telo."""
    assignments = np.full(
        len(midpoints),
        None,
        dtype=object,
    )

    chromosome_boundaries = COMPARTMENT_BOUNDARIES[
        COMPARTMENT_BOUNDARIES["chromosome"] == chromosome
    ]

    for boundary in chromosome_boundaries.itertuples(index=False):
        in_interval = (
            (midpoints >= boundary.start)
            & (midpoints <= boundary.end)
        )

        if pd.notna(assignments[in_interval]).any():
            raise ValueError(
                f"Overlapping compartment definitions detected for {chromosome}."
            )

        assignments[in_interval] = boundary.compartment

    return assignments


def safe_spearman(
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
) -> tuple[int, float, float]:
    """Return n, Spearman rho, and P value after finite-value filtering."""
    valid = data[
        [x_column, y_column]
    ].replace(
        [np.inf, -np.inf],
        np.nan,
    ).dropna()

    if len(valid) < 3:
        return len(valid), np.nan, np.nan

    x = valid[x_column].to_numpy(dtype=float)
    y = valid[y_column].to_numpy(dtype=float)

    if np.std(x) == 0 or np.std(y) == 0:
        return len(valid), np.nan, np.nan

    result = spearmanr(
        x,
        y,
        alternative="two-sided",
    )

    return (
        len(valid),
        float(result.statistic),
        float(result.pvalue),
    )


def package_version(package_name: str) -> str:
    """Return an installed package version."""
    try:
        return importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


# ==============================================================================
# 4. Open the three bigWigs and identify common chromosome arms
# ==============================================================================

for input_file in (
    H3K9_SPERMATOCYTE_BIGWIG,
    H3K9_WING_BIGWIG,
    TE_EXPRESSION_BIGWIG,
):
    require_file(input_file)

spermatocyte_bigwig = pyBigWig.open(
    str(H3K9_SPERMATOCYTE_BIGWIG)
)
wing_bigwig = pyBigWig.open(
    str(H3K9_WING_BIGWIG)
)
TE_bigwig = pyBigWig.open(
    str(TE_EXPRESSION_BIGWIG)
)

try:
    spermatocyte_map = chromosome_map(
        spermatocyte_bigwig
    )
    wing_map = chromosome_map(
        wing_bigwig
    )
    TE_map = chromosome_map(
        TE_bigwig
    )

    available_chromosomes = [
        chromosome
        for chromosome in TARGET_CHROMOSOMES
        if (
            chromosome in spermatocyte_map
            and chromosome in wing_map
            and chromosome in TE_map
        )
    ]

    if not available_chromosomes:
        raise ValueError(
            "None of the five target chromosome arms is shared by all "
            "three bigWigs."
        )

    missing_chromosomes = [
        chromosome
        for chromosome in TARGET_CHROMOSOMES
        if chromosome not in available_chromosomes
    ]

    if missing_chromosomes:
        print(
            "Warning: the following target chromosomes are absent from at "
            "least one bigWig and will be skipped: "
            + ", ".join(missing_chromosomes),
            file=sys.stderr,
        )

    # ==========================================================================
    # 5. Build the 1-kb genome-wide table
    #
    # pyBigWig's nBins implementation reproduces the original workflow:
    # each full chromosome interval is divided into ceil(length / 1 kb)
    # approximately equal-width bins.
    # ==========================================================================

    chromosome_tables: list[pd.DataFrame] = []

    for chromosome in available_chromosomes:
        spermatocyte_name = spermatocyte_map[
            chromosome
        ]
        wing_name = wing_map[
            chromosome
        ]
        TE_name = TE_map[
            chromosome
        ]

        chromosome_length = min(
            int(
                spermatocyte_bigwig.chroms()[
                    spermatocyte_name
                ]
            ),
            int(
                wing_bigwig.chroms()[
                    wing_name
                ]
            ),
            int(
                TE_bigwig.chroms()[
                    TE_name
                ]
            ),
        )

        n_bins = int(
            math.ceil(
                chromosome_length
                / BIN_SIZE_BP
            )
        )

        spermatocyte_values = normalize_stats(
            spermatocyte_bigwig.stats(
                spermatocyte_name,
                0,
                chromosome_length,
                type="mean",
                nBins=n_bins,
            )
        )

        wing_values = normalize_stats(
            wing_bigwig.stats(
                wing_name,
                0,
                chromosome_length,
                type="mean",
                nBins=n_bins,
            )
        )

        TE_values = normalize_stats(
            TE_bigwig.stats(
                TE_name,
                0,
                chromosome_length,
                type="mean",
                nBins=n_bins,
            )
        )

        n_observed = min(
            len(spermatocyte_values),
            len(wing_values),
            len(TE_values),
        )

        spermatocyte_values = spermatocyte_values[
            :n_observed
        ]
        wing_values = wing_values[
            :n_observed
        ]
        TE_values = TE_values[
            :n_observed
        ]

        bin_start = (
            np.arange(
                n_observed,
                dtype=int,
            )
            * BIN_SIZE_BP
            + 1
        )

        bin_end = np.minimum(
            bin_start
            + BIN_SIZE_BP
            - 1,
            chromosome_length,
        )

        bin_midpoint = (
            bin_start
            + bin_end
        ) // 2

        compartment = assign_compartments(
            chromosome=chromosome,
            midpoints=bin_midpoint,
        )

        chromosome_tables.append(
            pd.DataFrame(
                {
                    "chr": chromosome,
                    "bin_start": bin_start,
                    "bin_end": bin_end,
                    "bin_mid": bin_midpoint,
                    "h3k9_sperm": spermatocyte_values,
                    "h3k9_wing": wing_values,
                    "te_sperm_cpm": TE_values,
                    "region": compartment,
                }
            )
        )

finally:
    spermatocyte_bigwig.close()
    wing_bigwig.close()
    TE_bigwig.close()


binned_data = pd.concat(
    chromosome_tables,
    ignore_index=True,
)

unassigned_bins = binned_data[
    binned_data["region"].isna()
].copy()

unassigned_bins.to_csv(
    OUTPUT_DIR
    / "H3K9_delta_TE_unassigned_1kb_bins.csv",
    index=False,
)

binned_data = binned_data[
    binned_data["region"].isin(
        COMPARTMENT_ORDER
    )
].copy()

if binned_data.empty:
    raise ValueError(
        "No 1-kb bins were assigned to Eu, PCH, or Telo."
    )

binned_data["region"] = pd.Categorical(
    binned_data["region"],
    categories=COMPARTMENT_ORDER,
    ordered=True,
)

binned_data["h3k9_delta"] = (
    binned_data["h3k9_sperm"]
    - binned_data["h3k9_wing"]
)

binned_data["h3k9_log2fc"] = np.log2(
    (
        binned_data["h3k9_sperm"]
        + H3K9_PSEUDOCOUNT
    )
    /
    (
        binned_data["h3k9_wing"]
        + H3K9_PSEUDOCOUNT
    )
)

binned_data.to_csv(
    OUTPUT_DIR
    / "H3K9_delta_TE_binned_1kb.csv",
    index=False,
)


# ==============================================================================
# 6. Spearman correlations between delta H3K9me2 and TE expression
# ==============================================================================

correlation_rows: list[dict[str, object]] = []

correlation_scopes = [
    ("Genome-wide", binned_data)
]

correlation_scopes.extend(
    (
        compartment,
        binned_data[
            binned_data["region"]
            == compartment
        ],
    )
    for compartment in COMPARTMENT_ORDER
)

for scope, scope_data in correlation_scopes:
    n_bins, rho, p_value = safe_spearman(
        data=scope_data,
        x_column="h3k9_delta",
        y_column="te_sperm_cpm",
    )

    correlation_rows.append(
        {
            "scope": scope,
            "n_bins": n_bins,
            "spearman_rho": rho,
            "p_value": p_value,
        }
    )

correlation_table = pd.DataFrame(
    correlation_rows
)

correlation_table.to_csv(
    OUTPUT_DIR
    / "H3K9_delta_TE_Spearman_correlations_1kb.csv",
    index=False,
)


# ==============================================================================
# 7. Eu and PCH ventile summaries
# ==============================================================================

ventile_tables: list[pd.DataFrame] = []

for compartment in ("Eu", "PCH"):
    compartment_data = binned_data[
        binned_data["region"]
        == compartment
    ].copy()

    finite_mask = (
        np.isfinite(
            compartment_data["h3k9_delta"]
        )
        & np.isfinite(
            compartment_data["te_sperm_cpm"]
        )
    )

    compartment_data = compartment_data[
        finite_mask
    ].copy()

    if len(compartment_data) < N_VENTILES:
        print(
            f"Warning: {compartment} contains fewer than "
            f"{N_VENTILES} valid bins; ventiles were not calculated.",
            file=sys.stderr,
        )
        continue

    compartment_data["H3K9_delta_ventile"] = (
        pd.qcut(
            compartment_data["h3k9_delta"],
            q=N_VENTILES,
            labels=False,
            duplicates="drop",
        )
        + 1
    )

    ventile_summary = (
        compartment_data
        .groupby(
            "H3K9_delta_ventile",
            observed=True,
        )
        .agg(
            n_bins=(
                "h3k9_delta",
                "size",
            ),
            minimum_delta_H3K9me2=(
                "h3k9_delta",
                "min",
            ),
            maximum_delta_H3K9me2=(
                "h3k9_delta",
                "max",
            ),
            mean_delta_H3K9me2=(
                "h3k9_delta",
                "mean",
            ),
            median_delta_H3K9me2=(
                "h3k9_delta",
                "median",
            ),
            mean_TE_CPM=(
                "te_sperm_cpm",
                "mean",
            ),
            median_TE_CPM=(
                "te_sperm_cpm",
                "median",
            ),
        )
        .reset_index()
    )

    ventile_summary.insert(
        0,
        "compartment",
        compartment,
    )

    ventile_tables.append(
        ventile_summary
    )

if ventile_tables:
    ventile_table = pd.concat(
        ventile_tables,
        ignore_index=True,
    )
else:
    ventile_table = pd.DataFrame(
        columns=[
            "compartment",
            "H3K9_delta_ventile",
            "n_bins",
            "minimum_delta_H3K9me2",
            "maximum_delta_H3K9me2",
            "mean_delta_H3K9me2",
            "median_delta_H3K9me2",
            "mean_TE_CPM",
            "median_TE_CPM",
        ]
    )

ventile_table.to_csv(
    OUTPUT_DIR
    / "H3K9_delta_TE_ventiles_Eu_PCH_1kb.csv",
    index=False,
)


# ==============================================================================
# 8. Save parameters and software versions
# ==============================================================================

analysis_parameters = pd.DataFrame(
    {
        "parameter": [
            "Spermatocyte H3K9me2 bigWig",
            "Wing H3K9me2 bigWig",
            "TE-expression bigWig",
            "Nominal bin size",
            "H3K9me2 pseudocount",
            "Target chromosomes",
            "Delta definition",
            "Correlation",
            "Ventiles",
        ],
        "value": [
            str(H3K9_SPERMATOCYTE_BIGWIG),
            str(H3K9_WING_BIGWIG),
            str(TE_EXPRESSION_BIGWIG),
            BIN_SIZE_BP,
            H3K9_PSEUDOCOUNT,
            ", ".join(TARGET_CHROMOSOMES),
            "Spermatocyte RPGC - wing RPGC",
            "Two-sided Spearman correlation using all finite bins",
            (
                f"{N_VENTILES} equal-frequency groups calculated "
                "independently within Eu and PCH"
            ),
        ],
    }
)

analysis_parameters.to_csv(
    OUTPUT_DIR
    / "analysis_parameters.csv",
    index=False,
)

software_versions = pd.DataFrame(
    {
        "software": [
            "Python",
            "numpy",
            "pandas",
            "pyBigWig",
            "scipy",
        ],
        "version": [
            sys.version.replace("\n", " "),
            package_version("numpy"),
            package_version("pandas"),
            package_version("pyBigWig"),
            package_version("scipy"),
        ],
    }
)

software_versions.to_csv(
    OUTPUT_DIR
    / "software_versions.csv",
    index=False,
)


# ==============================================================================
# 9. Final summary
# ==============================================================================

print("\nH3K9me2–TE 1-kb correlation analysis complete.")
print(f"Results were written to: {OUTPUT_DIR}")

print("\nSpearman correlations:")
print(
    correlation_table.to_string(
        index=False
    )
)

print("\nNumber of retained 1-kb bins:")
print(
    binned_data[
        "region"
    ].value_counts(
        sort=False
    ).to_string()
)
