from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import openpyxl


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = Path.home() / "fragilemap_sc"

BREAK_FILE = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "human_break_regions_enriched.csv"
)

FEATURE_WORKBOOK = (
    PROJECT_DIR
    / "data"
    / "raw"
    / "41467_2026_76451_MOESM4_ESM.xlsx"
)

OUTPUT_FILE = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "u2os_break_regions_multifeature.csv"
)


# ============================================================
# U2OS feature columns in Supplementary Data 2
#
# 1-based spreadsheet columns:
#
# 8  RT slope / Untreated
# 12 MiDAS-seq
# 13 lateS/G2/M-seq
# 14 A/B compartment
# 15 GRO-seq
# 16 FANCD2-seq
# 17 DRIP-seq
# 18 G-quadruplex
# 19 NS-seq
#
# Converted below to zero-based Python indices.
# ============================================================

FEATURE_COLUMNS = {
    "rt_slope_untreated": 7,
    "midas": 11,
    "lateS_G2M": 12,
    "ab_compartment": 13,
    "gro_seq": 14,
    "fancd2": 15,
    "drip_seq": 16,
    "g_quadruplex": 17,
    "ns_seq": 18
}


def normalize_chrom(value):

    value = str(value).strip()

    if value.lower().startswith("chr"):
        return "chr" + value[3:]

    return "chr" + value


def is_number(value):

    try:
        float(value)
        return True

    except (TypeError, ValueError):
        return False


# ============================================================
# Load U2OS break regions
# ============================================================

breaks = pd.read_csv(
    BREAK_FILE
)

breaks = breaks[
    breaks["cell_type"] == "U2OS"
].copy()

breaks["chrom"] = (
    breaks["chrom"]
    .astype(str)
    .map(normalize_chrom)
)

assert len(breaks) == 195


# ============================================================
# Load U2OS Supplementary Data 2
# ============================================================

print("\nU2OS MULTI-FEATURE INTEGRATION")
print("=" * 80)

wb = openpyxl.load_workbook(
    FEATURE_WORKBOOK,
    read_only=True,
    data_only=True
)

ws = wb["U2OS"]


# ============================================================
# Store feature bins by chromosome
# ============================================================

feature_bins = defaultdict(list)

n_genomic_rows = 0

for row in ws.iter_rows(
    values_only=True
):

    if len(row) < 19:
        continue

    chrom = row[0]
    start = row[1]
    end = row[2]

    if chrom is None:
        continue

    chrom_text = str(chrom).strip()

    if not chrom_text.lower().startswith(
        "chr"
    ):
        continue

    if (
        not is_number(start)
        or not is_number(end)
    ):
        continue

    chrom_text = normalize_chrom(
        chrom_text
    )

    start = int(float(start))
    end = int(float(end))

    feature_values = {}

    for feature, col_idx in (
        FEATURE_COLUMNS.items()
    ):

        value = row[col_idx]

        if (
            value is not None
            and is_number(value)
        ):
            feature_values[
                feature
            ] = float(value)

        else:
            feature_values[
                feature
            ] = np.nan

    feature_bins[
        chrom_text
    ].append(
        (
            start,
            end,
            feature_values
        )
    )

    n_genomic_rows += 1


wb.close()


# Sort genomic bins
for chrom in feature_bins:

    feature_bins[
        chrom
    ].sort(
        key=lambda x: x[0]
    )


print(
    "Genomic rows loaded:",
    n_genomic_rows
)


# ============================================================
# Overlap-weighted mean
# ============================================================

def weighted_feature_mean(
    chrom,
    break_start,
    break_end,
    feature
):

    bins = feature_bins.get(
        chrom,
        []
    )

    weighted_sum = 0.0
    total_overlap = 0

    for (
        bin_start,
        bin_end,
        values
    ) in bins:

        if bin_start >= break_end:
            break

        if bin_end <= break_start:
            continue

        value = values[
            feature
        ]

        if pd.isna(value):
            continue

        overlap_start = max(
            break_start,
            bin_start
        )

        overlap_end = min(
            break_end,
            bin_end
        )

        overlap_bp = (
            overlap_end
            - overlap_start
        )

        if overlap_bp <= 0:
            continue

        weighted_sum += (
            value
            * overlap_bp
        )

        total_overlap += (
            overlap_bp
        )

    if total_overlap == 0:
        return np.nan

    return (
        weighted_sum
        / total_overlap
    )


# ============================================================
# Attach all genomic features
# ============================================================

for feature in FEATURE_COLUMNS:

    values = []

    for _, row in breaks.iterrows():

        value = weighted_feature_mean(
            row["chrom"],
            int(row["start"]),
            int(row["end"]),
            feature
        )

        values.append(value)

    breaks[feature] = values


# ============================================================
# Save
# ============================================================

breaks.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# Coverage summary
# ============================================================

print("\nBREAK-LEVEL FEATURE COVERAGE")
print("=" * 80)

coverage = {}

for feature in FEATURE_COLUMNS:

    observed = int(
        breaks[feature]
        .notna()
        .sum()
    )

    missing = int(
        breaks[feature]
        .isna()
        .sum()
    )

    coverage[
        feature
    ] = observed

    print(
        f"{feature:<22} "
        f"{observed:>3}/195 "
        f"({observed / 195 * 100:5.1f}%) "
        f"| missing={missing}"
    )


# ============================================================
# Validation
# ============================================================

expected_coverage = {
    "rt_slope_untreated": 195,
    "midas": 195,
    "lateS_G2M": 195,
    "ab_compartment": 193,
    "gro_seq": 194,
    "fancd2": 195,
    "drip_seq": 195,
    "g_quadruplex": 195,
    "ns_seq": 195
}


assert len(breaks) == 195

for (
    feature,
    expected
) in expected_coverage.items():

    observed = coverage[
        feature
    ]

    assert observed == expected, (
        f"{feature}: expected "
        f"{expected} non-missing values, "
        f"found {observed}"
    )


print("\nValidation: PASS")
print(
    "Integrated U2OS break regions:",
    len(breaks)
)

print("\nSaved:")
print(OUTPUT_FILE)
