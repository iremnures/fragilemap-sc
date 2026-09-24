from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = Path.home() / "fragilemap_sc"

INPUT_FILE = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "human_break_regions_enriched.csv"
)

TABLE_DIR = (
    PROJECT_DIR
    / "results"
    / "tables"
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# hg19 chromosome lengths
# ============================================================

CHR_LENGTHS = {
    "chr1": 249250621,
    "chr2": 243199373,
    "chr3": 198022430,
    "chr4": 191154276,
    "chr5": 180915260,
    "chr6": 171115067,
    "chr7": 159138663,
    "chr8": 146364022,
    "chr9": 141213431,
    "chr10": 135534747,
    "chr11": 135006516,
    "chr12": 133851895,
    "chr13": 115169878,
    "chr14": 107349540,
    "chr15": 102531392,
    "chr16": 90354753,
    "chr17": 81195210,
    "chr18": 78077248,
    "chr19": 59128983,
    "chr20": 63025520,
    "chr21": 48129895,
    "chr22": 51304566,
    "chrX": 155270560,
    "chrY": 59373566
}

CHR_ORDER = list(CHR_LENGTHS.keys())

CELL_ORDER = [
    "U2OS",
    "RPE1",
    "BJ"
]


# ============================================================
# Load data
# ============================================================

df = pd.read_csv(INPUT_FILE)

assert len(df) == 644

unexpected = sorted(
    set(df["chrom"])
    - set(CHR_LENGTHS)
)

assert not unexpected, (
    f"Unexpected chromosome labels: {unexpected}"
)


# ============================================================
# Count break regions per chromosome
# ============================================================

rows = []

for cell in CELL_ORDER:

    sub = df[
        df["cell_type"] == cell
    ]

    for chrom in CHR_ORDER:

        count = int(
            (sub["chrom"] == chrom).sum()
        )

        length_bp = CHR_LENGTHS[chrom]

        density_per_100mb = (
            count
            / length_bp
            * 100_000_000
        )

        rows.append({
            "cell_type": cell,
            "chromosome": chrom,
            "break_region_count": count,
            "chromosome_length_bp": length_bp,
            "density_per_100mb": density_per_100mb
        })


long_df = pd.DataFrame(rows)

long_output = (
    TABLE_DIR
    / "chromosome_break_density_long.csv"
)

long_df.to_csv(
    long_output,
    index=False
)


# ============================================================
# Wide table
# ============================================================

wide_df = (
    long_df.pivot(
        index="chromosome",
        columns="cell_type",
        values="density_per_100mb"
    )
    .reindex(CHR_ORDER)
    .reset_index()
)

wide_df = wide_df[
    [
        "chromosome",
        "U2OS",
        "RPE1",
        "BJ"
    ]
]

wide_output = (
    TABLE_DIR
    / "chromosome_break_density_wide.csv"
)

wide_df.to_csv(
    wide_output,
    index=False
)


# ============================================================
# Correlations — all primary chromosomes
# ============================================================

pairs = [
    ("U2OS", "RPE1"),
    ("U2OS", "BJ"),
    ("RPE1", "BJ")
]

corr_rows = []

for cell1, cell2 in pairs:

    rho, p = spearmanr(
        wide_df[cell1],
        wide_df[cell2]
    )

    corr_rows.append({
        "dataset": "All primary chromosomes",
        "comparison": f"{cell1} vs {cell2}",
        "n_chromosomes": 24,
        "spearman_rho": rho,
        "p_value": p
    })


corr_df = pd.DataFrame(
    corr_rows
)

corr_output = (
    TABLE_DIR
    / "chromosome_density_correlations.csv"
)

corr_df.to_csv(
    corr_output,
    index=False
)


# ============================================================
# Autosomes-only sensitivity
# ============================================================

autosomes = [
    f"chr{i}"
    for i in range(1, 23)
]

auto_df = wide_df[
    wide_df["chromosome"].isin(
        autosomes
    )
].copy()

sensitivity_rows = []

for cell1, cell2 in pairs:

    rho, p = spearmanr(
        auto_df[cell1],
        auto_df[cell2]
    )

    sensitivity_rows.append({
        "dataset": "Autosomes only",
        "comparison": f"{cell1} vs {cell2}",
        "n_chromosomes": 22,
        "spearman_rho": rho,
        "p_value": p
    })


sensitivity_df = pd.DataFrame(
    sensitivity_rows
)

sensitivity_output = (
    TABLE_DIR
    / "chromosome_density_correlation_sensitivity.csv"
)

sensitivity_df.to_csv(
    sensitivity_output,
    index=False
)


# ============================================================
# Validation
# ============================================================

print("\nCHROMOSOME-LEVEL BREAK-DENSITY ANALYSIS")
print("=" * 80)

for cell in CELL_ORDER:

    total_regions = int(
        long_df.loc[
            long_df["cell_type"] == cell,
            "break_region_count"
        ].sum()
    )

    print(
        f"{cell:<5} total regions = {total_regions}"
    )

expected_totals = {
    "U2OS": 195,
    "RPE1": 284,
    "BJ": 165
}

for cell, expected in expected_totals.items():

    observed = int(
        long_df.loc[
            long_df["cell_type"] == cell,
            "break_region_count"
        ].sum()
    )

    assert observed == expected


print("\nTOP DENSITY CHROMOSOMES")
print("=" * 80)

for cell in CELL_ORDER:

    top5 = (
        long_df[
            long_df["cell_type"] == cell
        ]
        .sort_values(
            "density_per_100mb",
            ascending=False
        )
        .head(5)
    )

    print(f"\n{cell}")

    for _, row in top5.iterrows():

        print(
            f"{row['chromosome']:<5} "
            f"{row['density_per_100mb']:.2f}"
        )


print("\nALL PRIMARY CHROMOSOMES")
print("=" * 80)

for _, row in corr_df.iterrows():

    print(
        f"{row['comparison']:<15} "
        f"rho={row['spearman_rho']:.4f} "
        f"p={row['p_value']:.6f}"
    )


print("\nAUTOSOMES-ONLY SENSITIVITY")
print("=" * 80)

for _, row in sensitivity_df.iterrows():

    print(
        f"{row['comparison']:<15} "
        f"rho={row['spearman_rho']:.4f} "
        f"p={row['p_value']:.6f}"
    )


# Previously established core values
expected_all = {
    "U2OS vs RPE1": (0.0391, 0.855948),
    "U2OS vs BJ": (-0.1000, 0.641920),
    "RPE1 vs BJ": (0.2509, 0.236933)
}

for comparison, (
    expected_rho,
    expected_p
) in expected_all.items():

    row = corr_df.loc[
        corr_df["comparison"]
        == comparison
    ].iloc[0]

    assert np.isclose(
        row["spearman_rho"],
        expected_rho,
        atol=0.001
    )

    assert np.isclose(
        row["p_value"],
        expected_p,
        atol=0.001
    )


expected_auto = {
    "U2OS vs RPE1": (-0.1474, 0.512799),
    "U2OS vs BJ": (-0.3552, 0.104797),
    "RPE1 vs BJ": (0.1101, 0.625698)
}

for comparison, (
    expected_rho,
    expected_p
) in expected_auto.items():

    row = sensitivity_df.loc[
        sensitivity_df["comparison"]
        == comparison
    ].iloc[0]

    assert np.isclose(
        row["spearman_rho"],
        expected_rho,
        atol=0.001
    )

    assert np.isclose(
        row["p_value"],
        expected_p,
        atol=0.001
    )


print("\nValidation: PASS")

print("\nSaved:")
print(long_output)
print(wide_output)
print(corr_output)
print(sensitivity_output)
