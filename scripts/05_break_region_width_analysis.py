from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


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

OUTPUT_FILE = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "recurrent_vs_single_width_tests.csv"
)


# ============================================================
# Load data
# ============================================================

df = pd.read_csv(INPUT_FILE)


def parse_bool(series):
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map({
            "true": True,
            "false": False
        })
    )


df["recurrent"] = parse_bool(
    df["recurrent"]
)

assert len(df) == 644
assert df["recurrent"].isna().sum() == 0
assert (df["region_width_kb"] > 0).all()


CELL_ORDER = [
    "U2OS",
    "RPE1",
    "BJ"
]


# ============================================================
# Rank-biserial effect size
# ============================================================

def rank_biserial_from_u(
    u,
    n1,
    n2
):
    """
    Positive values indicate larger values
    in group 1.
    """
    return (
        (2 * u)
        / (n1 * n2)
        - 1
    )


# ============================================================
# Recurrent vs single width
# ============================================================

rows = []

for cell in CELL_ORDER:

    subset = df[
        df["cell_type"] == cell
    ]

    recurrent = subset.loc[
        subset["recurrent"],
        "region_width_kb"
    ].dropna()

    single = subset.loc[
        ~subset["recurrent"],
        "region_width_kb"
    ].dropna()

    u, p = mannwhitneyu(
        recurrent,
        single,
        alternative="two-sided"
    )

    effect = rank_biserial_from_u(
        u,
        len(recurrent),
        len(single)
    )

    rows.append({
        "cell_type": cell,

        "single_n":
            len(single),

        "recurrent_n":
            len(recurrent),

        "median_width_single_kb":
            single.median(),

        "median_width_recurrent_kb":
            recurrent.median(),

        "mean_width_single_kb":
            single.mean(),

        "mean_width_recurrent_kb":
            recurrent.mean(),

        "rank_biserial_effect":
            effect,

        "mann_whitney_u":
            u,

        "raw_p":
            p
    })


results = pd.DataFrame(rows)

results[
    "holm_adjusted_p"
] = multipletests(
    results["raw_p"],
    method="holm"
)[1]


results.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# Validation
# ============================================================

print("\nBREAK-REGION WIDTH ANALYSIS")
print("=" * 80)

for _, row in results.iterrows():

    print(
        f"{row['cell_type']:<5} | "
        f"single n={int(row['single_n']):<3} | "
        f"recurrent n={int(row['recurrent_n']):<2} | "
        f"median "
        f"{row['median_width_single_kb']:.2f} -> "
        f"{row['median_width_recurrent_kb']:.2f} kb | "
        f"RBC={row['rank_biserial_effect']:.4f} | "
        f"Holm p={row['holm_adjusted_p']:.3e}"
    )


# Expected group sizes
expected_n = {
    "U2OS": (151, 44),
    "RPE1": (210, 74),
    "BJ": (123, 42)
}

for cell, (
    expected_single,
    expected_recurrent
) in expected_n.items():

    r = results.loc[
        results["cell_type"] == cell
    ].iloc[0]

    assert int(
        r["single_n"]
    ) == expected_single

    assert int(
        r["recurrent_n"]
    ) == expected_recurrent


# Previously established medians
expected_medians = {
    "U2OS": (268.83, 746.72),
    "RPE1": (435.50, 1201.49),
    "BJ": (298.20, 704.41)
}

for cell, (
    expected_single,
    expected_recurrent
) in expected_medians.items():

    r = results.loc[
        results["cell_type"] == cell
    ].iloc[0]

    assert np.isclose(
        r["median_width_single_kb"],
        expected_single,
        atol=0.02
    )

    assert np.isclose(
        r["median_width_recurrent_kb"],
        expected_recurrent,
        atol=0.02
    )


# Previously established rank-biserial effects
expected_effects = {
    "U2OS": 0.6141,
    "RPE1": 0.6568,
    "BJ": 0.6703
}

for cell, expected in (
    expected_effects.items()
):

    observed = results.loc[
        results["cell_type"] == cell,
        "rank_biserial_effect"
    ].iloc[0]

    assert np.isclose(
        observed,
        expected,
        atol=0.002
    )


# All three comparisons should remain
# significant after Holm correction.
assert (
    results["holm_adjusted_p"] < 0.001
).all()


print("\nValidation: PASS")
print("Saved:", OUTPUT_FILE)
