from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import (
    chi2_contingency,
    kruskal,
    mannwhitneyu,
    spearmanr
)
from statsmodels.stats.multitest import multipletests


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = Path.home() / "fragilemap_sc"

INPUT_FILE = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "human_break_regions_rt_classified.csv"
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
# Load data
# ============================================================

df = pd.read_csv(INPUT_FILE)

RT_COLUMN = "mean_pseudobulk_rt_untreated"

CELL_ORDER = [
    "U2OS",
    "RPE1",
    "BJ"
]

RT_ORDER = [
    "Early",
    "Mid",
    "Late"
]


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
assert df[RT_COLUMN].isna().sum() == 0
assert df["recurrent"].isna().sum() == 0


# ============================================================
# Helper: rank-biserial effect size
# ============================================================

def rank_biserial_from_u(
    u_statistic,
    n_group1,
    n_group2
):
    """
    Positive value means group 1 tends to have
    larger values than group 2.
    """

    return (
        (2 * u_statistic)
        / (n_group1 * n_group2)
        - 1
    )


# ============================================================
# 1. RT CLASS COUNTS
# ============================================================

rt_counts = (
    df.groupby(
        ["cell_type", "rt_class"]
    )
    .size()
    .unstack(fill_value=0)
    .reindex(
        index=CELL_ORDER,
        columns=RT_ORDER
    )
)

rt_counts.to_csv(
    TABLE_DIR
    / "rt_class_counts.csv"
)


# ============================================================
# 2. OVERALL RT CLASS CHI-SQUARE
# ============================================================

chi2, p, dof, expected = (
    chi2_contingency(
        rt_counts.values
    )
)

overall_rt_class = pd.DataFrame([{
    "test": "Chi-square",
    "comparison": "U2OS vs RPE1 vs BJ",
    "chi_square": chi2,
    "df": dof,
    "p_value": p,
    "minimum_expected_count": expected.min()
}])

overall_rt_class.to_csv(
    TABLE_DIR
    / "rt_class_overall_test.csv",
    index=False
)


# ============================================================
# 3. PAIRWISE RT CLASS CHI-SQUARE
# ============================================================

pairs = [
    ("U2OS", "RPE1"),
    ("U2OS", "BJ"),
    ("RPE1", "BJ")
]

pairwise_class_rows = []

for cell1, cell2 in pairs:

    subtable = rt_counts.loc[
        [cell1, cell2],
        RT_ORDER
    ]

    chi2_pair, p_pair, dof_pair, expected_pair = (
        chi2_contingency(
            subtable.values
        )
    )

    pairwise_class_rows.append({
        "comparison": f"{cell1} vs {cell2}",
        "chi_square": chi2_pair,
        "df": dof_pair,
        "raw_p": p_pair,
        "minimum_expected_count": (
            expected_pair.min()
        )
    })


pairwise_class = pd.DataFrame(
    pairwise_class_rows
)

pairwise_class["holm_adjusted_p"] = (
    multipletests(
        pairwise_class["raw_p"],
        method="holm"
    )[1]
)

pairwise_class.to_csv(
    TABLE_DIR
    / "rt_class_pairwise_tests.csv",
    index=False
)


# ============================================================
# 4. CONTINUOUS RT SUMMARY
# ============================================================

continuous_summary_rows = []

for cell in CELL_ORDER:

    values = df.loc[
        df["cell_type"] == cell,
        RT_COLUMN
    ].dropna()

    continuous_summary_rows.append({
        "cell_type": cell,
        "n": len(values),
        "mean_rt": values.mean(),
        "median_rt": values.median(),
        "std_rt": values.std(),
        "iqr_rt": (
            values.quantile(0.75)
            - values.quantile(0.25)
        )
    })


continuous_summary = pd.DataFrame(
    continuous_summary_rows
)

continuous_summary.to_csv(
    TABLE_DIR
    / "rt_continuous_summary.csv",
    index=False
)


# ============================================================
# 5. KRUSKAL-WALLIS
# ============================================================

groups = [
    df.loc[
        df["cell_type"] == cell,
        RT_COLUMN
    ].dropna()
    for cell in CELL_ORDER
]

kruskal_h, kruskal_p = (
    kruskal(*groups)
)

pd.DataFrame([{
    "test": "Kruskal-Wallis",
    "H": kruskal_h,
    "p_value": kruskal_p
}]).to_csv(
    TABLE_DIR
    / "rt_continuous_overall_test.csv",
    index=False
)


# ============================================================
# 6. PAIRWISE CONTINUOUS RT
# ============================================================

pairwise_rt_rows = []

for cell1, cell2 in pairs:

    x = df.loc[
        df["cell_type"] == cell1,
        RT_COLUMN
    ].dropna()

    y = df.loc[
        df["cell_type"] == cell2,
        RT_COLUMN
    ].dropna()

    u, p_pair = mannwhitneyu(
        x,
        y,
        alternative="two-sided"
    )

    effect = rank_biserial_from_u(
        u,
        len(x),
        len(y)
    )

    pairwise_rt_rows.append({
        "comparison": f"{cell1} vs {cell2}",
        "n_1": len(x),
        "n_2": len(y),
        "median_1": x.median(),
        "median_2": y.median(),
        "rank_biserial_effect": effect,
        "raw_p": p_pair
    })


pairwise_rt = pd.DataFrame(
    pairwise_rt_rows
)

pairwise_rt["holm_adjusted_p"] = (
    multipletests(
        pairwise_rt["raw_p"],
        method="holm"
    )[1]
)

pairwise_rt.to_csv(
    TABLE_DIR
    / "rt_pairwise_tests.csv",
    index=False
)


# ============================================================
# 7. RECURRENCE × RT CLASS
# ============================================================

recurrence_class_rows = []

for cell in CELL_ORDER:

    subset = df[
        df["cell_type"] == cell
    ].copy()

    table = pd.crosstab(
        subset["recurrent"],
        subset["rt_class"]
    ).reindex(
        index=[False, True],
        columns=RT_ORDER,
        fill_value=0
    )

    chi2_r, p_r, dof_r, expected_r = (
        chi2_contingency(
            table.values
        )
    )

    recurrence_class_rows.append({
        "cell_type": cell,
        "chi_square": chi2_r,
        "df": dof_r,
        "raw_p": p_r,
        "minimum_expected_count": (
            expected_r.min()
        ),
        "single_early": table.loc[
            False, "Early"
        ],
        "single_mid": table.loc[
            False, "Mid"
        ],
        "single_late": table.loc[
            False, "Late"
        ],
        "recurrent_early": table.loc[
            True, "Early"
        ],
        "recurrent_mid": table.loc[
            True, "Mid"
        ],
        "recurrent_late": table.loc[
            True, "Late"
        ]
    })


recurrence_class = pd.DataFrame(
    recurrence_class_rows
)

recurrence_class["holm_adjusted_p"] = (
    multipletests(
        recurrence_class["raw_p"],
        method="holm"
    )[1]
)

recurrence_class.to_csv(
    TABLE_DIR
    / "recurrence_rt_class_tests.csv",
    index=False
)


# ============================================================
# 8. RECURRENCE × CONTINUOUS RT
# ============================================================

recurrence_continuous_rows = []

for cell in CELL_ORDER:

    subset = df[
        df["cell_type"] == cell
    ]

    recurrent = subset.loc[
        subset["recurrent"],
        RT_COLUMN
    ].dropna()

    single = subset.loc[
        ~subset["recurrent"],
        RT_COLUMN
    ].dropna()

    u, p_r = mannwhitneyu(
        recurrent,
        single,
        alternative="two-sided"
    )

    effect = rank_biserial_from_u(
        u,
        len(recurrent),
        len(single)
    )

    recurrence_continuous_rows.append({
        "cell_type": cell,
        "single_n": len(single),
        "recurrent_n": len(recurrent),
        "single_median_rt": (
            single.median()
        ),
        "recurrent_median_rt": (
            recurrent.median()
        ),
        "rank_biserial_effect": (
            effect
        ),
        "raw_p": p_r
    })


recurrence_continuous = pd.DataFrame(
    recurrence_continuous_rows
)

recurrence_continuous[
    "holm_adjusted_p"
] = multipletests(
    recurrence_continuous["raw_p"],
    method="holm"
)[1]

recurrence_continuous.to_csv(
    TABLE_DIR
    / "recurrence_rt_continuous_tests.csv",
    index=False
)


# ============================================================
# 9. WIDTH × RT SPEARMAN
# ============================================================

width_rt_rows = []

for cell in CELL_ORDER:

    subset = df[
        df["cell_type"] == cell
    ].dropna(
        subset=[
            "region_width_kb",
            RT_COLUMN
        ]
    )

    rho, p_s = spearmanr(
        subset["region_width_kb"],
        subset[RT_COLUMN]
    )

    width_rt_rows.append({
        "cell_type": cell,
        "n": len(subset),
        "spearman_rho": rho,
        "p_value": p_s
    })


width_rt = pd.DataFrame(
    width_rt_rows
)

width_rt.to_csv(
    TABLE_DIR
    / "width_rt_spearman.csv",
    index=False
)


# ============================================================
# VALIDATION
# ============================================================

print("\nREPLICATION-TIMING STATISTICAL ANALYSIS")
print("=" * 75)

print(
    f"Overall RT class chi-square: "
    f"{chi2:.4f}"
)

print(
    f"Overall RT class p-value: "
    f"{p:.6f}"
)

print(
    f"\nContinuous RT Kruskal H: "
    f"{kruskal_h:.4f}"
)

print(
    f"Continuous RT Kruskal p: "
    f"{kruskal_p:.6f}"
)

print(
    "\nRECURRENCE × RT CLASS"
)

for _, row in recurrence_class.iterrows():

    print(
        f"{row['cell_type']:<5} | "
        f"raw p={row['raw_p']:.6f} | "
        f"Holm p="
        f"{row['holm_adjusted_p']:.6f}"
    )


# Validate previously established results
assert np.isclose(
    chi2,
    10.7489,
    atol=0.001
)

assert np.isclose(
    p,
    0.029536,
    atol=0.00001
)

assert np.isclose(
    kruskal_h,
    3.864703198936,
    atol=1e-9
)

assert np.isclose(
    kruskal_p,
    0.144807269079,
    atol=1e-9
)

u2os_result = recurrence_class.loc[
    recurrence_class["cell_type"]
    == "U2OS"
].iloc[0]

assert np.isclose(
    u2os_result[
        "holm_adjusted_p"
    ],
    0.01596,
    atol=0.0001
)

print("\nValidation: PASS")

print("\nGenerated tables:")

for path in [
    "rt_class_counts.csv",
    "rt_class_overall_test.csv",
    "rt_class_pairwise_tests.csv",
    "rt_continuous_summary.csv",
    "rt_continuous_overall_test.csv",
    "rt_pairwise_tests.csv",
    "recurrence_rt_class_tests.csv",
    "recurrence_rt_continuous_tests.csv",
    "width_rt_spearman.csv"
]:
    print(" -", TABLE_DIR / path)
