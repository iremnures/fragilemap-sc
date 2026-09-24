from pathlib import Path
import math

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm


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

OUTPUT_DISTANCE = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "human_break_regions_overlap_distance.csv"
)

TABLE_DIR = (
    PROJECT_DIR
    / "results"
    / "tables"
)

OUTPUT_TESTS = (
    TABLE_DIR
    / "overlap_edge_distance_tests.csv"
)

OUTPUT_LOGISTIC = (
    TABLE_DIR
    / "width_adjusted_overlap_logistic.csv"
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True
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


CELL_ORDER = [
    "U2OS",
    "RPE1",
    "BJ"
]


# ============================================================
# Interval edge distance
# ============================================================

def interval_edge_distance(
    start1,
    end1,
    start2,
    end2
):
    """
    Distance between two genomic intervals.

    Returns 0 when intervals overlap.
    Otherwise returns the gap between their nearest edges.
    """

    overlap = (
        max(start1, start2)
        < min(end1, end2)
    )

    if overlap:
        return 0

    if end1 <= start2:
        return start2 - end1

    return start1 - end2


# ============================================================
# Find nearest break region in another human cell line
# ============================================================

results = []

for idx, row in df.iterrows():

    candidates = df[
        (df["cell_type"] != row["cell_type"])
        &
        (df["chrom"] == row["chrom"])
    ]

    nearest_distance = None
    nearest_cell_line = None

    for _, other in candidates.iterrows():

        distance = interval_edge_distance(
            row["start"],
            row["end"],
            other["start"],
            other["end"]
        )

        if (
            nearest_distance is None
            or distance < nearest_distance
        ):
            nearest_distance = distance
            nearest_cell_line = other["cell_type"]

            if nearest_distance == 0:
                # We already know this region overlaps
                # another cell line.
                break

    # All primary human chromosomes represented here
    # should have a comparison region.
    if nearest_distance is None:
        nearest_distance = np.nan

    new_row = row.to_dict()

    new_row[
        "nearest_other_cell_line"
    ] = nearest_cell_line

    new_row[
        "nearest_edge_distance_bp"
    ] = nearest_distance

    new_row[
        "nearest_edge_distance_kb"
    ] = (
        nearest_distance / 1000
        if pd.notna(nearest_distance)
        else np.nan
    )

    new_row[
        "overlap_other_line"
    ] = (
        nearest_distance == 0
        if pd.notna(nearest_distance)
        else False
    )

    new_row[
        "within_1mb_other_line"
    ] = (
        nearest_distance <= 1_000_000
        if pd.notna(nearest_distance)
        else False
    )

    results.append(new_row)


distance_df = pd.DataFrame(results)

distance_df.to_csv(
    OUTPUT_DISTANCE,
    index=False
)


# ============================================================
# Helper functions
# ============================================================

def rank_biserial_from_u(
    u,
    n1,
    n2
):
    return (
        (2 * u)
        / (n1 * n2)
        - 1
    )


def odds_ratio_ci(
    a,
    b,
    c,
    d
):
    """
    Wald 95% CI for odds ratio.

    Table:
                 Event   No event
    recurrent      a        b
    single         c        d
    """

    values = [a, b, c, d]

    # Continuity correction only if required.
    if any(v == 0 for v in values):
        a, b, c, d = [
            v + 0.5
            for v in values
        ]

    odds_ratio = (
        a * d
        / (b * c)
    )

    se = math.sqrt(
        1/a
        + 1/b
        + 1/c
        + 1/d
    )

    log_or = math.log(
        odds_ratio
    )

    lower = math.exp(
        log_or - 1.96 * se
    )

    upper = math.exp(
        log_or + 1.96 * se
    )

    return (
        odds_ratio,
        lower,
        upper
    )


# ============================================================
# Edge-distance / overlap statistics
# ============================================================

test_rows = []

for cell in CELL_ORDER:

    sub = distance_df[
        distance_df["cell_type"] == cell
    ].copy()

    single = sub.loc[
        ~sub["recurrent"],
        "nearest_edge_distance_kb"
    ].dropna()

    recurrent = sub.loc[
        sub["recurrent"],
        "nearest_edge_distance_kb"
    ].dropna()

    # Single first so a positive effect means
    # single regions tend to be farther away,
    # i.e. recurrent regions are closer.
    u, distance_p = mannwhitneyu(
        single,
        recurrent,
        alternative="two-sided"
    )

    distance_effect = (
        rank_biserial_from_u(
            u,
            len(single),
            len(recurrent)
        )
    )

    # --------------------------------------------------------
    # Direct overlap
    # --------------------------------------------------------

    rec_overlap = int(
        sub.loc[
            sub["recurrent"],
            "overlap_other_line"
        ].sum()
    )

    rec_total = int(
        sub["recurrent"].sum()
    )

    rec_no_overlap = (
        rec_total
        - rec_overlap
    )

    single_overlap = int(
        sub.loc[
            ~sub["recurrent"],
            "overlap_other_line"
        ].sum()
    )

    single_total = int(
        (~sub["recurrent"]).sum()
    )

    single_no_overlap = (
        single_total
        - single_overlap
    )

    overlap_table = [
        [
            rec_overlap,
            rec_no_overlap
        ],
        [
            single_overlap,
            single_no_overlap
        ]
    ]

    fisher_or, overlap_p = (
        fisher_exact(
            overlap_table,
            alternative="two-sided"
        )
    )

    _, overlap_ci_low, overlap_ci_high = (
        odds_ratio_ci(
            rec_overlap,
            rec_no_overlap,
            single_overlap,
            single_no_overlap
        )
    )

    # --------------------------------------------------------
    # Within 1 Mb
    # --------------------------------------------------------

    rec_1mb = int(
        sub.loc[
            sub["recurrent"],
            "within_1mb_other_line"
        ].sum()
    )

    rec_not_1mb = (
        rec_total
        - rec_1mb
    )

    single_1mb = int(
        sub.loc[
            ~sub["recurrent"],
            "within_1mb_other_line"
        ].sum()
    )

    single_not_1mb = (
        single_total
        - single_1mb
    )

    one_mb_table = [
        [
            rec_1mb,
            rec_not_1mb
        ],
        [
            single_1mb,
            single_not_1mb
        ]
    ]

    one_mb_or, one_mb_p = (
        fisher_exact(
            one_mb_table,
            alternative="two-sided"
        )
    )

    _, one_mb_ci_low, one_mb_ci_high = (
        odds_ratio_ci(
            rec_1mb,
            rec_not_1mb,
            single_1mb,
            single_not_1mb
        )
    )

    test_rows.append({
        "cell_type": cell,

        "single_n": len(single),
        "recurrent_n": len(recurrent),

        "single_median_edge_distance_kb":
            single.median(),

        "recurrent_median_edge_distance_kb":
            recurrent.median(),

        "distance_rank_biserial_effect":
            distance_effect,

        "distance_raw_p":
            distance_p,

        "single_overlap_n":
            single_overlap,

        "recurrent_overlap_n":
            rec_overlap,

        "overlap_odds_ratio":
            fisher_or,

        "overlap_ci_low":
            overlap_ci_low,

        "overlap_ci_high":
            overlap_ci_high,

        "overlap_raw_p":
            overlap_p,

        "single_within_1mb_n":
            single_1mb,

        "recurrent_within_1mb_n":
            rec_1mb,

        "within_1mb_odds_ratio":
            one_mb_or,

        "within_1mb_ci_low":
            one_mb_ci_low,

        "within_1mb_ci_high":
            one_mb_ci_high,

        "within_1mb_raw_p":
            one_mb_p
    })


tests = pd.DataFrame(
    test_rows
)


# Separate Holm correction families
tests[
    "distance_holm_p"
] = multipletests(
    tests["distance_raw_p"],
    method="holm"
)[1]

tests[
    "overlap_holm_p"
] = multipletests(
    tests["overlap_raw_p"],
    method="holm"
)[1]

tests[
    "within_1mb_holm_p"
] = multipletests(
    tests["within_1mb_raw_p"],
    method="holm"
)[1]


tests.to_csv(
    OUTPUT_TESTS,
    index=False
)


# ============================================================
# Width-adjusted logistic regression
#
# Outcome:
#     overlap with another cell line
#
# Predictors:
#     recurrent status
#     standardized log10(region width)
# ============================================================

logistic_rows = []

for cell in CELL_ORDER:

    sub = distance_df[
        distance_df["cell_type"] == cell
    ].copy()

    sub[
        "recurrent_binary"
    ] = sub["recurrent"].astype(int)

    sub[
        "overlap_binary"
    ] = sub[
        "overlap_other_line"
    ].astype(int)

    sub[
        "log10_width"
    ] = np.log10(
        sub["region_width_kb"]
    )

    sub[
        "log10_width_z"
    ] = (
        sub["log10_width"]
        - sub["log10_width"].mean()
    ) / sub["log10_width"].std(
        ddof=0
    )

    X = sub[[
        "recurrent_binary",
        "log10_width_z"
    ]]

    X = sm.add_constant(
        X
    )

    y = sub[
        "overlap_binary"
    ]

    model = sm.GLM(
        y,
        X,
        family=sm.families.Binomial()
    ).fit()

    ci = model.conf_int()

    recurrence_beta = (
        model.params[
            "recurrent_binary"
        ]
    )

    width_beta = (
        model.params[
            "log10_width_z"
        ]
    )

    logistic_rows.append({
        "cell_line": cell,

        "n_regions": len(sub),

        "adjusted_or_recurrent":
            np.exp(recurrence_beta),

        "recurrent_ci_low":
            np.exp(
                ci.loc[
                    "recurrent_binary", 0
                ]
            ),

        "recurrent_ci_high":
            np.exp(
                ci.loc[
                    "recurrent_binary", 1
                ]
            ),

        "recurrent_raw_p":
            model.pvalues[
                "recurrent_binary"
            ],

        "width_or_per_1sd_log10":
            np.exp(width_beta),

        "width_ci_low":
            np.exp(
                ci.loc[
                    "log10_width_z", 0
                ]
            ),

        "width_ci_high":
            np.exp(
                ci.loc[
                    "log10_width_z", 1
                ]
            ),

        "width_raw_p":
            model.pvalues[
                "log10_width_z"
            ]
    })


logistic = pd.DataFrame(
    logistic_rows
)

logistic[
    "recurrent_holm_p"
] = multipletests(
    logistic["recurrent_raw_p"],
    method="holm"
)[1]

logistic[
    "width_holm_p"
] = multipletests(
    logistic["width_raw_p"],
    method="holm"
)[1]

logistic.to_csv(
    OUTPUT_LOGISTIC,
    index=False
)


# ============================================================
# Validation
# ============================================================

print("\nCROSS-CELL-LINE OVERLAP ANALYSIS")
print("=" * 80)

print(
    "Total regions:",
    len(distance_df)
)

print(
    "Regions overlapping another cell line:",
    int(
        distance_df[
            "overlap_other_line"
        ].sum()
    )
)

print(
    "Regions within 1 Mb of another cell line:",
    int(
        distance_df[
            "within_1mb_other_line"
        ].sum()
    )
)


print("\nEDGE DISTANCE / OVERLAP")
print("=" * 80)

for _, row in tests.iterrows():

    print(
        f"{row['cell_type']:<5} | "
        f"edge median single="
        f"{row['single_median_edge_distance_kb']:.2f} kb | "
        f"recurrent="
        f"{row['recurrent_median_edge_distance_kb']:.2f} kb | "
        f"distance Holm="
        f"{row['distance_holm_p']:.6f}"
    )

    print(
        f"      overlap single="
        f"{int(row['single_overlap_n'])} | "
        f"recurrent="
        f"{int(row['recurrent_overlap_n'])} | "
        f"OR="
        f"{row['overlap_odds_ratio']:.3f} | "
        f"Holm="
        f"{row['overlap_holm_p']:.6f}"
    )


print("\nWIDTH-ADJUSTED LOGISTIC REGRESSION")
print("=" * 80)

for _, row in logistic.iterrows():

    print(
        f"{row['cell_line']:<5} | "
        f"recurrence OR="
        f"{row['adjusted_or_recurrent']:.3f} "
        f"("
        f"{row['recurrent_ci_low']:.2f}-"
        f"{row['recurrent_ci_high']:.2f}"
        f") | "
        f"Holm p="
        f"{row['recurrent_holm_p']:.6f}"
    )


# Core structural validation
assert len(distance_df) == 644

assert (
    int(
        distance_df[
            "overlap_other_line"
        ].sum()
    )
    == 147
)

assert (
    int(
        distance_df[
            "within_1mb_other_line"
        ].sum()
    )
    == 287
)


# Validate key overlap counts
expected_overlap = {
    "U2OS": (18, 20),
    "RPE1": (34, 25),
    "BJ": (33, 17)
}

for cell, (
    expected_single,
    expected_recurrent
) in expected_overlap.items():

    row = tests.loc[
        tests["cell_type"] == cell
    ].iloc[0]

    assert int(
        row["single_overlap_n"]
    ) == expected_single

    assert int(
        row["recurrent_overlap_n"]
    ) == expected_recurrent


# Validate central adjusted OR results
expected_adjusted_or = {
    "U2OS": 6.1885,
    "RPE1": 1.6132,
    "BJ": 1.3867
}

for cell, expected_or in (
    expected_adjusted_or.items()
):

    observed = logistic.loc[
        logistic["cell_line"] == cell,
        "adjusted_or_recurrent"
    ].iloc[0]

    assert np.isclose(
        observed,
        expected_or,
        atol=0.02
    )


print("\nValidation: PASS")

print("\nSaved:")
print(OUTPUT_DISTANCE)
print(OUTPUT_TESTS)
print(OUTPUT_LOGISTIC)
