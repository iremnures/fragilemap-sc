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
    / "u2os_break_regions_multifeature.csv"
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

OUTPUT_TESTS = (
    TABLE_DIR
    / "u2os_recurrence_multifeature_tests.csv"
)

OUTPUT_BOOTSTRAP = (
    TABLE_DIR
    / "u2os_multifeature_effects_with_ci.csv"
)


# ============================================================
# Feature families
# ============================================================

PRIMARY_FEATURES = [
    "midas",
    "fancd2",
    "lateS_G2M",
    "g_quadruplex",
    "ns_seq",
    "gro_seq",
    "drip_seq"
]

SUPPORTING_FEATURES = [
    "ab_compartment",
    "rt_slope_untreated"
]

FEATURE_LABELS = {
    "midas": "MiDAS-seq",
    "fancd2": "FANCD2-seq",
    "lateS_G2M": "lateS/G2/M-seq",
    "g_quadruplex": "G-quadruplex",
    "ns_seq": "NS-seq",
    "gro_seq": "GRO-seq",
    "drip_seq": "DRIP-seq",
    "ab_compartment": "A/B compartment",
    "rt_slope_untreated": "RT slope"
}


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

assert len(df) == 195
assert df["recurrent"].isna().sum() == 0
assert int(df["recurrent"].sum()) == 44
assert int((~df["recurrent"]).sum()) == 151


# ============================================================
# Effect size
# ============================================================

def rank_biserial_from_u(
    u,
    n1,
    n2
):
    """
    Positive values mean group 1 tends to have
    higher values than group 2.
    """

    return (
        (2 * u)
        / (n1 * n2)
        - 1
    )


# ============================================================
# Mann-Whitney tests
# ============================================================

rows = []

all_features = (
    PRIMARY_FEATURES
    + SUPPORTING_FEATURES
)

for feature in all_features:

    recurrent = df.loc[
        df["recurrent"],
        feature
    ].dropna()

    single = df.loc[
        ~df["recurrent"],
        feature
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

    family = (
        "Primary"
        if feature in PRIMARY_FEATURES
        else "Supporting"
    )

    rows.append({
        "feature": feature,
        "label": FEATURE_LABELS[feature],
        "family": family,

        "single_n": len(single),
        "recurrent_n": len(recurrent),

        "single_median":
            single.median(),

        "recurrent_median":
            recurrent.median(),

        "mann_whitney_u":
            u,

        "rank_biserial_effect":
            effect,

        "raw_p":
            p
    })


results = pd.DataFrame(rows)


# ============================================================
# Holm correction WITHIN predefined feature families
# ============================================================

results["holm_p"] = np.nan

for family in [
    "Primary",
    "Supporting"
]:

    mask = (
        results["family"]
        == family
    )

    adjusted = multipletests(
        results.loc[
            mask,
            "raw_p"
        ],
        method="holm"
    )[1]

    results.loc[
        mask,
        "holm_p"
    ] = adjusted


results["significant_holm_0.05"] = (
    results["holm_p"] < 0.05
)

results.to_csv(
    OUTPUT_TESTS,
    index=False
)


# ============================================================
# Bootstrap rank-biserial 95% CI
#
# Resample recurrent and single groups independently.
# ============================================================

RNG = np.random.default_rng(20260924)

N_BOOT = 2000

bootstrap_rows = []


def rank_biserial_from_samples(
    group1,
    group2
):

    u, _ = mannwhitneyu(
        group1,
        group2,
        alternative="two-sided"
    )

    return rank_biserial_from_u(
        u,
        len(group1),
        len(group2)
    )


for feature in all_features:

    recurrent = (
        df.loc[
            df["recurrent"],
            feature
        ]
        .dropna()
        .to_numpy()
    )

    single = (
        df.loc[
            ~df["recurrent"],
            feature
        ]
        .dropna()
        .to_numpy()
    )

    observed = rank_biserial_from_samples(
        recurrent,
        single
    )

    boot_effects = np.empty(
        N_BOOT,
        dtype=float
    )

    for i in range(N_BOOT):

        recurrent_boot = RNG.choice(
            recurrent,
            size=len(recurrent),
            replace=True
        )

        single_boot = RNG.choice(
            single,
            size=len(single),
            replace=True
        )

        boot_effects[i] = (
            rank_biserial_from_samples(
                recurrent_boot,
                single_boot
            )
        )

    ci_low, ci_high = np.percentile(
        boot_effects,
        [2.5, 97.5]
    )

    row = results.loc[
        results["feature"]
        == feature
    ].iloc[0]

    bootstrap_rows.append({
        "feature": feature,
        "label": FEATURE_LABELS[feature],
        "family": row["family"],

        "effect": observed,

        "ci_low":
            ci_low,

        "ci_high":
            ci_high,

        "raw_p":
            row["raw_p"],

        "holm_p":
            row["holm_p"],

        "significant_holm_0.05":
            row["significant_holm_0.05"]
    })


bootstrap_df = pd.DataFrame(
    bootstrap_rows
)

bootstrap_df.to_csv(
    OUTPUT_BOOTSTRAP,
    index=False
)


# ============================================================
# Validation
# ============================================================

print("\nU2OS MULTI-FEATURE RECURRENCE ANALYSIS")
print("=" * 95)

for _, row in results.iterrows():

    print(
        f"{row['label']:<22} | "
        f"{row['family']:<10} | "
        f"single median="
        f"{row['single_median']:.4f} | "
        f"recurrent median="
        f"{row['recurrent_median']:.4f} | "
        f"RBC="
        f"{row['rank_biserial_effect']:+.4f} | "
        f"Holm p="
        f"{row['holm_p']:.6g}"
    )


print("\nBOOTSTRAP EFFECT-SIZE ANALYSIS")
print("=" * 95)

for _, row in bootstrap_df.iterrows():

    print(
        f"{row['label']:<22} "
        f"effect={row['effect']:>7.3f} "
        f"95% CI=("
        f"{row['ci_low']:>6.3f}, "
        f"{row['ci_high']:>6.3f}) "
        f"Holm p="
        f"{row['holm_p']:.6g}"
    )


# ============================================================
# Validate established central results
# ============================================================

expected_effects = {
    "midas": 0.4232,
    "fancd2": 0.3992,
    "lateS_G2M": 0.3110,
    "g_quadruplex": -0.1553,
    "ns_seq": -0.1321,
    "gro_seq": 0.0712,
    "drip_seq": 0.0024,
    "ab_compartment": -0.1650,
    "rt_slope_untreated": -0.0804
}

for feature, expected in (
    expected_effects.items()
):

    observed = results.loc[
        results["feature"]
        == feature,
        "rank_biserial_effect"
    ].iloc[0]

    assert np.isclose(
        observed,
        expected,
        atol=0.002
    ), (
        f"{feature}: effect mismatch "
        f"{observed} vs {expected}"
    )


# Expected Holm-significant primary features
significant_features = set(
    results.loc[
        results["significant_holm_0.05"],
        "feature"
    ]
)

expected_significant = {
    "midas",
    "fancd2",
    "lateS_G2M"
}

assert (
    significant_features
    == expected_significant
), (
    "Unexpected Holm-significant feature set: "
    f"{significant_features}"
)


# Validate central adjusted p-values
expected_holm = {
    "midas": 0.000138952,
    "fancd2": 0.000343719,
    "lateS_G2M": 0.00861381
}

for feature, expected in (
    expected_holm.items()
):

    observed = results.loc[
        results["feature"]
        == feature,
        "holm_p"
    ].iloc[0]

    assert np.isclose(
        observed,
        expected,
        atol=1e-5
    )


print("\nValidation: PASS")

print(
    "\nHolm-significant features:",
    ", ".join(
        FEATURE_LABELS[x]
        for x in [
            "midas",
            "fancd2",
            "lateS_G2M"
        ]
    )
)

print("\nSaved:")
print(OUTPUT_TESTS)
print(OUTPUT_BOOTSTRAP)
