from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path.home() / "fragilemap_sc"

DATA_DIR = PROJECT_DIR / "data" / "processed"
TABLE_DIR = PROJECT_DIR / "results" / "tables"

OUTPUT_FILE = TABLE_DIR / "final_qc_report.csv"

checks = []


def add_check(name, passed, observed, expected):
    checks.append({
        "check": name,
        "status": "PASS" if passed else "FAIL",
        "observed": str(observed),
        "expected": str(expected)
    })


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


# ============================================================
# 1. Main break-region dataset
# ============================================================

main = pd.read_csv(
    DATA_DIR / "human_break_regions_enriched.csv"
)

main["recurrent_bool"] = parse_bool(
    main["recurrent"]
)

add_check(
    "Total human break regions",
    len(main) == 644,
    len(main),
    644
)

expected_counts = {
    "U2OS": 195,
    "RPE1": 284,
    "BJ": 165
}

for cell, expected in expected_counts.items():

    observed = int(
        (main["cell_type"] == cell).sum()
    )

    add_check(
        f"{cell} unique break regions",
        observed == expected,
        observed,
        expected
    )


expected_events = {
    "U2OS": 274,
    "RPE1": 434,
    "BJ": 244
}

for cell, expected in expected_events.items():

    observed = int(
        main.loc[
            main["cell_type"] == cell,
            "frequency"
        ].sum()
    )

    add_check(
        f"{cell} total break events",
        observed == expected,
        observed,
        expected
    )


duplicates = int(
    main.duplicated(
        subset=[
            "cell_type",
            "chrom",
            "start",
            "end"
        ]
    ).sum()
)

add_check(
    "Duplicate genomic intervals",
    duplicates == 0,
    duplicates,
    0
)


invalid_coordinates = int(
    (main["start"] >= main["end"]).sum()
)

add_check(
    "Invalid coordinates",
    invalid_coordinates == 0,
    invalid_coordinates,
    0
)


recurrence_mismatch = int(
    (
        (main["frequency"] > 1)
        != main["recurrent_bool"]
    ).sum()
)

add_check(
    "Recurrent flag consistency",
    recurrence_mismatch == 0,
    recurrence_mismatch,
    0
)


width_mismatch = int(
    (
        np.abs(
            (main["end"] - main["start"])
            - main["region_width_bp"]
        ) > 0
    ).sum()
)

add_check(
    "Region width consistency",
    width_mismatch == 0,
    width_mismatch,
    0
)


# ============================================================
# 2. RT integration
# ============================================================

rt = pd.read_csv(
    DATA_DIR
    / "human_break_regions_rt_classified.csv"
)

add_check(
    "RT dataset row count",
    len(rt) == 644,
    len(rt),
    644
)

missing_rt = int(
    rt[
        "mean_pseudobulk_rt_untreated"
    ].isna().sum()
)

add_check(
    "Missing RT values",
    missing_rt == 0,
    missing_rt,
    0
)

observed_classes = set(
    rt["rt_class"].dropna()
)

expected_classes = {
    "Early",
    "Mid",
    "Late"
}

add_check(
    "RT classes",
    observed_classes == expected_classes,
    sorted(observed_classes),
    sorted(expected_classes)
)


# ============================================================
# 3. Cross-line overlap
# ============================================================

overlap = pd.read_csv(
    DATA_DIR
    / "human_break_regions_overlap_distance.csv"
)

add_check(
    "Cross-line dataset row count",
    len(overlap) == 644,
    len(overlap),
    644
)

overlap["overlap_bool"] = parse_bool(
    overlap["overlap_other_line"]
)

negative_distance = int(
    (
        overlap["nearest_edge_distance_kb"]
        < 0
    ).sum()
)

add_check(
    "Negative edge distances",
    negative_distance == 0,
    negative_distance,
    0
)

overlap_mismatch = int(
    (
        overlap["overlap_bool"]
        &
        (
            overlap[
                "nearest_edge_distance_bp"
            ] != 0
        )
    ).sum()
)

add_check(
    "Overlap implies edge distance = 0",
    overlap_mismatch == 0,
    overlap_mismatch,
    0
)

add_check(
    "Total cross-line overlaps",
    int(overlap["overlap_bool"].sum()) == 147,
    int(overlap["overlap_bool"].sum()),
    147
)


# ============================================================
# 4. U2OS multifeature dataset
# ============================================================

u2os = pd.read_csv(
    DATA_DIR
    / "u2os_break_regions_multifeature.csv"
)

add_check(
    "U2OS multifeature row count",
    len(u2os) == 195,
    len(u2os),
    195
)

expected_feature_coverage = {
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

for feature, expected in (
    expected_feature_coverage.items()
):

    observed = int(
        u2os[feature].notna().sum()
    )

    add_check(
        f"{feature} coverage",
        observed == expected,
        observed,
        expected
    )


# ============================================================
# 5. Key result-table checks
# ============================================================

feature_results = pd.read_csv(
    TABLE_DIR
    / "u2os_recurrence_multifeature_tests.csv"
)

sig_features = set(
    feature_results.loc[
        feature_results["holm_p"] < 0.05,
        "feature"
    ]
)

expected_sig = {
    "midas",
    "fancd2",
    "lateS_G2M"
}

add_check(
    "U2OS significant feature set",
    sig_features == expected_sig,
    sorted(sig_features),
    sorted(expected_sig)
)


logistic = pd.read_csv(
    TABLE_DIR
    / "width_adjusted_overlap_logistic.csv"
)

u2os_or = float(
    logistic.loc[
        logistic["cell_line"] == "U2OS",
        "adjusted_or_recurrent"
    ].iloc[0]
)

add_check(
    "U2OS width-adjusted overlap OR",
    np.isclose(
        u2os_or,
        6.1885,
        atol=0.02
    ),
    round(u2os_or, 4),
    "approximately 6.19"
)


# ============================================================
# Final report
# ============================================================

qc = pd.DataFrame(checks)

qc.to_csv(
    OUTPUT_FILE,
    index=False
)

print("\nFINAL PIPELINE QC")
print("=" * 90)

for _, row in qc.iterrows():

    symbol = (
        "✓"
        if row["status"] == "PASS"
        else "✗"
    )

    print(
        f"{symbol} {row['status']:<4} | "
        f"{row['check']:<45} | "
        f"observed={row['observed']}"
    )


n_pass = int(
    (qc["status"] == "PASS").sum()
)

n_fail = int(
    (qc["status"] == "FAIL").sum()
)

print("\n" + "=" * 90)

print("PASS:", n_pass)
print("FAIL:", n_fail)

if n_fail == 0:
    print(
        "\nQC STATUS: ALL PIPELINE CHECKS PASSED"
    )
else:
    raise RuntimeError(
        f"Final QC failed: {n_fail} check(s)"
    )

print("\nSaved:", OUTPUT_FILE)
