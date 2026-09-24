from pathlib import Path
import pandas as pd


PROJECT_DIR = Path.home() / "fragilemap_sc"

TABLE_DIR = PROJECT_DIR / "results" / "tables"
DOCS_DIR = PROJECT_DIR / "docs"

DOCS_DIR.mkdir(parents=True, exist_ok=True)

MASTER_FILE = (
    TABLE_DIR
    / "master_results_summary_final.csv"
)

ML_FILE = (
    TABLE_DIR
    / "u2os_ml_performance_summary.csv"
)

OUTPUT_CSV = (
    TABLE_DIR
    / "master_results_summary_final.csv"
)

OUTPUT_MD = (
    DOCS_DIR
    / "master_results_summary_final.md"
)


# ============================================================
# Load existing statistical summary
# ============================================================

master = pd.read_csv(MASTER_FILE)

# Remove old predictive-modeling rows if script is rerun
master = master[
    master["analysis"]
    != "Predictive modeling"
].copy()


# ============================================================
# Load repeated-CV ML summary
# ============================================================

ml = pd.read_csv(ML_FILE)

models = [
    "Genomic features only",
    "Genomic features + width"
]


def get_mean(model, metric):

    return float(
        ml.loc[
            (ml["model"] == model)
            & (ml["metric"] == metric),
            "mean"
        ].iloc[0]
    )


ml_rows = []

for model in models:

    roc = get_mean(
        model,
        "roc_auc"
    )

    ap = get_mean(
        model,
        "average_precision"
    )

    bal = get_mean(
        model,
        "balanced_accuracy"
    )

    if model == "Genomic features only":

        interpretation = (
            "Genomic features alone provide moderate "
            "discrimination of recurrent U2OS break regions."
        )

        direction = (
            "Predictive signal from genomic features"
        )

    else:

        interpretation = (
            "Adding region width substantially improves "
            "classification performance; this model is treated "
            "as a sensitivity analysis because region width may "
            "partly reflect break-region definition and merging."
        )

        direction = (
            "Improved performance after adding region width"
        )

    ml_rows.append({
        "analysis":
            "Predictive modeling",

        "comparison":
            model,

        "test":
            (
                "L2-regularized logistic regression; "
                "20x repeated stratified 5-fold CV"
            ),

        "effect_or_statistic":
            (
                f"ROC-AUC={roc:.3f}; "
                f"AP={ap:.3f}; "
                f"balanced accuracy={bal:.3f}"
            ),

        "p_value":
            "NA",

        "p_value_type":
            "Not applicable",

        "direction":
            direction,

        "interpretation":
            interpretation
    })


master = pd.concat(
    [
        master,
        pd.DataFrame(ml_rows)
    ],
    ignore_index=True
)


# ============================================================
# Save CSV
# ============================================================

master.to_csv(
    OUTPUT_CSV,
    index=False
)


# ============================================================
# Save Markdown
# ============================================================

with open(OUTPUT_MD, "w") as f:

    f.write(
        "# Master Results Summary\n\n"
    )

    f.write(
        "Major statistical and predictive-modeling findings "
        "from the FragileMap-SC analysis.\n\n"
    )

    f.write(
        "| Analysis | Comparison | Effect / statistic | "
        "p-value | p-value type | Interpretation |\n"
    )

    f.write(
        "|---|---|---|---:|---|---|\n"
    )

    for _, row in master.iterrows():

        f.write(
            f"| {row['analysis']} "
            f"| {row['comparison']} "
            f"| {row['effect_or_statistic']} "
            f"| {row['p_value']} "
            f"| {row['p_value_type']} "
            f"| {row['interpretation']} |\n"
        )


# ============================================================
# Validation
# ============================================================

genomic_roc = get_mean(
    "Genomic features only",
    "roc_auc"
)

genomic_ap = get_mean(
    "Genomic features only",
    "average_precision"
)

width_roc = get_mean(
    "Genomic features + width",
    "roc_auc"
)

width_ap = get_mean(
    "Genomic features + width",
    "average_precision"
)

assert abs(
    genomic_roc - 0.696
) < 0.01

assert abs(
    genomic_ap - 0.530
) < 0.01

assert abs(
    width_roc - 0.870
) < 0.01

assert abs(
    width_ap - 0.720
) < 0.01

assert (
    master["analysis"]
    == "Predictive modeling"
).sum() == 2


print("\nFINAL PROJECT SUMMARY")
print("=" * 80)

print(
    master.tail(2)[
        [
            "comparison",
            "effect_or_statistic",
            "interpretation"
        ]
    ].to_string(index=False)
)

print("\nValidation: PASS")

print("\nSaved:")
print(OUTPUT_CSV)
print(OUTPUT_MD)
