from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

from sklearn.model_selection import (
    RepeatedStratifiedKFold,
    StratifiedKFold,
    cross_validate,
    cross_val_predict
)

from sklearn.metrics import (
    make_scorer,
    recall_score,
    roc_auc_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix
)


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

OUTPUT_CV = (
    TABLE_DIR
    / "u2os_ml_cross_validation.csv"
)

OUTPUT_SUMMARY = (
    TABLE_DIR
    / "u2os_ml_performance_summary.csv"
)

OUTPUT_OOF = (
    TABLE_DIR
    / "u2os_ml_oof_predictions.csv"
)

OUTPUT_COEF = (
    TABLE_DIR
    / "u2os_ml_coefficients.csv"
)


# ============================================================
# Data
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

y = df["recurrent"].astype(int)


# ============================================================
# Feature sets
# ============================================================

GENOMIC_FEATURES = [
    "midas",
    "fancd2",
    "lateS_G2M",
    "g_quadruplex",
    "ns_seq",
    "gro_seq",
    "drip_seq",
    "ab_compartment",
    "rt_slope_untreated"
]


# Model A
X_features = df[
    GENOMIC_FEATURES
].copy()


# Model B
X_width = X_features.copy()

X_width[
    "log10_region_width_kb"
] = np.log10(
    df["region_width_kb"]
)


# ============================================================
# Pipeline
# ============================================================

def make_model():

    return Pipeline([
        (
            "imputer",
            SimpleImputer(
                strategy="median"
            )
        ),
        (
            "scaler",
            StandardScaler()
        ),
        (
            "logistic",
            LogisticRegression(
                penalty="l2",
                C=1.0,
                class_weight="balanced",
                solver="liblinear",
                max_iter=5000,
                random_state=20260924
            )
        )
    ])


# ============================================================
# Repeated stratified CV
#
# 5 folds × 20 repeats = 100 test folds
# ============================================================

cv = RepeatedStratifiedKFold(
    n_splits=5,
    n_repeats=20,
    random_state=20260924
)


specificity_scorer = make_scorer(
    recall_score,
    pos_label=0
)


scoring = {
    "roc_auc": "roc_auc",
    "average_precision": "average_precision",
    "balanced_accuracy": "balanced_accuracy",
    "sensitivity": make_scorer(
        recall_score,
        pos_label=1
    ),
    "specificity": specificity_scorer
}


models = {
    "Genomic features only": X_features,
    "Genomic features + width": X_width
}


# ============================================================
# CV evaluation
# ============================================================

cv_rows = []

for model_name, X in models.items():

    scores = cross_validate(
        make_model(),
        X,
        y,
        cv=cv,
        scoring=scoring,
        return_train_score=False,
        n_jobs=1
    )

    n_folds = len(
        scores["test_roc_auc"]
    )

    for i in range(n_folds):

        cv_rows.append({
            "model": model_name,
            "fold": i + 1,

            "roc_auc":
                scores[
                    "test_roc_auc"
                ][i],

            "average_precision":
                scores[
                    "test_average_precision"
                ][i],

            "balanced_accuracy":
                scores[
                    "test_balanced_accuracy"
                ][i],

            "sensitivity":
                scores[
                    "test_sensitivity"
                ][i],

            "specificity":
                scores[
                    "test_specificity"
                ][i]
        })


cv_df = pd.DataFrame(
    cv_rows
)

cv_df.to_csv(
    OUTPUT_CV,
    index=False
)


# ============================================================
# Summary
# ============================================================

metrics = [
    "roc_auc",
    "average_precision",
    "balanced_accuracy",
    "sensitivity",
    "specificity"
]

summary_rows = []

for model_name in models:

    sub = cv_df[
        cv_df["model"]
        == model_name
    ]

    for metric in metrics:

        values = sub[metric]

        summary_rows.append({
            "model": model_name,
            "metric": metric,
            "mean": values.mean(),
            "std": values.std(),
            "median": values.median(),
            "q025": values.quantile(0.025),
            "q975": values.quantile(0.975)
        })


summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    OUTPUT_SUMMARY,
    index=False
)


# ============================================================
# One fixed 5-fold split for out-of-fold predictions
#
# Used only for interpretable confusion matrix / OOF prediction
# table. Repeated CV above is the main performance estimate.
# ============================================================

oof_cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=20260924
)


oof_frames = []

for model_name, X in models.items():

    probabilities = cross_val_predict(
        make_model(),
        X,
        y,
        cv=oof_cv,
        method="predict_proba",
        n_jobs=1
    )[:, 1]

    predictions = (
        probabilities >= 0.5
    ).astype(int)

    temp = pd.DataFrame({
        "model": model_name,
        "row_id": np.arange(
            len(df)
        ),
        "true_recurrent": y,
        "predicted_probability":
            probabilities,
        "predicted_class":
            predictions
    })

    oof_frames.append(temp)


oof = pd.concat(
    oof_frames,
    ignore_index=True
)

oof.to_csv(
    OUTPUT_OOF,
    index=False
)


# ============================================================
# Full-data standardized coefficients
#
# Descriptive only.
# Predictive performance comes from CV.
# ============================================================

coef_rows = []

for model_name, X in models.items():

    model = make_model()

    model.fit(
        X,
        y
    )

    coefs = (
        model.named_steps[
            "logistic"
        ].coef_[0]
    )

    for feature, coef in zip(
        X.columns,
        coefs
    ):

        coef_rows.append({
            "model": model_name,
            "feature": feature,
            "standardized_log_odds_coefficient":
                coef,
            "odds_ratio_per_1sd":
                np.exp(coef)
        })


coef_df = pd.DataFrame(
    coef_rows
)

coef_df.to_csv(
    OUTPUT_COEF,
    index=False
)


# ============================================================
# Console output
# ============================================================

prevalence = y.mean()

print("\nU2OS INTERPRETABLE ML")
print("=" * 85)

print(
    f"Samples: {len(y)} "
    f"(single={int((y == 0).sum())}, "
    f"recurrent={int((y == 1).sum())})"
)

print(
    f"Recurrent prevalence: "
    f"{prevalence:.3f}"
)

print(
    f"Random-ranking AP baseline: "
    f"{prevalence:.3f}"
)


print("\nREPEATED 5-FOLD CROSS-VALIDATION")
print("=" * 85)

for model_name in models:

    print(
        f"\n{model_name}"
    )

    sub = summary[
        summary["model"]
        == model_name
    ]

    for metric in metrics:

        row = sub[
            sub["metric"]
            == metric
        ].iloc[0]

        print(
            f"{metric:<20} "
            f"{row['mean']:.3f} "
            f"+/- {row['std']:.3f}"
        )


# ============================================================
# OOF diagnostics
# ============================================================

print("\nOUT-OF-FOLD DIAGNOSTICS")
print("=" * 85)

for model_name in models:

    sub = oof[
        oof["model"]
        == model_name
    ]

    y_true = sub[
        "true_recurrent"
    ].to_numpy()

    prob = sub[
        "predicted_probability"
    ].to_numpy()

    pred = sub[
        "predicted_class"
    ].to_numpy()

    roc = roc_auc_score(
        y_true,
        prob
    )

    ap = average_precision_score(
        y_true,
        prob
    )

    bal = balanced_accuracy_score(
        y_true,
        pred
    )

    tn, fp, fn, tp = (
        confusion_matrix(
            y_true,
            pred,
            labels=[0, 1]
        ).ravel()
    )

    sensitivity = (
        tp / (tp + fn)
    )

    specificity = (
        tn / (tn + fp)
    )

    print(
        f"\n{model_name}"
    )

    print(
        f"ROC-AUC:           {roc:.3f}"
    )

    print(
        f"Average precision: {ap:.3f}"
    )

    print(
        f"Balanced accuracy: {bal:.3f}"
    )

    print(
        f"Sensitivity:       {sensitivity:.3f}"
    )

    print(
        f"Specificity:       {specificity:.3f}"
    )

    print(
        "Confusion matrix: "
        f"TN={tn}, FP={fp}, "
        f"FN={fn}, TP={tp}"
    )


# ============================================================
# Basic validation
# ============================================================

assert (
    len(
        cv_df[
            cv_df["model"]
            == "Genomic features only"
        ]
    )
    == 100
)

assert (
    len(
        cv_df[
            cv_df["model"]
            == "Genomic features + width"
        ]
    )
    == 100
)

assert (
    set(coef_df["model"])
    == set(models.keys())
)

assert (
    summary[
        "mean"
    ].between(
        0,
        1
    ).all()
)


print("\nValidation: PASS")

print("\nSaved:")
print(OUTPUT_CV)
print(OUTPUT_SUMMARY)
print(OUTPUT_OOF)
print(OUTPUT_COEF)
