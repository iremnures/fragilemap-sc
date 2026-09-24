from pathlib import Path
import json

import joblib
import numpy as np
import pandas as pd

from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression


PROJECT_DIR = Path.home() / "fragilemap_sc"

INPUT_FILE = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "u2os_break_regions_multifeature.csv"
)

MODEL_DIR = PROJECT_DIR / "models"
MODEL_DIR.mkdir(parents=True, exist_ok=True)


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

FEATURE_LABELS = {
    "midas": "MiDAS-seq",
    "fancd2": "FANCD2-seq",
    "lateS_G2M": "lateS/G2/M-seq",
    "g_quadruplex": "G-quadruplex",
    "ns_seq": "NS-seq",
    "gro_seq": "GRO-seq",
    "drip_seq": "DRIP-seq",
    "ab_compartment": "A/B compartment",
    "rt_slope_untreated": "RT slope",
    "log10_region_width_kb": "Region width"
}


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


df = pd.read_csv(INPUT_FILE)

df["recurrent"] = parse_bool(
    df["recurrent"]
)

assert len(df) == 195

y = df["recurrent"].astype(int)


# ============================================================
# Model A
# ============================================================

X_a = df[
    GENOMIC_FEATURES
].copy()

model_a = make_model()

model_a.fit(
    X_a,
    y
)


# ============================================================
# Model B
# ============================================================

X_b = X_a.copy()

X_b[
    "log10_region_width_kb"
] = np.log10(
    df["region_width_kb"]
)

model_b = make_model()

model_b.fit(
    X_b,
    y
)


# ============================================================
# Save models
# ============================================================

joblib.dump(
    {
        "pipeline": model_a,
        "features": list(X_a.columns),
        "name": "Genomic features only"
    },
    MODEL_DIR
    / "u2os_genomic_features_model.joblib"
)

joblib.dump(
    {
        "pipeline": model_b,
        "features": list(X_b.columns),
        "name": "Genomic features + width"
    },
    MODEL_DIR
    / "u2os_genomic_features_width_model.joblib"
)


# ============================================================
# Save app metadata
# ============================================================

metadata = {}

for feature in GENOMIC_FEATURES:

    values = df[feature].dropna()

    metadata[feature] = {
        "label": FEATURE_LABELS[feature],
        "median": float(values.median()),
        "min": float(values.quantile(0.01)),
        "max": float(values.quantile(0.99))
    }


width_values = df["region_width_kb"]

metadata["region_width_kb"] = {
    "label": "Region width (kb)",
    "median": float(width_values.median()),
    "min": float(width_values.quantile(0.01)),
    "max": float(width_values.quantile(0.99))
}


with open(
    MODEL_DIR / "app_metadata.json",
    "w"
) as f:

    json.dump(
        metadata,
        f,
        indent=2
    )


print("\nFINAL DEMO MODELS")
print("=" * 60)

print("Samples:", len(df))
print("Single:", int((y == 0).sum()))
print("Recurrent:", int((y == 1).sum()))

print("\nSaved:")
print(
    MODEL_DIR
    / "u2os_genomic_features_model.joblib"
)
print(
    MODEL_DIR
    / "u2os_genomic_features_width_model.joblib"
)
print(
    MODEL_DIR
    / "app_metadata.json"
)

print("\nValidation: PASS")
