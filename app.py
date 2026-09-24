from pathlib import Path
import json

import joblib
import numpy as np
import pandas as pd
import streamlit as st


# ------------------------------------------------------------
# Page
# ------------------------------------------------------------

st.set_page_config(
    page_title="U2OS Break Recurrence Explorer",
    page_icon="🧬",
    layout="wide"
)

st.markdown("""
<style>
.block-container {
    max-width: 1100px;
    padding-top: 2.2rem;
    padding-bottom: 3rem;
}

h1, h2, h3 {
    color: #243B53;
}

h1 {
    font-weight: 700;
    letter-spacing: -0.02em;
}

div[data-testid="stMetric"] {
    background: #FAFAFA;
    border: 1px solid #E5E7EB;
    padding: 18px;
    border-radius: 10px;
}

div.stButton > button {
    background-color: #8F2948;
    color: white;
    border: none;
    border-radius: 7px;
    padding: 0.65rem 1.3rem;
    font-weight: 600;
}

div.stButton > button:hover {
    background-color: #74223C;
    color: white;
}

.note-box {
    background: #F8FAFC;
    color: #334155;
    border-left: 4px solid #8F2948;
    padding: 14px 16px;
    border-radius: 6px;
    margin: 1rem 0 1.6rem 0;
}

.small-muted {
    color: #667085;
    font-size: 0.92rem;
}
</style>
""", unsafe_allow_html=True)


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

PROJECT_DIR = Path(__file__).resolve().parent
MODEL_DIR = PROJECT_DIR / "models"


# ------------------------------------------------------------
# Load models
# ------------------------------------------------------------

@st.cache_resource
def load_models():
    genomic = joblib.load(
        MODEL_DIR / "u2os_genomic_features_model.joblib"
    )

    genomic_width = joblib.load(
        MODEL_DIR / "u2os_genomic_features_width_model.joblib"
    )

    return genomic, genomic_width


@st.cache_data
def load_metadata():
    with open(MODEL_DIR / "app_metadata.json") as f:
        return json.load(f)


try:
    genomic_model, width_model = load_models()
    metadata = load_metadata()

except FileNotFoundError:
    st.error(
        "Model files are missing. Run "
        "`python scripts/12_fit_final_ml_models.py` first."
    )
    st.stop()


# ------------------------------------------------------------
# Header
# ------------------------------------------------------------

st.title("U2OS Break Recurrence Explorer")

st.markdown(
    """
    Explore how replication-stress and genomic features relate to
    recurrent chromosome break regions in **U2OS cells**.
    """
)

st.markdown(
    """
    <div class="note-box">
    This is an exploratory research model trained on 195 U2OS break regions.
    The output is a model score, not a calibrated biological or clinical probability.
    </div>
    """,
    unsafe_allow_html=True
)


# ------------------------------------------------------------
# Model selection
# ------------------------------------------------------------

st.subheader("Model")

model_choice = st.radio(
    "Feature set",
    [
        "Genomic features only",
        "Genomic features + region width"
    ],
    horizontal=True
)

if model_choice == "Genomic features only":
    model_bundle = genomic_model
else:
    model_bundle = width_model


# ------------------------------------------------------------
# Feature inputs
# ------------------------------------------------------------

st.subheader("Region features")

st.markdown(
    """
    <div class="small-muted">
    Inputs start at the median values observed in the U2OS dataset.
    Replace them with values from a region you want to explore.
    </div>
    """,
    unsafe_allow_html=True
)

feature_labels = {
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

features = [
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

inputs = {}

cols = st.columns(3)

for i, feature in enumerate(features):
    with cols[i % 3]:
        inputs[feature] = st.number_input(
            feature_labels[feature],
            value=float(metadata[feature]["median"]),
            step=0.05,
            format="%.4f"
        )


region_width = None

if model_choice == "Genomic features + region width":

    st.markdown("---")

    region_width = st.number_input(
        "Region width (kb)",
        min_value=1.0,
        value=float(metadata["region_width_kb"]["median"]),
        step=25.0,
        format="%.2f"
    )


# ------------------------------------------------------------
# Prediction
# ------------------------------------------------------------

st.markdown("")

if st.button("Estimate recurrence score"):

    row = {
        feature: inputs[feature]
        for feature in features
    }

    if model_choice == "Genomic features + region width":
        row["log10_region_width_kb"] = np.log10(region_width)

    X = pd.DataFrame([row])
    X = X[model_bundle["features"]]

    pipeline = model_bundle["pipeline"]

    score = float(
        pipeline.predict_proba(X)[0, 1]
    )

    pattern = (
        "Recurrent-like"
        if score >= 0.5
        else "Single-occurrence-like"
    )

    st.markdown("---")
    st.subheader("Result")

    left, right = st.columns([1, 2])

    with left:
        st.metric(
            "Recurrence score",
            f"{score:.1%}"
        )

        st.markdown(
            f"**Pattern:** {pattern}"
        )

    with right:
        st.progress(score)

        st.markdown(
            """
            <div class="small-muted">
            Higher values indicate that the supplied feature profile
            more closely resembles recurrent regions in the fitted U2OS dataset.
            </div>
            """,
            unsafe_allow_html=True
        )


    # --------------------------------------------------------
    # Feature contributions
    # --------------------------------------------------------

    st.subheader("Feature contributions")

    imputer = pipeline.named_steps["imputer"]
    scaler = pipeline.named_steps["scaler"]
    logistic = pipeline.named_steps["logistic"]

    imputed = imputer.transform(X)
    standardized = scaler.transform(imputed)[0]

    coefficients = logistic.coef_[0]
    contributions = standardized * coefficients

    pretty_names = {
        **feature_labels,
        "log10_region_width_kb": "Region width"
    }

    contribution_df = pd.DataFrame({
        "Feature": [
            pretty_names[x]
            for x in model_bundle["features"]
        ],
        "Contribution": contributions
    })

    contribution_df["Direction"] = np.where(
        contribution_df["Contribution"] >= 0,
        "Raises score",
        "Lowers score"
    )

    contribution_df["absolute"] = (
        contribution_df["Contribution"].abs()
    )

    contribution_df = (
        contribution_df
        .sort_values(
            "absolute",
            ascending=False
        )
        .drop(columns="absolute")
    )

    contribution_df["Contribution"] = (
        contribution_df["Contribution"].round(3)
    )

    st.dataframe(
        contribution_df,
        hide_index=True,
        use_container_width=True
    )


# ------------------------------------------------------------
# Model information
# ------------------------------------------------------------

st.markdown("---")

st.subheader("About the model")

st.markdown(
    """
    This explorer uses **L2-regularized logistic regression**.

    - **Genomic features only:** nine genomic and replication-stress features.
    - **Genomic features + region width:** the same feature set plus break-region width.

    Predictive performance was evaluated separately using
    **20 repeats of stratified 5-fold cross-validation**.

    The models used by this interface are fitted to all 195 U2OS regions
    and are intended for exploratory use.
    """
)

st.caption(
    "FragileMap-SC · exploratory genomic analysis"
)
