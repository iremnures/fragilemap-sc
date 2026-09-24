from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = Path.home() / "fragilemap_sc"

INPUT_FILE = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "u2os_ml_cross_validation.csv"
)

OUTPUT_PNG = (
    PROJECT_DIR
    / "results"
    / "figures"
    / "final"
    / "main"
    / "figure_07_u2os_ml_performance.png"
)

OUTPUT_SVG = (
    PROJECT_DIR
    / "results"
    / "figures"
    / "final"
    / "main"
    / "figure_07_u2os_ml_performance.svg"
)


# ============================================================
# Load data
# ============================================================

df = pd.read_csv(INPUT_FILE)

metrics = [
    ("roc_auc", "ROC-AUC"),
    ("average_precision", "Average precision"),
    ("balanced_accuracy", "Balanced accuracy")
]

models = [
    "Genomic features only",
    "Genomic features + width"
]


# ============================================================
# Mean and SD across repeated CV folds
# ============================================================

summary = (
    df
    .groupby("model")
    [[metric for metric, _ in metrics]]
    .agg(["mean", "std"])
)


# ============================================================
# Style
# ============================================================

NAVY = "#243B53"
BURGUNDY = "#972E4F"
LIGHT_GRAY = "#C8C8C8"

fig, ax = plt.subplots(
    figsize=(11.2, 7.6)
)

x = np.arange(len(metrics))

bar_width = 0.28
offset = 0.16

positions = {
    "Genomic features only":
        x - offset,

    "Genomic features + width":
        x + offset
}

colors = {
    "Genomic features only":
        LIGHT_GRAY,

    "Genomic features + width":
        BURGUNDY
}


# ============================================================
# Bars
# ============================================================

for model in models:

    means = np.array([
        summary.loc[
            model,
            (metric, "mean")
        ]
        for metric, _ in metrics
    ])

    stds = np.array([
        summary.loc[
            model,
            (metric, "std")
        ]
        for metric, _ in metrics
    ])

    xpos = positions[model]

    # Bars
    ax.bar(
        xpos,
        means,
        width=bar_width,
        color=colors[model],
        edgecolor="#222222",
        linewidth=0.8,
        zorder=2
    )

    # Error bars
    ax.errorbar(
        xpos,
        means,
        yerr=stds,
        fmt="none",
        ecolor="#111111",
        elinewidth=1.8,
        capsize=5,
        capthick=1.8,
        zorder=3
    )

    # Value labels ABOVE error bars
    for xp, mean, sd in zip(
        xpos,
        means,
        stds
    ):

        ax.text(
            xp,
            mean + sd + 0.025,
            f"{mean:.2f}",
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="bold",
            color="#111111",
            zorder=5
        )


# ============================================================
# Axes
# ============================================================

ax.set_xticks(x)

ax.set_xticklabels(
    [label for _, label in metrics],
    fontsize=12
)

ax.set_ylabel(
    "Repeated cross-validation performance",
    fontsize=13
)

ax.set_ylim(
    0,
    1.05
)

ax.set_yticks(
    np.arange(
        0,
        1.01,
        0.2
    )
)

ax.tick_params(
    axis="y",
    labelsize=11
)

ax.tick_params(
    axis="x",
    pad=8
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.spines["left"].set_linewidth(1.0)
ax.spines["bottom"].set_linewidth(1.0)

ax.grid(False)


# ============================================================
# Title
# ============================================================

fig.suptitle(
    "Region width improves recurrence classification in U2OS",
    x=0.5,
    y=0.955,
    ha="center",
    fontsize=19,
    fontweight="bold",
    color=NAVY
)


# ============================================================
# Legend
# ============================================================

legend_handles = [
    Patch(
        facecolor=LIGHT_GRAY,
        edgecolor="#222222",
        label="Genomic features only"
    ),

    Patch(
        facecolor=BURGUNDY,
        edgecolor="#222222",
        label="Genomic features + width"
    )
]

fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.905),
    ncol=2,
    frameon=False,
    fontsize=11.5,
    handlelength=1.8,
    columnspacing=2.0
)


# ============================================================
# Footnote
# ============================================================

fig.text(
    0.5,
    0.045,
    (
        "Mean ± SD across 20 repeats of stratified "
        "5-fold cross-validation (100 test folds)."
    ),
    ha="center",
    va="center",
    fontsize=10,
    color="#333333"
)


# ============================================================
# Layout
# ============================================================

plt.subplots_adjust(
    left=0.11,
    right=0.97,
    top=0.79,
    bottom=0.16
)


# ============================================================
# Save
# ============================================================

plt.savefig(
    OUTPUT_PNG,
    dpi=300
)

plt.savefig(
    OUTPUT_SVG
)

plt.close()


print("Saved:")
print(OUTPUT_PNG)
print(OUTPUT_SVG)
