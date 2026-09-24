#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$HOME/fragilemap_sc"
SCRIPT_DIR="$PROJECT_DIR/scripts"

echo
echo "============================================================"
echo " FragileMap-SC reproducible analysis pipeline"
echo "============================================================"
echo

run_step () {
    echo
    echo "------------------------------------------------------------"
    echo "Running: $1"
    echo "------------------------------------------------------------"
    python3 "$SCRIPT_DIR/$1"
}

run_step "01_prepare_break_regions.py"
run_step "02_integrate_replication_timing.py"
run_step "03_replication_timing_analysis.py"
run_step "04_crossline_overlap_analysis.py"
run_step "05_break_region_width_analysis.py"
run_step "06_chromosome_landscape.py"
run_step "07_integrate_u2os_multifeature.py"
run_step "08_u2os_multifeature_analysis.py"

# Core QC before exploratory ML
run_step "09_final_qc.py"

# Interpretable ML
run_step "10_u2os_logistic_ml.py"

# Final ML figure
run_step "11_plot_ml_performance.py"

# Fit final exploratory models used by the Streamlit app
run_step "12_fit_final_ml_models.py"

# Final project-level results summary
run_step "13_finalize_project_summary.py"

echo
echo "============================================================"
echo " PIPELINE COMPLETED SUCCESSFULLY"
echo "============================================================"
echo
echo "Processed data:"
echo "  $PROJECT_DIR/data/processed/"
echo
echo "Statistical tables:"
echo "  $PROJECT_DIR/results/tables/"
echo
echo "Final figures:"
echo "  $PROJECT_DIR/results/figures/final/"
echo
