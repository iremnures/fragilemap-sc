import csv
import openpyxl
from pathlib import Path

# ============================================================
# Paths
# ============================================================

PROJECT_DIR = Path.home() / "fragilemap_sc"

INPUT_FILE = (
    PROJECT_DIR
    / "data"
    / "raw"
    / "41467_2026_76451_MOESM3_ESM.xlsx"
)

OUTPUT_FILE = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "human_break_regions_enriched.csv"
)

# ============================================================
# Cell-line mapping
# ============================================================

SHEET_MAP = {
    "U2OS": "U2OS",
    "hTERT-RPE1": "RPE1",
    "BJ-hTERT": "BJ"
}

# ============================================================
# Load supplementary workbook
# ============================================================

wb = openpyxl.load_workbook(
    INPUT_FILE,
    read_only=True,
    data_only=True
)

rows_written = 0
summary = {}

# ============================================================
# Parse break regions
# ============================================================

with open(OUTPUT_FILE, "w", newline="") as outfile:

    writer = csv.writer(outfile)

    writer.writerow([
        "cell_type",
        "chrom",
        "start",
        "end",
        "frequency",
        "recurrent",
        "region_width_bp",
        "region_width_kb"
    ])

    for sheet_name, cell_type in SHEET_MAP.items():

        ws = wb[sheet_name]

        n_regions = 0
        total_events = 0
        recurrent_regions = 0

        # Rows 1–3 contain supplementary-table metadata.
        # Actual data begin at row 4.
        for row in ws.iter_rows(
            min_row=4,
            values_only=True
        ):

            chrom, start, end, frequency = row[:4]

            if chrom is None:
                continue

            start = int(start)
            end = int(end)
            frequency = int(frequency)

            recurrent = frequency > 1

            # Interval width.
            # No +1 is added because coordinate convention is
            # retained exactly as represented in the source intervals.
            width_bp = end - start
            width_kb = width_bp / 1000

            writer.writerow([
                cell_type,
                chrom,
                start,
                end,
                frequency,
                recurrent,
                width_bp,
                round(width_kb, 2)
            ])

            n_regions += 1
            total_events += frequency

            if recurrent:
                recurrent_regions += 1

            rows_written += 1

        summary[cell_type] = {
            "regions": n_regions,
            "events": total_events,
            "recurrent": recurrent_regions
        }

wb.close()

# ============================================================
# Validation
# ============================================================

expected = {
    "U2OS": {
        "regions": 195,
        "events": 274
    },
    "RPE1": {
        "regions": 284,
        "events": 434
    },
    "BJ": {
        "regions": 165,
        "events": 244
    }
}

print("\nBREAK-REGION DATA PREPARATION")
print("=" * 60)

for cell_type in ["U2OS", "RPE1", "BJ"]:

    observed = summary[cell_type]

    print(f"\n{cell_type}")
    print("Unique break regions:", observed["regions"])
    print("Total break events:", observed["events"])
    print("Recurrent regions:", observed["recurrent"])

    assert (
        observed["regions"]
        == expected[cell_type]["regions"]
    ), f"{cell_type}: unexpected region count"

    assert (
        observed["events"]
        == expected[cell_type]["events"]
    ), f"{cell_type}: unexpected event count"

assert rows_written == 644, "Unexpected total row count"

print("\nValidation: PASS")
print("Total human break regions:", rows_written)
print("Saved:", OUTPUT_FILE)
