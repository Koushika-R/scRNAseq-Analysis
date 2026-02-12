#!/bin/bash

# === Configuration ===
# Mapping of GSM IDs to their sample codes
declare -A samples=(
    [GSM8561110]="S_1907_004567"
    [GSM8561111]="S_1907_004614"
    [GSM8561112]="S_1907_004802"
    [GSM8561113]="S_1907_004896"
    [GSM8561114]="S_1907_005225"
    [GSM8561115]="S_2006_004078"
    [GSM8561116]="S_2007_002809"
    [GSM8561117]="S_2007_002950"
    [GSM8561118]="S_2007_002997"
    [GSM8561119]="S_2007_003044"
    [GSM8561120]="S_2007_003091"
    [GSM8561121]="S_2007_003138"
    [GSM8561122]="S_2007_003793"
    [GSM8561123]="S_2007_003840"
    [GSM8561124]="S_2007_003934"
    [GSM8561125]="S_2007_004028"
    [GSM8561126]="S_2007_004216"
    [GSM8561127]="S_2007_004263"
    [GSM8561128]="S_2103_004019"
    [GSM8561129]="S_2103_004028"
    [GSM8561130]="S_2103_004037"
    [GSM8561131]="S_2103_004046"
    [GSM8561132]="S_2103_004064"
    [GSM8561133]="S_2103_004073"
    [GSM8561134]="S_2103_004091"
    [GSM8561135]="S_2103_004100"
    [GSM8561136]="S_2103_004109"
    [GSM8561137]="S_2103_004118"
    [GSM8561138]="S_2103_004127"
    [GSM8561139]="S_2103_004145"
    [GSM8561140]="S_2103_004154"
    [GSM8561141]="S_2103_004163"
    [GSM8561142]="S_2103_004181"
    [GSM8561143]="S_2107_023520"
    [GSM8561144]="S_2107_023538"
    [GSM8561145]="S_2107_023547"
    [GSM8561146]="S_2107_023556"
    [GSM8561147]="S_2107_023592"
    [GSM8561148]="S_2107_023610"
    [GSM8561149]="S_2107_023619"
)

# Base URL for GEO sample supplemental files
BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples"

# Files to download per GSM
FILES=(
    "barcodes.tsv.gz"
    "barcodes_processed.tsv.gz"
    "features.tsv.gz"
    "features_processed.tsv.gz"
    "matrix.mtx.gz"
    "matrix_processed.mtx.gz"
)

# === Download Process ===
# Create main output directory
BASE_DIR="~GSE279086_Practice/raw files"
mkdir -p "$BASE_DIR"
cd "$BASE_DIR" || exit 1

for GSM in "${!samples[@]}"; do
    SAMPLE="${samples[$GSM]}"
    SUBDIR="${GSM:0:7}nnn"
    DEST_DIR="${GSM}"
    mkdir -p "$DEST_DIR"

    echo " Downloading for $GSM ($SAMPLE)..."
    for FILE in "${FILES[@]}"; do
        FILENAME="${GSM}_${SAMPLE}_${FILE}"
        URL="${BASE_URL}/${SUBDIR}/${GSM}/suppl/${FILENAME}"

        if [[ ! -f "${DEST_DIR}/${FILENAME}" ]]; then
            echo "  → $FILENAME"
            # <-- use curl instead of wget -->
            curl -L -o "${DEST_DIR}/${FILENAME}" "$URL" || echo " Failed to download $FILENAME"
        else
            echo " Already exists: $FILENAME"
        fi
    done
    echo ""
done

echo " All downloads completed successfully!"
