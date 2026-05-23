#!/usr/bin/env bash
# Batch run AGAT statistics, AGAT comparison, and gffcompare.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS="${ANALYSIS_DIR:-$SCRIPT_DIR}"
PYTHON_BIN="${PYTHON:-python}"

usage() {
    cat <<'EOF'
Usage: bash run_analyses.sh [--analysis-dir DIR] [SPECIES_ID ...]

Runs AGAT statistics, AGAT pair comparison, and gffcompare for configured
species. If no species IDs are provided, all species in species.json are used.
EOF
}

SELECTED_SPECIES=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --analysis-dir)
            ANALYSIS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            SELECTED_SPECIES+=("$1")
            shift
            ;;
    esac
done

if [[ ${#SELECTED_SPECIES[@]} -eq 0 ]]; then
    mapfile -t SELECTED_SPECIES < <(
        PYTHONPATH="$SCRIPT_DIR" "$PYTHON_BIN" -c 'from analysis_config import SPECIES_IDS; print("\n".join(SPECIES_IDS))'
    )
fi

if [[ ${#SELECTED_SPECIES[@]} -gt 0 ]]; then
    SPECIES_ARG=$(IFS=,; printf '%s' "${SELECTED_SPECIES[*]}")
    PYTHONPATH="$SCRIPT_DIR" "$PYTHON_BIN" - "$SPECIES_ARG" <<'PY'
import sys
from analysis_config import SPECIES_IDS

selected = [item for item in sys.argv[1].split(",") if item]
unknown = sorted(set(selected) - set(SPECIES_IDS))
if unknown:
    raise SystemExit(f"ERROR: unknown or excluded species ID(s): {', '.join(unknown)}")
PY
fi

STATS_DIR="$ANALYSIS/stats"
COMPARE_DIR="$ANALYSIS/compare"
TCOMPARE_DIR="$ANALYSIS/tcompare"
LOG_DIR="$ANALYSIS/logs"
mkdir -p "$STATS_DIR" "$COMPARE_DIR" "$TCOMPARE_DIR" "$LOG_DIR"

find_annotation() {
    local sp=$1
    local state=$2
    local candidate

    for candidate in \
        "$ANALYSIS/$sp.$state.gff" \
        "$ANALYSIS/$sp.$state.gff3" \
        "$ANALYSIS/$sp.$state.gff.gz" \
        "$ANALYSIS/$sp.$state.gff3.gz"
    do
        if [[ -f "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

require_annotation() {
    local sp=$1
    local state=$2
    local path

    if ! path=$(find_annotation "$sp" "$state"); then
        echo "ERROR: missing $state annotation for $sp in $ANALYSIS" >&2
        exit 1
    fi
    printf '%s\n' "$path"
}

run_stat() {
    local sp=$1
    local state=$2
    local gff=$3
    local out="$STATS_DIR/$sp.$state"
    local log="$LOG_DIR/$sp.$state.agat_stats.log"

    echo "  [$sp] $state statistics"
    if agat_sp_statistics.pl --gff "$gff" -o "$out" > "$log" 2>&1; then
        grep -E "(Job done|ERROR|Error|error)" "$log" || true
    else
        echo "ERROR: AGAT statistics failed for $sp $state. See $log" >&2
        tail -40 "$log" >&2 || true
        exit 1
    fi
}

echo "============================================"
echo "Step 1: AGAT statistics"
echo "============================================"

for sp in "${SELECTED_SPECIES[@]}"; do
    before=$(require_annotation "$sp" "before")
    after=$(require_annotation "$sp" "after")

    echo ""
    echo "--- $sp ---"
    run_stat "$sp" "before" "$before" &
    run_stat "$sp" "after" "$after" &
    wait
done

echo ""
echo "============================================"
echo "Step 2: AGAT comparison"
echo "============================================"

for sp in "${SELECTED_SPECIES[@]}"; do
    before=$(require_annotation "$sp" "before")
    after=$(require_annotation "$sp" "after")
    out="$COMPARE_DIR/$sp"
    log="$LOG_DIR/$sp.agat_compare.log"

    echo "--- $sp ---"
    mkdir -p "$out"
    if agat_sp_compare_two_annotations.pl --gff1 "$before" --gff2 "$after" -o "$out" > "$log" 2>&1; then
        grep -E "(Job done|ERROR|Error|error|Summary)" "$log" || true
    else
        echo "ERROR: AGAT comparison failed for $sp. See $log" >&2
        tail -40 "$log" >&2 || true
        exit 1
    fi
done

echo ""
echo "============================================"
echo "Step 3: gffcompare"
echo "============================================"

for sp in "${SELECTED_SPECIES[@]}"; do
    before=$(require_annotation "$sp" "before")
    after=$(require_annotation "$sp" "after")
    out="$TCOMPARE_DIR/$sp"
    log="$LOG_DIR/$sp.gffcompare.log"

    echo "--- $sp ---"
    if gffcompare -r "$before" -o "$out" "$after" > "$log" 2>&1; then
        tail -3 "$log"
    else
        echo "ERROR: gffcompare failed for $sp. See $log" >&2
        tail -40 "$log" >&2 || true
        exit 1
    fi
done

echo ""
echo "============================================"
echo "All analyses complete"
echo "============================================"
