#!/bin/bash
# run_radiative_transfer.sh
# Usage: ./run_radiative_transfer.sh [start_id] [end_id]

START_ID=${1:-1}
END_ID=${2:-150}

RESULTS_DIR="/media/admin1/DATA/Input_Params_bashscript/Taurus_input_code/results"
RT_EXECUTABLE="/media/admin1/DATA/Input_Params_bashscript/Taurus_input_code/radiative_transfer/radiative_transfer"
LOG_DIR="/media/admin1/DATA/Input_Params_bashscript/Taurus_input_code/logs_rt"

mkdir -p "$LOG_DIR"

COMPLETED=0
FAILED=0
TOTAL=$((END_ID - START_ID + 1))

echo "=========================================="
echo "Running radiative transfer for IDs $START_ID to $END_ID"
echo "=========================================="

for disk_id in $(seq $START_ID $END_ID); do
    DISK_DIR="${RESULTS_DIR}/disk_${disk_id}"

    echo ""
    echo "Processing disk ID: $disk_id ($(($disk_id - $START_ID + 1))/$TOTAL)"
    echo "------------------------------------------"

    if [ ! -d "$DISK_DIR" ]; then
        echo "✗ Directory $DISK_DIR not found, skipping..."
        FAILED=$((FAILED + 1))
        continue
    fi

    START_TIME=$(date +%s)
    (cd "$DISK_DIR" && $RT_EXECUTABLE) > "${LOG_DIR}/rt_disk_${disk_id}.log" 2>&1
    EXIT_CODE=$?
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))

    if [ $EXIT_CODE -eq 0 ]; then
        echo "✓ Completed in ${ELAPSED}s"
        COMPLETED=$((COMPLETED + 1))
    else
        echo "✗ Failed with exit code $EXIT_CODE"
        echo "  Check log: ${LOG_DIR}/rt_disk_${disk_id}.log"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "Completed: $COMPLETED"
echo "Failed:    $FAILED"
echo "Total:     $TOTAL"
echo "=========================================="

if [ $FAILED -gt 0 ]; then
    exit 1
fi
