#!/bin/bash
#
# Script to run radiative transfer for disks that completed disk model setup
# Usage: ./run_radiative_transfer.sh [start_id] [end_id] [rt_program_dir]
#

# Default values
START_ID=${1:-1}
END_ID=${2:-150}
RT_PROGRAM_DIR=${3:-"../radiative_transfer"}  # Path to radiative transfer program directory
DISK_RESULTS_DIR="results"
LOG_DIR="logs_rt"
RT_EXECUTABLE="radiative_transfer"

# Create log directory for RT runs
mkdir -p "$LOG_DIR"

# Get absolute paths
MAIN_DIR=$(pwd)
RT_PROGRAM_DIR=$(cd "$RT_PROGRAM_DIR" && pwd)

# Check if RT executable exists
if [ ! -f "${RT_PROGRAM_DIR}/${RT_EXECUTABLE}" ]; then
    echo "Error: Radiative transfer executable not found at ${RT_PROGRAM_DIR}/${RT_EXECUTABLE}"
    echo "Please compile the radiative transfer code first:"
    echo "  cd ${RT_PROGRAM_DIR}"
    echo "  make clean"
    echo "  make"
    exit 1
fi

# Count total disks to process
TOTAL_DISKS=$((END_ID - START_ID + 1))
echo "=========================================="
echo "Running radiative transfer for IDs $START_ID to $END_ID"
echo "Total disks to process: $TOTAL_DISKS"
echo "RT program directory: $RT_PROGRAM_DIR"
echo "=========================================="

# Loop through disk IDs
COMPLETED=0
FAILED=0
SKIPPED=0

for disk_id in $(seq $START_ID $END_ID); do
    echo ""
    echo "Processing disk ID: $disk_id ($(($disk_id - $START_ID + 1))/$TOTAL_DISKS)"
    echo "------------------------------------------"
    
    DISK_DIR="${MAIN_DIR}/${DISK_RESULTS_DIR}/disk_${disk_id}"
    
    # Check if disk directory exists (i.e., disk model completed successfully)
    if [ ! -d "$DISK_DIR" ]; then
        echo "⚠ Skipping: Disk directory not found (disk model may have failed)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    # Check if required input files exist
    if [ ! -f "${DISK_DIR}/amr_grid.inp" ] || [ ! -f "${DISK_DIR}/dust_density.inp" ]; then
        echo "⚠ Skipping: Required input files missing (disk model may have failed)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    # Change to disk directory
    cd "$DISK_DIR"
    
    # Run radiative transfer
    START_TIME=$(date +%s)
    "${RT_PROGRAM_DIR}/${RT_EXECUTABLE}" > "${MAIN_DIR}/${LOG_DIR}/disk_${disk_id}_rt.log" 2>&1
    EXIT_CODE=$?
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    # Return to main directory
    cd "${MAIN_DIR}"
    
    # Check if successful
    if [ $EXIT_CODE -eq 0 ]; then
        echo "✓ Completed successfully in ${ELAPSED}s"
        COMPLETED=$((COMPLETED + 1))
    else
        echo "✗ Failed with exit code $EXIT_CODE"
        echo "  Check log: ${LOG_DIR}/disk_${disk_id}_rt.log"
        FAILED=$((FAILED + 1))
    fi
done

# Summary
echo ""
echo "=========================================="
echo "RADIATIVE TRANSFER SUMMARY"
echo "=========================================="
echo "Completed: $COMPLETED"
echo "Failed:    $FAILED"
echo "Skipped:   $SKIPPED"
echo "Total:     $TOTAL_DISKS"
echo "=========================================="

# Exit with error if any failed
if [ $FAILED -gt 0 ]; then
    exit 1
fi
