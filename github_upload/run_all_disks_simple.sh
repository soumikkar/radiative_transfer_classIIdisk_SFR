#!/bin/bash
#
# Script to run disk model for multiple parameter sets
# Runs from main directory, moves output files after completion
# Usage: ./run_all_disks_simple.sh [start_id] [end_id]
#

# Default values
START_ID=${1:-1}
END_ID=${2:-150}
PARAM_FILE="disk_parameters.tsv"
EXECUTABLE="./disk_model"
OUTPUT_DIR="results"
LOG_DIR="logs"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "Error: Parameter file $PARAM_FILE not found!"
    exit 1
fi

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable $EXECUTABLE not found!"
    echo "Please compile the code first using: make"
    exit 1
fi

# Count total disks to process
TOTAL_DISKS=$((END_ID - START_ID + 1))
echo "=========================================="
echo "Running disk models for IDs $START_ID to $END_ID"
echo "Total disks to process: $TOTAL_DISKS"
echo "=========================================="

# Loop through disk IDs
COMPLETED=0
FAILED=0

for disk_id in $(seq $START_ID $END_ID); do
    echo ""
    echo "Processing disk ID: $disk_id ($(($disk_id - $START_ID + 1))/$TOTAL_DISKS)"
    echo "------------------------------------------"
    
    # Run the model from main directory
    START_TIME=$(date +%s)
    echo "$disk_id" | $EXECUTABLE > "${LOG_DIR}/disk_${disk_id}.log" 2>&1
    EXIT_CODE=$?
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    # Check if successful
    if [ $EXIT_CODE -eq 0 ]; then
        # Create disk-specific directory
        DISK_DIR="${OUTPUT_DIR}/disk_${disk_id}"
        mkdir -p "$DISK_DIR"
        
        # Move output files to disk directory
        # Adjust these patterns based on what your program outputs
        mv -f *.inp "$DISK_DIR/" 2>/dev/null
        mv -f *.out "$DISK_DIR/" 2>/dev/null
        mv -f *.dat "$DISK_DIR/" 2>/dev/null
        mv -f frequency.dat "$DISK_DIR/" 2>/dev/null
        mv -f wavelength_micron.inp "$DISK_DIR/" 2>/dev/null
        mv -f stars.inp "$DISK_DIR/" 2>/dev/null
        mv -f amr_grid.inp "$DISK_DIR/" 2>/dev/null
        mv -f dust_density.inp "$DISK_DIR/" 2>/dev/null
        mv -f dust_temperature.dat "$DISK_DIR/" 2>/dev/null
        
        # Add any other output file patterns your program creates
        
        echo "✓ Completed successfully in ${ELAPSED}s"
        COMPLETED=$((COMPLETED + 1))
    else
        echo "✗ Failed with exit code $EXIT_CODE"
        echo "  Check log: ${LOG_DIR}/disk_${disk_id}.log"
        FAILED=$((FAILED + 1))
    fi
done

# Summary
echo ""
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "Completed: $COMPLETED"
echo "Failed:    $FAILED"
echo "Total:     $TOTAL_DISKS"
echo "=========================================="

# Exit with error if any failed
if [ $FAILED -gt 0 ]; then
    exit 1
fi
