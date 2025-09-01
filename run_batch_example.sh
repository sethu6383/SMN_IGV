#!/bin/bash

# SMN CNV Pipeline - Advanced Batch Processing Script
# Handles large-scale batch processing with monitoring and error recovery

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config/config.json"
PIPELINE_SCRIPT="${SCRIPT_DIR}/smn_pipeline.py"

# Default values
BATCH_SIZE=20
MAX_PARALLEL=4
RETRY_FAILED=true
CLEANUP_TEMP=true
GENERATE_QC_REPORT=true

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

SMN CNV Pipeline Batch Processor

OPTIONS:
    -i, --input DIR          Input directory containing BAM files
    -o, --output DIR         Output directory for results
    -c, --config FILE        Configuration file (default: config/config.json)
    -b, --batch-size N       Number of samples per batch (default: 20)
    -p, --parallel N         Maximum parallel processes (default: 4)
    -r, --retry              Retry failed samples (default: true)
    -q, --qc-report         Generate QC report (default: true)
    --no-cleanup            Don't cleanup temporary files
    --dry-run               Show what would be processed without running
    -h, --help              Show this help message

EXAMPLES:
    # Basic batch processing
    $0 -i /data/bam_files -o /results/batch_001

    # Large batch with custom settings
    $0 -i /data/large_cohort -o /results/cohort_analysis \\
       -b 50 -p 8 --qc-report

    # Dry run to check files
    $0 -i /data/bam_files -o /results/test --dry-run

EOF
}

# Function to log messages
log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case $level in
        "INFO")
            echo -e "${BLUE}[INFO]${NC} ${timestamp} - ${message}"
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} ${timestamp} - ${message}"
            ;;
        "ERROR")
            echo -e "${RED}[ERROR]${NC} ${timestamp} - ${message}"
            ;;
        "SUCCESS")
            echo -e "${GREEN}[SUCCESS]${NC} ${timestamp} - ${message}"
            ;;
    esac
}

# Function to check prerequisites
check_prerequisites() {
    log_message "INFO" "Checking prerequisites..."
    
    # Check if Python script exists
    if [ ! -f "$PIPELINE_SCRIPT" ]; then
        log_message "ERROR" "Pipeline script not found: $PIPELINE_SCRIPT"
        exit 1
    fi
    
    # Check if config exists
    if [ ! -f "$CONFIG_FILE" ]; then
        log_message "ERROR" "Configuration file not found: $CONFIG_FILE"
        exit 1
    fi
    
    # Check Python dependencies
    if ! python3 -c "import pysam, numpy, pandas, scipy, sklearn, plotly, optuna" 2>/dev/null; then
        log_message "ERROR" "Missing Python dependencies. Run: pip3 install -r requirements.txt"
        exit 1
    fi
    
    # Check samtools
    if ! command -v samtools &> /dev/null; then
        log_message "WARN" "samtools not found in PATH"
    fi
    
    log_message "SUCCESS" "Prerequisites check completed"
}

# Function to validate BAM files
validate_bam_files() {
    local input_dir=$1
    local bam_files=()
    
    log_message "INFO" "Validating BAM files in $input_dir..."
    
    # Find BAM files
    while IFS= read -r -d '' file; do
        bam_files+=("$file")
    done < <(find "$input_dir" -name "*.bam" -type f -print0)
    
    if [ ${#bam_files[@]} -eq 0 ]; then
        log_message "ERROR" "No BAM files found in $input_dir"
        exit 1
    fi
    
    log_message "INFO" "Found ${#bam_files[@]} BAM files"
    
    # Validate each BAM file
    local valid_files=()
    local invalid_files=()
    
    for bam_file in "${bam_files[@]}"; do
        # Check if file is readable
        if [ ! -r "$bam_file" ]; then
            log_message "WARN" "Cannot read BAM file: $bam_file"
            invalid_files+=("$bam_file")
            continue
        fi
        
        # Check if indexed
        if [ ! -f "${bam_file}.bai" ] && [ ! -f "${bam_file%.*}.bai" ]; then
            log_message "WARN" "BAM file not indexed: $bam_file"
            log_message "INFO" "Attempting to index..."
            
            if command -v samtools &> /dev/null; then
                if samtools index "$bam_file"; then
                    log_message "SUCCESS" "Indexed: $bam_file"
                else
                    log_message "ERROR" "Failed to index: $bam_file"
                    invalid_files+=("$bam_file")
                    continue
                fi
            else
                log_message "WARN" "samtools not available for indexing"
            fi
        fi
        
        valid_files+=("$bam_file")
    done
    
    if [ ${#invalid_files[@]} -gt 0 ]; then
        log_message "WARN" "Invalid BAM files found: ${#invalid_files[@]}"
        for invalid_file in "${invalid_files[@]}"; do
            log_message "WARN" "  - $invalid_file"
        done
    fi
    
    if [ ${#valid_files[@]} -eq 0 ]; then
        log_message "ERROR" "No valid BAM files found"
        exit 1
    fi
    
    log_message "SUCCESS" "Validation completed. ${#valid_files[@]} valid BAM files"
    printf '%s\0' "${valid_files[@]}"
}

# Function to process batch
process_batch() {
    local batch_files=("$@")
    local batch_num=$((batch_counter++))
    local batch_output="${OUTPUT_DIR}/batch_${batch_num}"
    
    log_message "INFO" "Processing batch $batch_num with ${#batch_files[@]} files"
    
    # Create batch-specific output directory
    mkdir -p "$batch_output"
    
    # Create file list for this batch
    local batch_list="${batch_output}/file_list.txt"
    printf '%s\n' "${batch_files[@]}" > "$batch_list"
    
    # Run pipeline for this batch
    local batch_start_time=$(date +%s)
    
    if python3 "$PIPELINE_SCRIPT" \
        --config "$CONFIG_FILE" \
        --input "$batch_list" \
        --output "$batch_output" \
        --batch; then
        
        local batch_end_time=$(date +%s)
        local batch_duration=$((batch_end_time - batch_start_time))
        
        log_message "SUCCESS" "Batch $batch_num completed in ${batch_duration}s"
        
        # Record successful batch
        echo "batch_${batch_num},${#batch_files[@]},${batch_duration},success" >> "${OUTPUT_DIR}/batch_summary.csv"
        
        return 0
    else
        log_message "ERROR" "Batch $batch_num failed"
        
        # Record failed batch
        echo "batch_${batch_num},${#batch_files[@]},0,failed" >> "${OUTPUT_DIR}/batch_summary.csv"
        
        # Handle retry logic
        if [ "$RETRY_FAILED" = true ]; then
            log_message "INFO" "Retrying failed samples individually..."
            retry_failed_samples "${batch_files[@]}"
        fi
        
        return 1
    fi
}

# Function to retry failed samples individually
retry_failed_samples() {
    local failed_files=("$@")
    local retry_dir="${OUTPUT_DIR}/retries"
    
    mkdir -p "$retry_dir"
    
    for bam_file in "${failed_files[@]}"; do
        local sample_id=$(basename "$bam_file" .bam)
        local sample_output="${retry_dir}/${sample_id}"
        
        log_message "INFO" "Retrying sample: $sample_id"
        
        if python3 "$PIPELINE_SCRIPT" \
            --config "$CONFIG_FILE" \
            --input "$bam_file" \
            --output "$sample_output" \
            --sample-id "$sample_id"; then
            
            log_message "SUCCESS" "Retry successful: $sample_id"
        else
            log_message "ERROR" "Retry failed: $sample_id"
            echo "$sample_id,retry_failed" >> "${OUTPUT_DIR}/failed_samples.txt"
        fi
    done
}

# Function to monitor resources
monitor_resources() {
    local output_file="${OUTPUT_DIR}/resource_usage.log"
    
    while true; do
        local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        local cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d'%' -f1)
        local mem_usage=$(free | grep Mem | awk '{printf "%.1f", $3/$2 * 100.0}')
        local disk_usage=$(df "$OUTPUT_DIR" | tail -1 | awk '{print $5}' | cut -d'%' -f1)
        
        echo "$timestamp,$cpu_usage,$mem_usage,$disk_usage" >> "$output_file"
        
        sleep 60  # Monitor every minute
    done &
    
    monitor_pid=$!
    log_message "INFO" "Resource monitoring started (PID: $monitor_pid)"
}

# Function to generate final QC report
generate_qc_report() {
    log_message "INFO" "Generating QC report..."
    
    # Count results
    local total_samples=$(find "$OUTPUT_DIR" -name "*_summary.tsv" | wc -l)
    local failed_samples=0
    
    if [ -f "${OUTPUT_DIR}/failed_samples.txt" ]; then
        failed_samples=$(wc -l < "${OUTPUT_DIR}/failed_samples.txt")
    fi
    
    local success_rate=$(( (total_samples - failed_samples) * 100 / total_samples ))
    
    # Create QC report
    cat > "${OUTPUT_DIR}/qc_report.txt" << EOF
SMN CNV Pipeline - Batch QC Report
Generated: $(date)

=== BATCH SUMMARY ===
Total samples processed: $total_samples
Successful samples: $((total_samples - failed_samples))
Failed samples: $failed_samples
Success rate: ${success_rate}%

=== PROCESSING TIME ===
Start time: $start_time
End time: $(date)
Total duration: $(($(date +%s) - $(date -d "$start_time" +%s)))s

=== OUTPUT FILES ===
Batch reports: $(find "$OUTPUT_DIR" -name "batch_report_*.html" | wc -l)
Individual reports: $(find "$OUTPUT_DIR" -name "*_report.html" | wc -l)
TSV files: $(find "$OUTPUT_DIR" -name "*.tsv" | wc -l)
IGV snapshots: $(find "$OUTPUT_DIR" -name "*.png" | wc -l)

=== RESOURCE USAGE ===
Peak memory usage: $(tail -1 "${OUTPUT_DIR}/resource_usage.log" 2>/dev/null | cut -d',' -f3 || echo "N/A")%
Disk space used: $(du -sh "$OUTPUT_DIR" | cut -f1)

=== RECOMMENDATIONS ===
EOF

    # Add recommendations based on results
    if [ $success_rate -lt 90 ]; then
        echo "- Review failed samples in failed_samples.txt" >> "${OUTPUT_DIR}/qc_report.txt"
        echo "- Check BAM file quality and indexing" >> "${OUTPUT_DIR}/qc_report.txt"
    fi
    
    if [ $total_samples -gt 100 ]; then
        echo "- Consider updating ML thresholds with this large dataset" >> "${OUTPUT_DIR}/qc_report.txt"
        echo "- Review performance metrics in models/performance_history.json" >> "${OUTPUT_DIR}/qc_report.txt"
    fi
    
    echo "- Review consolidated results in batch_report_*.html" >> "${OUTPUT_DIR}/qc_report.txt"
    
    log_message "SUCCESS" "QC report generated: ${OUTPUT_DIR}/qc_report.txt"
}

# Function to cleanup temporary files
cleanup_temp_files() {
    if [ "$CLEANUP_TEMP" = true ]; then
        log_message "INFO" "Cleaning up temporary files..."
        
        # Remove temporary files older than 1 day
        find "$OUTPUT_DIR" -name "*.tmp" -mtime +1 -delete 2>/dev/null || true
        find "$OUTPUT_DIR" -name "*.bat" -mtime +1 -delete 2>/dev/null || true
        find "${SCRIPT_DIR}/temp" -name "*" -mtime +1 -delete 2>/dev/null || true
        
        log_message "INFO" "Temporary files cleaned"
    fi
}

# Function to handle script termination
cleanup_on_exit() {
    log_message "INFO" "Script terminating..."
    
    # Stop resource monitoring
    if [ ! -z "$monitor_pid" ]; then
        kill $monitor_pid 2>/dev/null || true
    fi
    
    # Cleanup temp files
    cleanup_temp_files
    
    # Generate final report if processing was started
    if [ ! -z "$processing_started" ]; then
        generate_qc_report
    fi
    
    log_message "INFO" "Cleanup completed"
}

# Set trap for cleanup
trap cleanup_on_exit EXIT

# Parse command line arguments
INPUT_DIR=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -b|--batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        -p|--parallel)
            MAX_PARALLEL="$2"
            shift 2
            ;;
        -r|--retry)
            RETRY_FAILED=true
            shift
            ;;
        -q|--qc-report)
            GENERATE_QC_REPORT=true
            shift
            ;;
        --no-cleanup)
            CLEANUP_TEMP=false
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log_message "ERROR" "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    log_message "ERROR" "Input and output directories are required"
    usage
    exit 1
fi

# Validate input directory
if [ ! -d "$INPUT_DIR" ]; then
    log_message "ERROR" "Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Start logging
start_time=$(date)
log_message "INFO" "Starting SMN CNV batch processing"
log_message "INFO" "Input directory: $INPUT_DIR"
log_message "INFO" "Output directory: $OUTPUT_DIR"
log_message "INFO" "Configuration: $CONFIG_FILE"
log_message "INFO" "Batch size: $BATCH_SIZE"
log_message "INFO" "Max parallel: $MAX_PARALLEL"

# Check prerequisites
check_prerequisites

# Validate BAM files
log_message "INFO" "Validating BAM files..."
mapfile -d '' valid_bam_files < <(validate_bam_files "$INPUT_DIR")

if [ ${#valid_bam_files[@]} -eq 0 ]; then
    log_message "ERROR" "No valid BAM files found"
    exit 1
fi

log_message "INFO" "Found ${#valid_bam_files[@]} valid BAM files"

# Dry run check
if [ "$DRY_RUN" = true ]; then
    log_message "INFO" "DRY RUN - Files that would be processed:"
    for bam_file in "${valid_bam_files[@]}"; do
        echo "  - $bam_file"
    done
    
    log_message "INFO" "Total batches: $(( (${#valid_bam_files[@]} + BATCH_SIZE - 1) / BATCH_SIZE ))"
    log_message "INFO" "Estimated runtime: $(python3 -c "
import sys
sys.path.append('src')
from utils import estimate_runtime
est = estimate_runtime(${#valid_bam_files[@]}, 30)
print(f\"{est['total_hours']:.1f} hours\")
")"
    exit 0
fi

# Initialize batch processing
processing_started=true
batch_counter=1
failed_batches=()
successful_batches=()

# Create batch summary header
echo "batch_id,sample_count,duration_seconds,status" > "${OUTPUT_DIR}/batch_summary.csv"

# Start resource monitoring
if command -v top &> /dev/null && command -v free &> /dev/null; then
    monitor_resources
fi

# Process files in batches
log_message "INFO" "Starting batch processing..."

for (( i=0; i<${#valid_bam_files[@]}; i+=BATCH_SIZE )); do
    # Get batch files
    batch_files=("${valid_bam_files[@]:i:BATCH_SIZE}")
    
    # Process batch
    if process_batch "${batch_files[@]}"; then
        successful_batches+=($((batch_counter-1)))
    else
        failed_batches+=($((batch_counter-1)))
    fi
    
    # Progress update
    local processed=$((i + ${#batch_files[@]}))
    local total=${#valid_bam_files[@]}
    local progress=$(( processed * 100 / total ))
    
    log_message "INFO" "Progress: $processed/$total samples ($progress%)"
    
    # Optional: pause between batches to prevent system overload
    if [ ${#batch_files[@]} -eq $BATCH_SIZE ] && [ $i -lt $((${#valid_bam_files[@]} - BATCH_SIZE)) ]; then
        sleep 10
    fi
done

# Stop resource monitoring
if [ ! -z "$monitor_pid" ]; then
    kill $monitor_pid 2>/dev/null || true
fi

# Final summary
log_message "INFO" "Batch processing completed"
log_message "INFO" "Successful batches: ${#successful_batches[@]}"
log_message "INFO" "Failed batches: ${#failed_batches[@]}"

if [ ${#failed_batches[@]} -gt 0 ]; then
    log_message "WARN" "Failed batch numbers: ${failed_batches[*]}"
fi

# Generate QC report
if [ "$GENERATE_QC_REPORT" = true ]; then
    generate_qc_report
fi

# Final consolidation
log_message "INFO" "Consolidating results..."

# Merge all TSV files
find "$OUTPUT_DIR" -name "*_summary.tsv" -exec cat {} \; | \
    awk 'NR==1 || !/^sample_id/' > "${OUTPUT_DIR}/consolidated_results.tsv"

# Count final results
total_processed=$(tail -n +2 "${OUTPUT_DIR}/consolidated_results.tsv" | wc -l)
high_risk_count=$(tail -n +2 "${OUTPUT_DIR}/consolidated_results.tsv" | grep -c "HIGH_RISK" || echo "0")
carrier_count=$(tail -n +2 "${OUTPUT_DIR}/consolidated_results.tsv" | grep -c "CARRIER" || echo "0")

log_message "SUCCESS" "Processing completed!"
log_message "INFO" "Final statistics:"
log_message "INFO" "  - Total samples processed: $total_processed"
log_message "INFO" "  - High risk samples: $high_risk_count"
log_message "INFO" "  - Carrier samples: $carrier_count"
log_message "INFO" "  - Success rate: $(( (total_processed * 100) / ${#valid_bam_files[@]} ))%"

# Performance summary
end_time=$(date +%s)
total_duration=$((end_time - $(date -d "$start_time" +%s)))
samples_per_minute=$(echo "scale=2; $total_processed * 60 / $total_duration" | bc -l 2>/dev/null || echo "N/A")

log_message "INFO" "Performance: $samples_per_minute samples/minute"

# Final recommendations
echo ""
log_message "INFO" "=== NEXT STEPS ==="
log_message "INFO" "1. Review consolidated results: ${OUTPUT_DIR}/consolidated_results.tsv"
log_message "INFO" "2. Check batch report: ${OUTPUT_DIR}/batch_report_*.html"
log_message "INFO" "3. Follow up on high-risk samples for clinical action"
log_message "INFO" "4. Consider updating ML thresholds if >50 new samples processed"

if [ $high_risk_count -gt 0 ]; then
    log_message "WARN" "⚠️  $high_risk_count high-risk samples detected - review urgently"
fi

if [ ${#failed_batches[@]} -gt 0 ]; then
    log_message "WARN" "⚠️  Some batches failed - check logs for details"
fi

log_message "SUCCESS" "🧬 SMN CNV batch processing completed successfully!"
