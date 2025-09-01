#!/bin/bash

#!/bin/bash

# SMN CNV Pipeline - Example Batch Processing Script
# Simple example for processing multiple BAM files

# Configuration - UPDATE THESE PATHS
BAM_DIR="/path/to/your/bam/files"
OUTPUT_DIR="output/$(date +%Y%m%d_%H%M%S)"
CONFIG_FILE="config/config.json"

# Processing options
BATCH_SIZE=20
PARALLEL_JOBS=4

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${GREEN}=== SMN CNV Pipeline - Batch Processing ===${NC}"
echo "Start time: $(date)"

# Validate inputs
if [ ! -d "$BAM_DIR" ]; then
    echo -e "${RED}Error: BAM directory not found: $BAM_DIR${NC}"
    echo "Please update BAM_DIR in this script"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo -e "${RED}Error: Configuration file not found: $CONFIG_FILE${NC}"
    exit 1
fi

# Count BAM files
BAM_COUNT=$(find "$BAM_DIR" -name "*.bam" | wc -l)
echo -e "${YELLOW}Found $BAM_COUNT BAM files in $BAM_DIR${NC}"

if [ $BAM_COUNT -eq 0 ]; then
    echo -e "${RED}No BAM files found in directory${NC}"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
echo -e "${YELLOW}Output directory: $OUTPUT_DIR${NC}"

# Estimate processing time
echo -e "${YELLOW}Estimated processing time: $(python3 -c "
import sys
sys.path.append('src')
try:
    from utils import estimate_runtime
    est = estimate_runtime($BAM_COUNT, 30)
    print(f'{est[\"total_hours\"]:.1f} hours')
except:
    print('~$(($BAM_COUNT * 3 / 60)) hours (rough estimate)')
")${NC}"

# Ask for confirmation
read -p "Continue with processing? (y/N): " confirm
if [[ ! $confirm =~ ^[Yy]$ ]]; then
    echo "Processing cancelled"
    exit 0
fi

# Run the pipeline
echo -e "${GREEN}Starting batch processing...${NC}"

python3 smn_pipeline.py \
    --config "$CONFIG_FILE" \
    --input "$BAM_DIR" \
    --output "$OUTPUT_DIR" \
    --batch

# Check if processing was successful
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✅ Batch processing completed successfully!${NC}"
    
    # Display results summary
    echo ""
    echo "=== RESULTS SUMMARY ==="
    
    # Count results
    HTML_REPORTS=$(find "$OUTPUT_DIR" -name "*.html" | wc -l)
    TSV_FILES=$(find "$OUTPUT_DIR" -name "*.tsv" | wc -l)
    SNAPSHOTS=$(find "$OUTPUT_DIR" -name "*.png" | wc -l)
    
    echo "HTML reports generated: $HTML_REPORTS"
    echo "TSV files created: $TSV_FILES"
    echo "IGV snapshots: $SNAPSHOTS"
    
    # Find main results files
    BATCH_REPORT=$(find "$OUTPUT_DIR" -name "batch_report_*.html" | head -1)
    BATCH_TSV=$(find "$OUTPUT_DIR" -name "batch_results_*.tsv" | head -1)
    
    if [ -f "$BATCH_REPORT" ]; then
        echo -e "\n${GREEN}📊 Main batch report: $BATCH_REPORT${NC}"
        echo "Open this file in a web browser to view results"
    fi
    
    if [ -f "$BATCH_TSV" ]; then
        echo -e "${GREEN}📋 Results table: $BATCH_TSV${NC}"
        
        # Quick results summary
        if command -v tail &> /dev/null && command -v cut &> /dev/null; then
            echo ""
            echo "=== QUICK RESULTS SUMMARY ==="
            
            # Count SMA risk categories
            HIGH_RISK=$(tail -n +2 "$BATCH_TSV" | cut -f2 | grep -c "HIGH_RISK" || echo "0")
            CARRIERS=$(tail -n +2 "$BATCH_TSV" | cut -f2 | grep -c "CARRIER" || echo "0")
            LOW_RISK=$(tail -n +2 "$BATCH_TSV" | cut -f2 | grep -c "LOW_RISK" || echo "0")
            UNCERTAIN=$(tail -n +2 "$BATCH_TSV" | cut -f2 | grep -c "UNCERTAIN" || echo "0")
            
            echo "High Risk (SMA): $HIGH_RISK samples"
            echo "Carriers: $CARRIERS samples"
            echo "Low Risk: $LOW_RISK samples"
            echo "Uncertain: $UNCERTAIN samples"
            
            if [ $HIGH_RISK -gt 0 ]; then
                echo -e "${RED}⚠️ URGENT: $HIGH_RISK high-risk samples detected${NC}"
                echo "Review these samples immediately for clinical action"
            fi
            
            if [ $CARRIERS -gt 0 ]; then
                echo -e "${YELLOW}ℹ️ $CARRIERS carrier samples detected${NC}"
                echo "Genetic counseling recommended for these samples"
            fi
        fi
    fi
    
    echo ""
    echo "=== NEXT STEPS ==="
    echo "1. Review the batch report in a web browser"
    echo "2. Follow up on high-risk and carrier samples"
    echo "3. Archive results to permanent storage"
    
    if [ $BAM_COUNT -gt 50 ]; then
        echo "4. Consider updating ML thresholds with this large dataset"
    fi
    
    echo ""
    echo "End time: $(date)"
    
else
    echo -e "${RED}❌ Batch processing failed${NC}"
    echo "Check the log files in $OUTPUT_DIR/logs/ for details"
    exit 1
fi
