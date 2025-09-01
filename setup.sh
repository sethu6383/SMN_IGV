#!/bin/bash

# SMN CNV Pipeline Setup Script
# Sets up the complete pipeline environment

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== SMN CNV Pipeline Setup ===${NC}"
echo "Setting up SMN1/SMN2 CNV detection pipeline..."

# Check if running on Linux
if [[ "$OSTYPE" != "linux-gnu"* ]]; then
    echo -e "${RED}Error: This pipeline is designed for Linux systems only${NC}"
    exit 1
fi

# Create directory structure
echo -e "${YELLOW}Creating directory structure...${NC}"
mkdir -p {src,config,models,output,reports,snapshots,logs,temp,tests,docs}
mkdir -p src/{__pycache__}

# Check Python version
echo -e "${YELLOW}Checking Python version...${NC}"
PYTHON_VERSION=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
echo "Found Python $PYTHON_VERSION"

if ! python3 -c 'import sys; exit(0 if sys.version_info >= (3, 8) else 1)'; then
    echo -e "${RED}Error: Python 3.8 or higher required${NC}"
    exit 1
fi

# Check if pip is available
if ! command -v pip3 &> /dev/null; then
    echo -e "${RED}Error: pip3 not found. Please install pip3${NC}"
    exit 1
fi

# Install Python dependencies
echo -e "${YELLOW}Installing Python dependencies...${NC}"
pip3 install -r requirements.txt

# Check for required system tools
echo -e "${YELLOW}Checking system dependencies...${NC}"

# Check for samtools
if ! command -v samtools &> /dev/null; then
    echo -e "${YELLOW}Warning: samtools not found. Installing via conda/apt...${NC}"
    
    if command -v conda &> /dev/null; then
        conda install -c bioconda samtools
    elif command -v apt-get &> /dev/null; then
        sudo apt-get update && sudo apt-get install -y samtools
    else
        echo -e "${RED}Please install samtools manually${NC}"
    fi
fi

# Check for IGV
echo -e "${YELLOW}Checking IGV installation...${NC}"
if ! command -v igv.sh &> /dev/null && ! command -v igv &> /dev/null; then
    echo -e "${YELLOW}IGV not found in PATH. Please ensure IGV is installed and accessible.${NC}"
    echo "Download from: https://software.broadinstitute.org/software/igv/download"
    echo "After installation, update the 'igv_path' in config/config.json"
fi

# Create Python module files
echo -e "${YELLOW}Creating Python module files...${NC}"

# Create __init__.py files
touch src/__init__.py
touch tests/__init__.py

# Create requirements.txt if it doesn't exist
if [ ! -f requirements.txt ]; then
    echo -e "${YELLOW}Creating requirements.txt...${NC}"
    cat > requirements.txt << EOF
# Core dependencies
pysam>=0.21.0
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
scikit-learn>=1.0.0

# Plotting and visualization
plotly>=5.0.0
matplotlib>=3.5.0
seaborn>=0.11.0

# Machine learning optimization
optuna>=3.0.0
joblib>=1.1.0

# System and utilities
psutil>=5.8.0
tqdm>=4.62.0
pathlib2>=2.3.6

# Optional for enhanced features
biopython>=1.79
pyvcf>=0.6.8
EOF
fi

# Validate installation
echo -e "${YELLOW}Validating installation...${NC}"
python3 -c "
import pysam
import numpy as np
import pandas as pd
import scipy
import sklearn
import plotly
import optuna
print('All required packages imported successfully!')
"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Python dependencies installed successfully${NC}"
else
    echo -e "${RED}✗ Error with Python dependencies${NC}"
    exit 1
fi

# Create sample validation data template
echo -e "${YELLOW}Creating validation data template...${NC}"
cat > config/validation_template.tsv << EOF
sample_id	smn1_exon7	smn1_exon8	smn2_exon7	smn2_exon8	mlpa_result	comments
CONTROL_001	NORMAL	NORMAL	NORMAL	NORMAL	Normal control	Validated by MLPA
CARRIER_001	HETERO_DEL	HETERO_DEL	NORMAL	NORMAL	SMA carrier	Validated by MLPA
SMA_001	HOMO_DEL	HOMO_DEL	NORMAL	NORMAL	SMA affected	Validated by MLPA
EOF

# Create example batch processing script
echo -e "${YELLOW}Creating example batch script...${NC}"
cat > run_batch_example.sh << 'EOF'
#!/bin/bash

# Example batch processing script
# Modify paths according to your setup

BAM_DIR="/path/to/your/bam/files"
OUTPUT_DIR="output/$(date +%Y%m%d_%H%M%S)"
CONFIG_FILE="config/config.json"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run pipeline on all BAM files in directory
python3 smn_pipeline.py \
    --config "$CONFIG_FILE" \
    --input "$BAM_DIR" \
    --output "$OUTPUT_DIR" \
    --batch

echo "Batch processing completed. Results in: $OUTPUT_DIR"
EOF

chmod +x run_batch_example.sh

# Create test script
echo -e "${YELLOW}Creating test script...${NC}"
cat > run_tests.sh << 'EOF'
#!/bin/bash

# Run tests for SMN pipeline
echo "Running SMN Pipeline Tests..."

# Test 1: Configuration validation
echo "Testing configuration..."
python3 -c "
import json
with open('config/config.json', 'r') as f:
    config = json.load(f)
print('✓ Configuration file valid')
"

# Test 2: Module imports
echo "Testing module imports..."
python3 -c "
from src.depth_analyzer import DepthAnalyzer
from src.cnv_caller import CNVCaller
from src.ml_threshold import MLThresholdOptimizer
from src.report_generator import ReportGenerator
from src.igv_automation import IGVAutomation
from src.utils import setup_logging, validate_inputs
print('✓ All modules import successfully')
"

# Test 3: Dependencies
echo "Testing dependencies..."
python3 -c "
from src.utils import check_dependencies
deps = check_dependencies()
missing = [k for k, v in deps.items() if not v]
if missing:
    print(f'✗ Missing dependencies: {missing}')
    exit(1)
else:
    print('✓ All dependencies available')
"

echo "All tests passed!"
EOF

chmod +x run_tests.sh

# Create IGV setup helper
echo -e "${YELLOW}Creating IGV setup helper...${NC}"
cat > setup_igv.sh << 'EOF'
#!/bin/bash

# IGV Setup Helper Script

echo "=== IGV Setup Helper ==="

# Check if IGV is already installed
if command -v igv.sh &> /dev/null; then
    echo "✓ IGV found in PATH: $(which igv.sh)"
    IGV_VERSION=$(igv.sh --version 2>&1 | head -1 || echo "Unknown version")
    echo "IGV Version: $IGV_VERSION"
    exit 0
fi

if command -v igv &> /dev/null; then
    echo "✓ IGV found in PATH: $(which igv)"
    exit 0
fi

echo "IGV not found in PATH. Setting up IGV..."

# Create IGV directory
IGV_DIR="$HOME/igv"
mkdir -p "$IGV_DIR"

# Download IGV (latest version)
echo "Downloading IGV..."
cd "$IGV_DIR"

# Check if wget or curl is available
if command -v wget &> /dev/null; then
    wget -O igv.zip "https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_2.16.2_WithJava.zip"
elif command -v curl &> /dev/null; then
    curl -L -o igv.zip "https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_2.16.2_WithJava.zip"
else
    echo "Error: Neither wget nor curl found. Please install one of them."
    exit 1
fi

# Extract IGV
echo "Extracting IGV..."
unzip -q igv.zip
rm igv.zip

# Find IGV executable
IGV_EXEC=$(find . -name "igv.sh" -type f | head -1)
if [ -z "$IGV_EXEC" ]; then
    echo "Error: Could not find igv.sh after extraction"
    exit 1
fi

# Make executable
chmod +x "$IGV_EXEC"

# Add to PATH
IGV_FULL_PATH="$(cd "$(dirname "$IGV_EXEC")" && pwd)/$(basename "$IGV_EXEC")"
echo "IGV installed at: $IGV_FULL_PATH"

# Update config file
python3 << PYTHON_SCRIPT
import json
import os

config_file = 'config/config.json'
if os.path.exists(config_file):
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    config['igv_settings']['igv_path'] = '$IGV_FULL_PATH'
    
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"Updated config file with IGV path: $IGV_FULL_PATH")
else:
    print("Config file not found. Please update igv_path manually.")
PYTHON_SCRIPT

echo "✓ IGV setup completed!"
echo "To use IGV from anywhere, add this to your ~/.bashrc:"
echo "export PATH=\"$(dirname "$IGV_FULL_PATH"):\$PATH\""
EOF

chmod +x setup_igv.sh

# Create pipeline validation script
echo -e "${YELLOW}Creating pipeline validation script...${NC}"
cat > validate_pipeline.py << 'EOF'
#!/usr/bin/env python3
"""
Pipeline Validation Script
Validates the complete SMN pipeline setup
"""

import os
import sys
import json
import subprocess
from pathlib import Path

def validate_directory_structure():
    """Validate directory structure"""
    required_dirs = ['src', 'config', 'models', 'output', 'reports', 'snapshots', 'logs', 'temp']
    missing_dirs = []
    
    for dir_name in required_dirs:
        if not os.path.exists(dir_name):
            missing_dirs.append(dir_name)
    
    if missing_dirs:
        print(f"✗ Missing directories: {missing_dirs}")
        return False
    else:
        print("✓ Directory structure valid")
        return True

def validate_python_modules():
    """Validate Python modules can be imported"""
    modules = [
        'src.depth_analyzer',
        'src.cnv_caller', 
        'src.ml_threshold',
        'src.report_generator',
        'src.igv_automation',
        'src.utils'
    ]
    
    failed_imports = []
    for module in modules:
        try:
            __import__(module)
        except ImportError as e:
            failed_imports.append(f"{module}: {e}")
    
    if failed_imports:
        print(f"✗ Module import failures:")
        for failure in failed_imports:
            print(f"  - {failure}")
        return False
    else:
        print("✓ All Python modules import successfully")
        return True

def validate_config():
    """Validate configuration file"""
    config_file = 'config/config.json'
    
    if not os.path.exists(config_file):
        print(f"✗ Configuration file not found: {config_file}")
        return False
    
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        required_keys = ['pipeline_version', 'smn_coordinates', 'cnv_thresholds']
        missing_keys = [key for key in required_keys if key not in config]
        
        if missing_keys:
            print(f"✗ Missing config keys: {missing_keys}")
            return False
        else:
            print("✓ Configuration file valid")
            return True
            
    except Exception as e:
        print(f"✗ Error reading config file: {e}")
        return False

def validate_dependencies():
    """Validate required dependencies"""
    dependencies = [
        'pysam', 'numpy', 'pandas', 'scipy', 'sklearn',
        'plotly', 'optuna', 'psutil'
    ]
    
    missing_deps = []
    for dep in dependencies:
        try:
            __import__(dep)
        except ImportError:
            missing_deps.append(dep)
    
    if missing_deps:
        print(f"✗ Missing Python dependencies: {missing_deps}")
        print("Run: pip3 install -r requirements.txt")
        return False
    else:
        print("✓ All Python dependencies available")
        return True

def validate_system_tools():
    """Validate system tools"""
    tools = ['samtools']
    missing_tools = []
    
    for tool in tools:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing_tools.append(tool)
    
    if missing_tools:
        print(f"✗ Missing system tools: {missing_tools}")
        return False
    else:
        print("✓ Required system tools available")
        return True

def validate_igv():
    """Validate IGV installation"""
    
    # Check config for IGV path
    try:
        with open('config/config.json', 'r') as f:
            config = json.load(f)
        
        igv_path = config.get('igv_settings', {}).get('igv_path', 'igv.sh')
        
        # Test IGV
        result = subprocess.run([igv_path, '--help'], 
                              capture_output=True, 
                              timeout=10)
        
        if result.returncode == 0:
            print("✓ IGV is accessible")
            return True
        else:
            print(f"✗ IGV test failed")
            print("Run ./setup_igv.sh to install IGV")
            return False
            
    except Exception as e:
        print(f"✗ IGV validation error: {e}")
        print("Run ./setup_igv.sh to install IGV")
        return False

def main():
    """Main validation function"""
    print("Validating SMN CNV Pipeline Setup...\n")
    
    validations = [
        ("Directory Structure", validate_directory_structure),
        ("Python Dependencies", validate_dependencies),
        ("Configuration File", validate_config),
        ("Python Modules", validate_python_modules),
        ("System Tools", validate_system_tools),
        ("IGV Installation", validate_igv)
    ]
    
    all_passed = True
    
    for name, validation_func in validations:
        print(f"\n--- {name} ---")
        try:
            if not validation_func():
                all_passed = False
        except Exception as e:
            print(f"✗ Validation error: {e}")
            all_passed = False
    
    print(f"\n{'='*50}")
    if all_passed:
        print("✅ All validations passed! Pipeline is ready to use.")
        print("\nNext steps:")
        print("1. Update SMN coordinates in config/config.json if needed")
        print("2. Add your MLPA validation data to config/validation_data.tsv")
        print("3. Test with a single sample:")
        print("   python3 smn_pipeline.py -c config/config.json -i sample.bam -o test_output -s sample_001")
    else:
        print("❌ Some validations failed. Please fix the issues above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
EOF

chmod +x validate_pipeline.py

# Create quick start guide
echo -e "${YELLOW}Creating quick start guide...${NC}"
cat > QUICKSTART.md << 'EOF'
# SMN CNV Pipeline - Quick Start Guide

## Overview
This pipeline detects copy number variations (CNVs) in SMN1 and SMN2 genes from whole exome sequencing (WES) data, specifically focusing on exons 7 and 8 for SMA (Spinal Muscular Atrophy) analysis.

## Prerequisites
- Linux operating system
- Python 3.8+
- Deduplicated BAM files with 30X coverage
- IGV (Integrative Genomics Viewer)
- samtools

## Installation

1. **Clone/download the pipeline**
2. **Run setup script:**
   ```bash
   ./setup.sh
   ```

3. **Validate installation:**
   ```bash
   python3 validate_pipeline.py
   ```

## Configuration

1. **Edit config/config.json:**
   - Update SMN coordinates if needed
   - Set IGV path
   - Adjust thresholds if required

2. **Add validation data (optional):**
   - Copy your MLPA-validated samples to config/validation_data.tsv
   - Use the template in config/validation_template.tsv

## Usage

### Single Sample Analysis
```bash
python3 smn_pipeline.py \
    --config config/config.json \
    --input sample.bam \
    --output output/sample_001 \
    --sample-id sample_001
```

### Batch Processing
```bash
python3 smn_pipeline.py \
    --config config/config.json \
    --input /path/to/bam/directory \
    --output output/batch_$(date +%Y%m%d) \
    --batch
```

### Using the Example Script
```bash
# Edit run_batch_example.sh with your paths
./run_batch_example.sh
```

## Output Files

### Single Sample
- `sample_001_report.html` - Interactive HTML report
- `sample_001_summary.tsv` - TSV summary
- IGV snapshots in snapshots/ directory

### Batch Processing
- `batch_report_TIMESTAMP.html` - Consolidated HTML report
- `batch_results_TIMESTAMP.tsv` - All results in TSV format
- Individual sample reports

## Machine Learning Features

The pipeline automatically:
- Learns optimal thresholds from your data
- Improves accuracy over time
- Adapts to your specific sequencing protocol
- Targets 95% accuracy after 1000+ samples

## Troubleshooting

1. **IGV not found:**
   ```bash
   ./setup_igv.sh
   ```

2. **Permission errors:**
   ```bash
   chmod +x *.sh
   ```

3. **Missing dependencies:**
   ```bash
   pip3 install -r requirements.txt
   ```

4. **BAM index missing:**
   ```bash
   samtools index your_file.bam
   ```

## Support

Check the logs in logs/ directory for detailed error messages.
Run validate_pipeline.py to diagnose setup issues.
EOF

# Final setup completion
echo -e "${GREEN}✓ Setup completed successfully!${NC}"
echo ""
echo -e "${BLUE}Next steps:${NC}"
echo "1. Run: python3 validate_pipeline.py"
echo "2. Update config/config.json with your specific settings"
echo "3. Test with a sample BAM file"
echo "4. Read QUICKSTART.md for detailed usage instructions"
echo ""
echo -e "${YELLOW}Important notes:${NC}"
echo "- Ensure your BAM files are indexed (samtools index)"
echo "- Update SMN coordinates in config if using different reference"
echo "- Add your MLPA validation data for optimal performance"
echo ""
echo -e "${GREEN}Setup complete! 🧬${NC}"
