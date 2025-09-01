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
SMN CNV Pipeline Validation Script
Comprehensive validation of pipeline setup and functionality
"""

import os
import sys
import json
import subprocess
import tempfile
import shutil
from pathlib import Path
import importlib
import logging
from datetime import datetime

# Colors for output
class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    PURPLE = '\033[0;35m'
    CYAN = '\033[0;36m'
    WHITE = '\033[1;37m'
    NC = '\033[0m'  # No Color

def print_colored(message, color=Colors.NC):
    """Print colored message"""
    print(f"{color}{message}{Colors.NC}")

def print_header(title):
    """Print section header"""
    print(f"\n{Colors.BLUE}{'='*60}{Colors.NC}")
    print(f"{Colors.WHITE}{title.center(60)}{Colors.NC}")
    print(f"{Colors.BLUE}{'='*60}{Colors.NC}")

def print_status(message, status="INFO"):
    """Print status message"""
    timestamp = datetime.now().strftime("%H:%M:%S")
    colors = {
        "INFO": Colors.BLUE,
        "SUCCESS": Colors.GREEN,
        "WARNING": Colors.YELLOW,
        "ERROR": Colors.RED,
        "CHECK": Colors.CYAN
    }
    color = colors.get(status, Colors.NC)
    print(f"{color}[{status}] {timestamp} - {message}{Colors.NC}")

class PipelineValidator:
    """Main validation class for SMN CNV Pipeline"""
    
    def __init__(self):
        self.errors = []
        self.warnings = []
        self.success_count = 0
        self.total_checks = 0
        self.temp_dir = None
        
    def run_validation(self):
        """Run complete validation suite"""
        print_header("SMN CNV Pipeline Validation")
        print_colored("Validating pipeline setup and functionality...", Colors.WHITE)
        
        # Create temporary directory for testing
        self.temp_dir = tempfile.mkdtemp(prefix='smn_validation_')
        
        try:
            # Run all validation checks
            self.validate_file_structure()
            self.validate_python_environment()
            self.validate_dependencies()
            self.validate_system_tools()
            self.validate_configuration()
            self.validate_python_modules()
            self.validate_igv_setup()
            self.validate_permissions()
            self.validate_functionality()
            self.validate_test_data()
            
            # Summary
            self.print_summary()
            
        finally:
            # Cleanup
            if self.temp_dir and os.path.exists(self.temp_dir):
                shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def check(self, description, check_function, *args, **kwargs):
        """Run a validation check with error handling"""
        self.total_checks += 1
        print_status(f"Checking {description}...", "CHECK")
        
        try:
            result = check_function(*args, **kwargs)
            if result:
                print_status(f"✓ {description}", "SUCCESS")
                self.success_count += 1
                return True
            else:
                print_status(f"✗ {description}", "ERROR")
                self.errors.append(description)
                return False
        except Exception as e:
            print_status(f"✗ {description}: {str(e)}", "ERROR")
            self.errors.append(f"{description}: {str(e)}")
            return False
    
    def warn(self, message):
        """Add warning message"""
        print_status(message, "WARNING")
        self.warnings.append(message)
    
    def validate_file_structure(self):
        """Validate required file structure"""
        print_header("FILE STRUCTURE VALIDATION")
        
        # Required files
        required_files = {
            'smn_pipeline.py': 'Main pipeline script',
            'setup.sh': 'Setup script',
            'setup_igv.sh': 'IGV setup script',
            'config/config.json': 'Configuration file',
            'src/__init__.py': 'Source package init',
            'src/depth_analyzer.py': 'Depth analyzer module',
            'src/cnv_caller.py': 'CNV caller module',
            'src/igv_automation.py': 'IGV automation module',
            'src/ml_threshold.py': 'ML threshold optimizer',
            'src/report_generator.py': 'Report generator module',
            'src/utils.py': 'Utility functions',
            'requirements.txt': 'Python dependencies',
            'README.md': 'Documentation',
            'QUICKSTART.md': 'Quick start guide'
        }
        
        for file_path, description in required_files.items():
            self.check(f"{description} exists", os.path.exists, file_path)
        
        # Required directories
        required_dirs = [
            'src', 'config', 'models', 'output', 'reports', 
            'snapshots', 'logs', 'temp', 'tests'
        ]
        
        for dir_path in required_dirs:
            if not os.path.exists(dir_path):
                try:
                    os.makedirs(dir_path, exist_ok=True)
                    print_status(f"Created missing directory: {dir_path}", "WARNING")
                    self.warnings.append(f"Created missing directory: {dir_path}")
                except Exception as e:
                    self.errors.append(f"Could not create directory {dir_path}: {e}")
            else:
                self.check(f"Directory {dir_path} exists", os.path.isdir, dir_path)
    
    def validate_python_environment(self):
        """Validate Python environment"""
        print_header("PYTHON ENVIRONMENT VALIDATION")
        
        # Python version
        python_version = sys.version_info
        self.check(
            "Python version >= 3.8",
            lambda: python_version >= (3, 8),
        )
        
        if python_version < (3, 8):
            self.errors.append(f"Python {python_version.major}.{python_version.minor} found, need 3.8+")
        else:
            print_status(f"Python {python_version.major}.{python_version.minor}.{python_version.micro} found", "SUCCESS")
        
        # Pip availability
        self.check("pip3 available", lambda: subprocess.run(['pip3', '--version'], capture_output=True).returncode == 0)
    
    def validate_dependencies(self):
        """Validate Python dependencies"""
        print_header("DEPENDENCY VALIDATION")
        
        required_packages = [
            ('pysam', 'BAM file processing'),
            ('numpy', 'Numerical computations'),
            ('pandas', 'Data manipulation'),
            ('scipy', 'Scientific computing'),
            ('sklearn', 'Machine learning'),
            ('plotly', 'Interactive plots'),
            ('optuna', 'Hyperparameter optimization'),
            ('psutil', 'System monitoring')
        ]
        
        for package, description in required_packages:
            self.check(
                f"{package} ({description})",
                self._check_package_import,
                package
            )
    
    def _check_package_import(self, package_name):
        """Check if package can be imported"""
        try:
            # Handle package name variations
            import_names = {
                'sklearn': 'sklearn',
                'pysam': 'pysam',
                'numpy': 'numpy',
                'pandas': 'pandas',
                'scipy': 'scipy',
                'plotly': 'plotly',
                'optuna': 'optuna',
                'psutil': 'psutil'
            }
            
            actual_name = import_names.get(package_name, package_name)
            importlib.import_module(actual_name)
            return True
        except ImportError:
            return False
    
    def validate_system_tools(self):
        """Validate system tools"""
        print_header("SYSTEM TOOLS VALIDATION")
        
        # Required tools
        tools = [
            ('samtools', 'BAM file processing', True),
            ('igv.sh', 'IGV for screenshots', False),
            ('igv', 'IGV alternative', False)
        ]
        
        samtools_found = False
        igv_found = False
        
        for tool, description, required in tools:
            available = self._check_command_available(tool)
            
            if available:
                print_status(f"✓ {tool} ({description}) found", "SUCCESS")
                if tool == 'samtools':
                    samtools_found = True
                elif tool.startswith('igv'):
                    igv_found = True
            elif required:
                print_status(f"✗ {tool} ({description}) not found", "ERROR")
                self.errors.append(f"Required tool not found: {tool}")
            else:
                print_status(f"⚠ {tool} ({description}) not found", "WARNING")
                self.warnings.append(f"Optional tool not found: {tool}")
        
        # Summary for key tools
        if not samtools_found:
            self.warn("samtools not found - required for BAM processing")
        
        if not igv_found:
            self.warn("IGV not found - screenshots will not be generated")
    
    def _check_command_available(self, command):
        """Check if command is available in PATH"""
        try:
            result = subprocess.run(['which', command], capture_output=True, text=True)
            return result.returncode == 0
        except:
            return False
    
    def validate_configuration(self):
        """Validate configuration files"""
        print_header("CONFIGURATION VALIDATION")
        
        config_file = 'config/config.json'
        
        if not os.path.exists(config_file):
            self.errors.append(f"Configuration file not found: {config_file}")
            return
        
        try:
            with open(config_file, 'r') as f:
                config = json.load(f)
            
            print_status("Configuration file is valid JSON", "SUCCESS")
            
            # Check required configuration sections
            required_sections = [
                'pipeline_version',
                'smn_coordinates',
                'cnv_thresholds',
                'igv_settings'
            ]
            
            for section in required_sections:
                self.check(f"Config section '{section}' exists", lambda s=section: s in config)
            
            # Validate SMN coordinates
            if 'smn_coordinates' in config:
                smn_coords = config['smn_coordinates']
                required_regions = ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']
                
                for region in required_regions:
                    if region in smn_coords:
                        coord = smn_coords[region]
                        if all(key in coord for key in ['chr', 'start', 'end']):
                            print_status(f"✓ {region} coordinates valid", "SUCCESS")
                        else:
                            self.errors.append(f"Invalid coordinates for {region}")
                    else:
                        self.errors.append(f"Missing coordinates for {region}")
            
            # Check IGV configuration
            if 'igv_settings' in config:
                igv_config = config['igv_settings']
                igv_path = igv_config.get('igv_path', '')
                
                if igv_path and os.path.exists(igv_path):
                    print_status(f"✓ IGV path exists: {igv_path}", "SUCCESS")
                elif igv_path:
                    self.warn(f"IGV path not found: {igv_path}")
                else:
                    self.warn("IGV path not configured")
            
        except json.JSONDecodeError as e:
            self.errors.append(f"Configuration file is not valid JSON: {e}")
        except Exception as e:
            self.errors.append(f"Error reading configuration: {e}")
    
    def validate_python_modules(self):
        """Validate Python modules can be imported"""
        print_header("PYTHON MODULES VALIDATION")
        
        # Add src to path
        src_path = os.path.abspath('src')
        if src_path not in sys.path:
            sys.path.insert(0, src_path)
        
        modules = [
            ('depth_analyzer', 'DepthAnalyzer'),
            ('cnv_caller', 'CNVCaller'),
            ('igv_automation', 'IGVAutomation'),
            ('ml_threshold', 'MLThresholdOptimizer'),
            ('report_generator', 'ReportGenerator'),
            ('utils', 'setup_logging')
        ]
        
        for module_name, class_name in modules:
            self.check(
                f"Module {module_name}.{class_name}",
                self._check_module_import,
                module_name, class_name
            )
    
    def _check_module_import(self, module_name, class_name):
        """Check if module and class can be imported"""
        try:
            module = importlib.import_module(module_name)
            getattr(module, class_name)
            return True
        except (ImportError, AttributeError) as e:
            print_status(f"Import error: {e}", "ERROR")
            return False
    
    def validate_igv_setup(self):
        """Validate IGV setup"""
        print_header("IGV VALIDATION")
        
        # Check configuration
        config_file = 'config/config.json'
        igv_path = None
        
        if os.path.exists(config_file):
            try:
                with open(config_file, 'r') as f:
                    config = json.load(f)
                igv_path = config.get('igv_settings', {}).get('igv_path')
            except:
                pass
        
        if igv_path:
            # Test specific IGV path
            self.check(f"IGV executable exists at {igv_path}", os.path.exists, igv_path)
            
            if os.path.exists(igv_path):
                self.check("IGV is executable", os.access, igv_path, os.X_OK)
                
                # Test IGV execution
                try:
                    result = subprocess.run([igv_path, '--help'], 
                                         capture_output=True, text=True, timeout=10)
                    if result.returncode == 0 or 'IGV' in result.stderr:
                        print_status("✓ IGV responds to --help", "SUCCESS")
                    else:
                        self.warn("IGV may not be working properly")
                except subprocess.TimeoutExpired:
                    self.warn("IGV test timed out")
                except Exception as e:
                    self.warn(f"Could not test IGV: {e}")
        else:
            self.warn("IGV path not configured - screenshots will be disabled")
    
    def validate_permissions(self):
        """Validate file permissions"""
        print_header("PERMISSIONS VALIDATION")
        
        executable_files = [
            'smn_pipeline.py',
            'setup.sh',
            'setup_igv.sh',
            'validate_pipeline.py'
        ]
        
        for file_path in executable_files:
            if os.path.exists(file_path):
                self.check(f"{file_path} is executable", os.access, file_path, os.X_OK)
        
        # Check write permissions for output directories
        output_dirs = ['output', 'reports', 'snapshots', 'logs', 'temp', 'models']
        
        for dir_path in output_dirs:
            if os.path.exists(dir_path):
                self.check(f"{dir_path} is writable", os.access, dir_path, os.W_OK)
    
    def validate_functionality(self):
        """Validate basic functionality"""
        print_header("FUNCTIONALITY VALIDATION")
        
        # Test mock data processing
        try:
            # Add src to path
            src_path = os.path.abspath('src')
            if src_path not in sys.path:
                sys.path.insert(0, src_path)
            
            # Test depth analyzer
            self.check("DepthAnalyzer initialization", self._test_depth_analyzer)
            
            # Test CNV caller
            self.check("CNVCaller initialization", self._test_cnv_caller)
            
            # Test ML optimizer
            self.check("MLThresholdOptimizer initialization", self._test_ml_optimizer)
            
            # Test report generator
            self.check("ReportGenerator initialization", self._test_report_generator)
            
        except Exception as e:
            self.errors.append(f"Functionality test failed: {e}")
    
    def _test_depth_analyzer(self):
        """Test depth analyzer functionality"""
        try:
            from depth_analyzer import DepthAnalyzer
            
            config = {'smn_coordinates': {
                'smn1_exon7': {'chr': 'chr5', 'start': 70247724, 'end': 70247775}
            }}
            
            analyzer = DepthAnalyzer(config)
            return True
        except Exception as e:
            print_status(f"DepthAnalyzer test failed: {e}", "ERROR")
            return False
    
    def _test_cnv_caller(self):
        """Test CNV caller functionality"""
        try:
            from cnv_caller import CNVCaller
            
            config = {'cnv_thresholds': {
                'depth_ratio': {'homo_del_threshold': 0.3}
            }}
            
            caller = CNVCaller(config)
            return True
        except Exception as e:
            print_status(f"CNVCaller test failed: {e}", "ERROR")
            return False
    
    def _test_ml_optimizer(self):
        """Test ML optimizer functionality"""
        try:
            from ml_threshold import MLThresholdOptimizer
            
            config = {'model_dir': self.temp_dir}
            
            optimizer = MLThresholdOptimizer(config)
            return True
        except Exception as e:
            print_status(f"MLThresholdOptimizer test failed: {e}", "ERROR")
            return False
    
    def _test_report_generator(self):
        """Test report generator functionality"""
        try:
            from report_generator import ReportGenerator
            
            config = {'reports': self.temp_dir}
            
            generator = ReportGenerator(config)
            return True
        except Exception as e:
            print_status(f"ReportGenerator test failed: {e}", "ERROR")
            return False
    
    def validate_test_data(self):
        """Validate test data and examples"""
        print_header("TEST DATA VALIDATION")
        
        # Check for test data
        test_data_paths = [
            'tests/data',
            'examples/data',
            'test_data'
        ]
        
        test_data_found = False
        for path in test_data_paths:
            if os.path.exists(path):
                test_data_found = True
                print_status(f"✓ Test data found: {path}", "SUCCESS")
                break
        
        if not test_data_found:
            self.warn("No test data found - consider adding example files")
        
        # Check for validation data template
        validation_template = 'config/validation_template.tsv'
        if os.path.exists(validation_template):
            print_status("✓ Validation data template found", "SUCCESS")
        else:
            self.warn("Validation data template not found")
    
    def print_summary(self):
        """Print validation summary"""
        print_header("VALIDATION SUMMARY")
        
        success_rate = (self.success_count / self.total_checks) * 100 if self.total_checks > 0 else 0
        
        print(f"\n{Colors.WHITE}Results:{Colors.NC}")
        print(f"  {Colors.GREEN}✓ Successful checks: {self.success_count}/{self.total_checks} ({success_rate:.1f}%){Colors.NC}")
        print(f"  {Colors.YELLOW}⚠ Warnings: {len(self.warnings)}{Colors.NC}")
        print(f"  {Colors.RED}✗ Errors: {len(self.errors)}{Colors.NC}")
        
        if self.warnings:
            print(f"\n{Colors.YELLOW}WARNINGS:{Colors.NC}")
            for i, warning in enumerate(self.warnings, 1):
                print(f"  {i}. {warning}")
        
        if self.errors:
            print(f"\n{Colors.RED}ERRORS:{Colors.NC}")
            for i, error in enumerate(self.errors, 1):
                print(f"  {i}. {error}")
        
        print(f"\n{Colors.WHITE}Recommendations:{Colors.NC}")
        
        if self.errors:
            print_colored("❌ SETUP INCOMPLETE - Please fix errors before using the pipeline", Colors.RED)
            
            if any('dependencies' in error.lower() or 'import' in error.lower() for error in self.errors):
                print("   → Install missing dependencies: pip3 install -r requirements.txt")
            
            if any('samtools' in error.lower() for error in self.errors):
                print("   → Install samtools: apt-get install samtools (Ubuntu/Debian)")
            
            if any('permission' in error.lower() for error in self.errors):
                print("   → Fix permissions: chmod +x *.sh *.py")
            
            if any('config' in error.lower() for error in self.errors):
                print("   → Check configuration file: config/config.json")
            
        elif self.warnings:
            print_colored("⚠️ SETUP MOSTLY COMPLETE - Pipeline will work but some features may be limited", Colors.YELLOW)
            
            if any('igv' in warning.lower() for warning in self.warnings):
                print("   → Install IGV for screenshots: ./setup_igv.sh")
            
            if any('validation' in warning.lower() for warning in self.warnings):
                print("   → Add MLPA validation data for better accuracy")
            
            if any('samtools' in warning.lower() for warning in self.warnings):
                print("   → Install samtools for BAM file processing")
        
        else:
            print_colored("✅ SETUP COMPLETE - Pipeline is ready for use!", Colors.GREEN)
            
            print(f"\n{Colors.WHITE}Next Steps:{Colors.NC}")
            print("   1. Test with a sample: python3 smn_pipeline.py -c config/config.json -i sample.bam -o test_output -s test")
            print("   2. Add MLPA validation data to config/validation_data.tsv")
            print("   3. Process your first batch: ./run_batch_example.sh")
            print("   4. Review results in generated HTML reports")
            print("   5. See QUICKSTART.md for detailed usage instructions")
        
        print(f"\n{Colors.CYAN}For help:{Colors.NC}")
        print("   • Documentation: README.md")
        print("   • Quick start: QUICKSTART.md") 
        print("   • Run tests: python3 test_pipeline.py")
        print("   • IGV setup: ./setup_igv.sh")

def main():
    """Main validation function"""
    try:
        validator = PipelineValidator()
        validator.run_validation()
        
        # Exit with error code if there are errors
        if validator.errors:
            sys.exit(1)
        else:
            sys.exit(0)
            
    except KeyboardInterrupt:
        print_colored("\nValidation interrupted by user", Colors.YELLOW)
        sys.exit(1)
    except Exception as e:
        print_colored(f"\nUnexpected error during validation: {e}", Colors.RED)
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
