#!/bin/bash

# IGV Setup Helper Script for SMN CNV Pipeline
# Automatically downloads and configures IGV for headless operation

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
IGV_VERSION="2.16.2"
IGV_BASE_URL="https://data.broadinstitute.org/igv/projects/downloads"
IGV_INSTALL_DIR="$HOME/igv"
CONFIG_FILE="config/config.json"

echo -e "${BLUE}=== IGV Setup Helper for SMN Pipeline ===${NC}"
echo "This script will download and configure IGV for automated screenshot generation"

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

# Check if IGV is already installed and working
check_existing_igv() {
    log_message "INFO" "Checking for existing IGV installation..."
    
    # Check common IGV locations
    IGV_PATHS=(
        "igv.sh"
        "igv"
        "/usr/local/bin/igv.sh"
        "/usr/local/igv/igv.sh"
        "/opt/igv/igv.sh"
        "$HOME/IGV_*/igv.sh"
        "$HOME/igv/IGV_*/igv.sh"
    )
    
    for igv_path in "${IGV_PATHS[@]}"; do
        if command -v "$igv_path" &> /dev/null || [ -f "$igv_path" ]; then
            log_message "INFO" "Found IGV at: $igv_path"
            
            # Test IGV
            if test_igv "$igv_path"; then
                log_message "SUCCESS" "IGV is already installed and working!"
                
                # Update config file
                update_config_file "$igv_path"
                
                echo ""
                echo -e "${GREEN}✅ IGV setup complete!${NC}"
                echo "IGV path: $igv_path"
                echo "Configuration updated in: $CONFIG_FILE"
                exit 0
            fi
        fi
    done
    
    log_message "WARN" "No working IGV installation found"
}

# Test IGV functionality
test_igv() {
    local igv_path=$1
    
    log_message "INFO" "Testing IGV at: $igv_path"
    
    # Try to get IGV version/help
    if timeout 10s "$igv_path" --help &>/dev/null; then
        log_message "SUCCESS" "IGV test passed"
        return 0
    elif timeout 10s "$igv_path" -h &>/dev/null; then
        log_message "SUCCESS" "IGV test passed"
        return 0
    else
        log_message "WARN" "IGV test failed for: $igv_path"
        return 1
    fi
}

# Download and install IGV
install_igv() {
    log_message "INFO" "Starting IGV installation..."
    
    # Create installation directory
    mkdir -p "$IGV_INSTALL_DIR"
    cd "$IGV_INSTALL_DIR"
    
    # Determine download URL based on system
    local igv_package=""
    local download_url=""
    
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        igv_package="IGV_Linux_${IGV_VERSION}_WithJava.zip"
        download_url="${IGV_BASE_URL}/${IGV_VERSION}/${igv_package}"
    else
        log_message "ERROR" "Unsupported operating system: $OSTYPE"
        log_message "ERROR" "This script supports Linux only"
        exit 1
    fi
    
    log_message "INFO" "Downloading IGV from: $download_url"
    
    # Download IGV
    if command -v wget &> /dev/null; then
        if ! wget -O "$igv_package" "$download_url"; then
            log_message "ERROR" "Failed to download IGV with wget"
            try_alternative_download "$igv_package" "$download_url"
        fi
    elif command -v curl &> /dev/null; then
        if ! curl -L -o "$igv_package" "$download_url"; then
            log_message "ERROR" "Failed to download IGV with curl"
            try_alternative_download "$igv_package" "$download_url"
        fi
    else
        log_message "ERROR" "Neither wget nor curl found. Please install one of them."
        exit 1
    fi
    
    # Verify download
    if [ ! -f "$igv_package" ] || [ ! -s "$igv_package" ]; then
        log_message "ERROR" "IGV download failed or file is empty"
        exit 1
    fi
    
    log_message "SUCCESS" "IGV downloaded successfully"
    
    # Extract IGV
    log_message "INFO" "Extracting IGV..."
    
    if ! unzip -q "$igv_package"; then
        log_message "ERROR" "Failed to extract IGV package"
        exit 1
    fi
    
    # Remove zip file
    rm "$igv_package"
    
    # Find IGV executable
    local igv_exec=$(find . -name "igv.sh" -type f | head -1)
    if [ -z "$igv_exec" ]; then
        log_message "ERROR" "Could not find igv.sh after extraction"
        exit 1
    fi
    
    # Make executable
    chmod +x "$igv_exec"
    
    # Get full path
    local igv_full_path="$(cd "$(dirname "$igv_exec")" && pwd)/$(basename "$igv_exec")"
    
    log_message "SUCCESS" "IGV installed at: $igv_full_path"
    
    # Test installation
    if test_igv "$igv_full_path"; then
        log_message "SUCCESS" "IGV installation verified!"
        
        # Update configuration
        update_config_file "$igv_full_path"
        
        # Create convenience symlink
        create_symlink "$igv_full_path"
        
        return 0
    else
        log_message "ERROR" "IGV installation verification failed"
        return 1
    fi
}

# Try alternative download methods
try_alternative_download() {
    local package=$1
    local url=$2
    
    log_message "WARN" "Trying alternative download methods..."
    
    # Try different versions
    local alt_versions=("2.16.1" "2.15.4" "2.14.1")
    
    for version in "${alt_versions[@]}"; do
        local alt_package="IGV_Linux_${version}_WithJava.zip"
        local alt_url="${IGV_BASE_URL}/${version}/${alt_package}"
        
        log_message "INFO" "Trying IGV version $version..."
        
        if command -v wget &> /dev/null; then
            if wget -O "$package" "$alt_url"; then
                log_message "SUCCESS" "Downloaded IGV version $version"
                return 0
            fi
        elif command -v curl &> /dev/null; then
            if curl -L -o "$package" "$alt_url"; then
                log_message "SUCCESS" "Downloaded IGV version $version"
                return 0
            fi
        fi
    done
    
    log_message "ERROR" "All download attempts failed"
    log_message "ERROR" "Please download IGV manually from: https://software.broadinstitute.org/software/igv/download"
    exit 1
}

# Update configuration file
update_config_file() {
    local igv_path=$1
    
    log_message "INFO" "Updating configuration file: $CONFIG_FILE"
    
    if [ ! -f "$CONFIG_FILE" ]; then
        log_message "WARN" "Configuration file not found: $CONFIG_FILE"
        log_message "INFO" "Please manually set igv_path to: $igv_path"
        return
    fi
    
    # Create backup
    cp "$CONFIG_FILE" "${CONFIG_FILE}.backup.$(date +%Y%m%d_%H%M%S)"
    
    # Update IGV path using Python (more reliable than sed for JSON)
    python3 << EOF
import json
import sys

try:
    with open('$CONFIG_FILE', 'r') as f:
        config = json.load(f)
    
    # Update IGV path
    if 'igv_settings' not in config:
        config['igv_settings'] = {}
    
    config['igv_settings']['igv_path'] = '$igv_path'
    
    with open('$CONFIG_FILE', 'w') as f:
        json.dump(config, f, indent=2)
    
    print("Configuration updated successfully")

except Exception as e:
    print(f"Error updating configuration: {e}")
    sys.exit(1)
EOF
    
    if [ $? -eq 0 ]; then
        log_message "SUCCESS" "Configuration file updated with IGV path"
    else
        log_message "ERROR" "Failed to update configuration file"
        log_message "INFO" "Please manually set igv_path to: $igv_path in $CONFIG_FILE"
    fi
}

# Create convenience symlink
create_symlink() {
    local igv_path=$1
    local link_path="$HOME/bin/igv.sh"
    
    # Create ~/bin if it doesn't exist
    mkdir -p "$HOME/bin"
    
    # Create symlink
    if ln -sf "$igv_path" "$link_path" 2>/dev/null; then
        log_message "SUCCESS" "Created symlink: $link_path"
        
        # Check if ~/bin is in PATH
        if [[ ":$PATH:" != *":$HOME/bin:"* ]]; then
            log_message "INFO" "To use 'igv.sh' from anywhere, add this to your ~/.bashrc:"
            echo "export PATH=\"\$HOME/bin:\$PATH\""
        fi
    else
        log_message "WARN" "Could not create symlink in $HOME/bin"
    fi
}

# Check Java installation
check_java() {
    log_message "INFO" "Checking Java installation..."
    
    if command -v java &> /dev/null; then
        local java_version=$(java -version 2>&1 | head -1)
        log_message "SUCCESS" "Java found: $java_version"
        return 0
    else
        log_message "WARN" "Java not found in PATH"
        log_message "INFO" "IGV includes Java, so this may not be an issue"
        return 1
    fi
}

# Check system requirements
check_requirements() {
    log_message "INFO" "Checking system requirements..."
    
    # Check operating system
    if [[ "$OSTYPE" != "linux-gnu"* ]]; then
        log_message "ERROR" "This script only supports Linux"
        exit 1
    fi
    
    # Check architecture
    local arch=$(uname -m)
    if [[ "$arch" != "x86_64" ]]; then
        log_message "WARN" "IGV is optimized for x86_64, detected: $arch"
    fi
    
    # Check disk space (need ~500MB for IGV)
    local available_space=$(df -BM . | tail -1 | awk '{print $4}' | sed 's/M//')
    if [ "$available_space" -lt 500 ]; then
        log_message "WARN" "Low disk space: ${available_space}MB available, 500MB+ recommended"
    fi
    
    # Check memory
    local total_mem=$(free -m | awk '/^Mem:/{print $2}')
    if [ "$total_mem" -lt 2048 ]; then
        log_message "WARN" "Low memory: ${total_mem}MB available, 2GB+ recommended for IGV"
    fi
    
    log_message "SUCCESS" "System requirements check completed"
}

# Create IGV preferences for batch processing
create_igv_preferences() {
    local igv_dir=$(dirname "$1")
    local prefs_file="$igv_dir/igv_batch.preferences"
    
    log_message "INFO" "Creating IGV preferences for batch processing..."
    
    cat > "$prefs_file" << 'EOF'
# IGV Preferences for SMN Pipeline Batch Processing
DEFAULT_GENOME=hg38
SHOW_SINGLE_TRACK_PANE=true
NORMALIZE_COVERAGE=true
MAX_DEPTH=1000
SHOW_CENTER_LINE=true
BACKGROUND_COLOR=255,255,255
SHOW_REGION_BARS=true
SHOW_JUNCTION_TRACK=false
SHOW_COVERAGE_TRACK=true
SAM_SHOW_DUPLICATES=false
SAM_FILTER_ALIGNMENTS=true
SAM_QUALITY_THRESHOLD=0
SAM_FLAG_UNMAPPED_PAIR=false
OVERLAY_TRACKS_TEXT=true
CHART_DRAW_TRACK_NAME=true
ENABLE_GOOGLE_MENU=false
SEND_USAGE_STATS=false
EOF
    
    log_message "SUCCESS" "IGV preferences created: $prefs_file"
}

# Display usage instructions
show_usage_instructions() {
    local igv_path=$1
    
    echo ""
    echo -e "${GREEN}🎉 IGV Setup Complete!${NC}"
    echo ""
    echo "=== IGV Installation Details ==="
    echo "IGV Path: $igv_path"
    echo "Version: IGV $IGV_VERSION"
    echo "Installation Directory: $IGV_INSTALL_DIR"
    echo ""
    echo "=== Testing IGV ==="
    echo "You can test IGV manually with:"
    echo "  $igv_path"
    echo ""
    echo "=== SMN Pipeline Integration ==="
    echo "IGV is now configured for the SMN pipeline."
    echo "The pipeline will automatically use IGV for screenshots."
    echo ""
    echo "=== Next Steps ==="
    echo "1. Test the pipeline with IGV:"
    echo "   python3 smn_pipeline.py -c config/config.json -i sample.bam -o test_output -s test"
    echo ""
    echo "2. Check IGV snapshots are generated in the output directory"
    echo ""
    echo "3. If you encounter issues, check:"
    echo "   - IGV memory settings in config/config.json"
    echo "   - Display/X11 forwarding if using SSH"
    echo "   - System memory availability"
    echo ""
    echo "=== Troubleshooting ==="
    echo "- Increase IGV memory: Edit 'igv_memory' in config/config.json"
    echo "- For headless systems: Ensure Xvfb is available or disable IGV"
    echo "- View logs: Check logs/ directory for IGV-related errors"
}

# Main execution
main() {
    echo "Starting IGV setup for SMN CNV Pipeline..."
    echo ""
    
    # Check system requirements
    check_requirements
    
    # Check for existing IGV
    check_existing_igv
    
    # Check Java (informational)
    check_java
    
    # Ask user confirmation
    echo ""
    echo "IGV will be downloaded and installed to: $IGV_INSTALL_DIR"
    echo "This will download ~200MB and require ~500MB disk space."
    read -p "Continue with IGV installation? (y/N): " confirm
    
    if [[ ! $confirm =~ ^[Yy]$ ]]; then
        echo "Installation cancelled"
        exit 0
    fi
    
    # Install IGV
    if install_igv; then
        local installed_igv=$(find "$IGV_INSTALL_DIR" -name "igv.sh" -type f | head -1)
        if [ -n "$installed_igv" ]; then
            # Get absolute path
            local abs_igv_path="$(cd "$(dirname "$installed_igv")" && pwd)/$(basename "$installed_igv")"
            
            # Create preferences
            create_igv_preferences "$abs_igv_path"
            
            # Show usage instructions
            show_usage_instructions "$abs_igv_path"
        fi
    else
        log_message "ERROR" "IGV installation failed"
        exit 1
    fi
}

# Handle command line arguments
case "$1" in
    --help|-h)
        cat << EOF
IGV Setup Helper for SMN CNV Pipeline

USAGE:
    $0 [OPTIONS]

OPTIONS:
    --help, -h          Show this help message
    --check            Check for existing IGV installation only
    --test PATH        Test IGV installation at specified path
    --uninstall        Remove IGV installation

EXAMPLES:
    $0                 # Install IGV automatically
    $0 --check         # Check for existing IGV
    $0 --test /usr/local/bin/igv.sh  # Test specific IGV installation

EOF
        exit 0
        ;;
    --check)
        check_existing_igv
        exit 0
        ;;
    --test)
        if [ -z "$2" ]; then
            echo "Error: --test requires IGV path"
            exit 1
        fi
        test_igv "$2"
        exit $?
        ;;
    --uninstall)
        if [ -d "$IGV_INSTALL_DIR" ]; then
            log_message "INFO" "Removing IGV installation: $IGV_INSTALL_DIR"
            rm -rf "$IGV_INSTALL_DIR"
            log_message "SUCCESS" "IGV uninstalled"
        else
            log_message "INFO" "No IGV installation found to remove"
        fi
        exit 0
        ;;
    "")
        # No arguments, run main installation
        main
        ;;
    *)
        echo "Unknown option: $1"
        echo "Use --help for usage information"
        exit 1
        ;;
esac
