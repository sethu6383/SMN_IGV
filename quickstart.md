# SMN CNV Pipeline - Quick Start Guide

Get up and running with SMN CNV detection in minutes!

## 🚀 Quick Setup (5 minutes)

### 1. Prerequisites Check
```bash
# Ensure you have:
# - Linux OS (Ubuntu 18.04+, CentOS 7+)
# - Python 3.8+
# - 8GB+ RAM
# - 10GB+ free disk space

python3 --version  # Should be 3.8+
```

### 2. Get the Pipeline
```bash
# Clone or download the repository
cd /your/working/directory
# [Extract/clone the SMN pipeline files here]

# Make scripts executable
chmod +x *.sh
chmod +x validate_pipeline.py
```

### 3. Automatic Setup
```bash
# Run the setup script (handles everything automatically)
./setup.sh

# Validate the installation
python3 validate_pipeline.py
```

**✅ If validation passes, you're ready to go!**

## 🧬 First Analysis (2 minutes)

### Test with a Single Sample
```bash
# Basic command structure
python3 smn_pipeline.py \
    --config config/config.json \
    --input your_sample.bam \
    --output results/test_run \
    --sample-id TEST_001

# Example with real paths
python3 smn_pipeline.py \
    -c config/config.json \
    -i /data/samples/patient_001.bam \
    -o output/patient_001 \
    -s patient_001
```

### Check Results
```bash
# View HTML report
firefox output/patient_001/patient_001_report.html

# Check TSV summary
cat output/patient_001/patient_001_summary.tsv
```

## 🔄 Batch Processing (5 minutes)

### Process Multiple Samples
```bash
# Using the advanced batch script
./run_batch_advanced.sh \
    --input /path/to/bam/directory \
    --output results/batch_$(date +%Y%m%d) \
    --batch-size 20 \
    --parallel 4

# Or using the main pipeline script
python3 smn_pipeline.py \
    --config config/config.json \
    --input /path/to/bam/directory \
    --output results/batch_analysis \
    --batch
```

### Monitor Progress
```bash
# Watch log files
tail -f output/batch_analysis/logs/smn_pipeline_*.log

# Check resource usage
htop
```

## 📊 Understanding Results

### Result Files Structure
```
output/
├── batch_report_20241201_143022.html    # 📈 Interactive batch report
├── batch_results_20241201_143022.tsv    # 📋 All results in table format
├── individual_reports/                   # 📁 Per-sample reports
│   ├── sample_001_report.html
│   ├── sample_002_report.html
│   └── ...
├── snapshots/                           # 🖼️ IGV screenshots
│   ├── sample_001_smn1_exon7.png
│   └── ...
└── logs/                               # 📝 Processing logs
```

### Key Result Columns (TSV)
- **`sma_risk`**: `HIGH_RISK`, `CARRIER`, `LOW_RISK`, or `UNCERTAIN`
- **`smn1_exon7_call`**: CNV call (`HOMO_DEL`, `HETERO_DEL`, `NORMAL`, `DUP`)
- **`smn1_exon7_confidence`**: Confidence score (0.0-1.0)
- **`depth_smn1_exon7`**: Mean sequencing depth

### Clinical Interpretation

| SMA Risk | Clinical Meaning | Action Required |
|----------|------------------|-----------------|
| **HIGH_RISK** | 🚨 Likely SMA affected | Urgent referral to genetics |
| **CARRIER** | ⚠️ SMA carrier | Genetic counseling recommended |
| **LOW_RISK** | ✅ Normal SMN1 | Routine reporting |
| **UNCERTAIN** | ❓ Unclear result | Additional testing needed |

## ⚙️ Configuration

### Essential Config Updates

Edit `config/config.json`:

```json
{
  "smn_coordinates": {
    "smn1_exon7": {"chr": "chr5", "start": 70247724, "end": 70247775},
    // Update these if using different reference
  },
  "igv_settings": {
    "igv_path": "/usr/local/bin/igv.sh",  // Update IGV path
    "igv_memory": "4g"
  },
  "validation_data_file": "config/my_mlpa_data.tsv"  // Add your MLPA data
}
```

### Add Validation Data (Recommended)

Create `config/my_mlpa_data.tsv`:
```tsv
sample_id	smn1_exon7	smn1_exon8	smn2_exon7	smn2_exon8	comments
CONTROL_001	NORMAL	NORMAL	NORMAL	NORMAL	Normal control
CARRIER_001	HETERO_DEL	HETERO_DEL	NORMAL	NORMAL	Known carrier
SMA_001	HOMO_DEL	HOMO_DEL	NORMAL	NORMAL	SMA patient
```

## 🔧 Common Issues & Quick Fixes

### ❌ "IGV not found"
```bash
# Install IGV automatically
./setup_igv.sh

# Or update config manually
nano config/config.json
# Set "igv_path": "/full/path/to/igv.sh"
```

### ❌ "BAM file not indexed"
```bash
# Index your BAM files
samtools index your_file.bam

# Or index all BAM files in directory
for bam in /path/to/bams/*.bam; do
    samtools index "$bam"
done
```

### ❌ "Python dependencies missing"
```bash
# Install missing packages
pip3 install -r requirements.txt

# Or install individual packages
pip3 install pysam numpy pandas scipy scikit-learn plotly optuna
```

### ❌ "Permission denied"
```bash
# Fix permissions
chmod +x *.sh
chmod +x *.py
```

### ❌ "Memory error"
```bash
# Reduce IGV memory in config/config.json
"igv_memory": "2g"

# Process smaller batches
./run_batch_advanced.sh --batch-size 10
```

## 📋 Quick Command Reference

### Single Sample Commands
```bash
# Minimal command
python3 smn_pipeline.py -c config/config.json -i sample.bam -o output -s sample_id

# With all options
python3 smn_pipeline.py \
    --config config/config.json \
    --input sample.bam \
    --output output/sample_001 \
    --sample-id sample_001

# Dry run (check without processing)
python3 smn_pipeline.py --config config/config.json --input sample.bam --output test --dry-run
```

### Batch Processing Commands
```bash
# Basic batch
python3 smn_pipeline.py -c config/config.json -i /bam/dir -o output --batch

# Advanced batch with options
./run_batch_advanced.sh -i /bam/dir -o results -b 30 -p 8 --qc-report

# Dry run for batch
./run_batch_advanced.sh -i /bam/dir -o test --dry-run
```

### Utility Commands
```bash
# Validate setup
python3 validate_pipeline.py

# Run tests
./run_tests.sh

# Check system requirements
python3 -c "from src.utils import check_system_requirements; print(check_system_requirements())"
```

## 🎯 Performance Tips

### For Best Results
1. **Use indexed BAM files** - 10x faster processing
2. **Process in batches of 20-30** - optimal resource usage
3. **Add MLPA validation data** - improves accuracy over time
4. **Use SSD storage** - faster I/O operations
5. **Ensure 30X+ coverage** - better CNV detection

### Resource Optimization
```bash
# For high-memory systems
./run_batch_advanced.sh -i /data -o results -b 50 -p 8

# For limited-memory systems
./run_batch_advanced.sh -i /data -o results -b 10 -p 2

# Monitor resources
htop  # Check CPU/memory usage
df -h  # Check disk space
```

## 🚀 Next Steps

### After First Successful Run
1. **Review HTML reports** - understand your results
2. **Add MLPA validation data** - improve accuracy
3. **Process larger batches** - leverage ML improvements
4. **Set up routine processing** - automate your workflow

### For Production Use
1. **Validate with known samples** - confirm accuracy
2. **Optimize batch sizes** - match your hardware
3. **Set up monitoring** - track performance over time
4. **Plan data management** - archive old results

### Advanced Features
```bash
# Force threshold update after new data
python3 smn_pipeline.py --config config/config.json --input data --output results --update-thresholds

# Generate performance report
python3 -c "
from src.ml_threshold import MLThresholdOptimizer
optimizer = MLThresholdOptimizer({'model_dir': 'models'})
print(optimizer.export_model_performance())
"

# Reset learning for fresh start
python3 -c "
from src.ml_threshold import MLThresholdOptimizer
optimizer = MLThresholdOptimizer({'model_dir': 'models'})
optimizer.reset_learning()
"
```

## 📞 Getting Help

### Debug Information
```bash
# Check pipeline status
python3 validate_pipeline.py

# View recent logs
tail -n 50 logs/smn_pipeline_*.log

# Check dependencies
python3 -c "from src import check_dependencies; print(check_dependencies())"
```

### Common Workflows

**Daily Clinical Processing:**
```bash
./run_batch_advanced.sh -i /daily/samples -o /results/$(date +%Y%m%d) -b 20 --qc-report
```

**Research Cohort Analysis:**
```bash
./run_batch_advanced.sh -i /research/cohort -o /analysis/cohort_study -b 50 --update-thresholds
```

**Quality Check Single Sample:**
```bash
python3 smn_pipeline.py -c config/config.json -i suspicious.bam -o qc_check -s suspicious_sample
```

---

**🎉 You're now ready to detect SMN CNVs efficiently!**

For detailed documentation, see [README.md](README.md)  
For troubleshooting, check the logs in `logs/` directory  
For advanced configuration, explore `config/config.json`
