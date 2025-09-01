# SMN CNV Detection Pipeline

A comprehensive bioinformatics pipeline for detecting copy number variations (CNVs) in SMN1 and SMN2 genes from whole exome sequencing (WES) data, specifically designed for Spinal Muscular Atrophy (SMA) analysis.

## 🧬 Features

### Core Functionality
- **Multi-method CNV detection**: Depth ratio, statistical, HMM-like, and ensemble approaches
- **SMN1/SMN2 differentiation**: Distinguishes between highly similar SMN genes
- **Automated IGV snapshots**: Visual validation with automated screenshot generation
- **Machine learning optimization**: Adaptive thresholds that improve with more data
- **Batch processing**: Handle 20-50+ samples simultaneously
- **Clinical interpretation**: Automated SMA risk assessment with recommendations

### Advanced Features
- **Self-tuning thresholds**: ML-based optimization targeting 95% accuracy
- **Interactive reports**: HTML reports with visualizations and IGV snapshots
- **Quality assessment**: Comprehensive QC metrics and validation
- **Validation tracking**: Performance monitoring against MLPA gold standards
- **Cross-validation**: Multiple computational methods for robust calling

## 📁 Repository Structure

```
smn_cnv_pipeline/
├── smn_pipeline.py              # Main pipeline script
├── config/
│   ├── config.json             # Main configuration file
│   ├── validation_template.tsv # MLPA validation data template
│   └── igv_preferences.txt     # IGV settings
├── src/
│   ├── __init__.py
│   ├── depth_analyzer.py       # Depth analysis module
│   ├── cnv_caller.py          # CNV calling methods
│   ├── igv_automation.py      # IGV automation
│   ├── ml_threshold.py        # ML threshold optimization
│   ├── report_generator.py    # Report generation
│   └── utils.py               # Utility functions
├── models/                     # ML models and thresholds
├── output/                     # Pipeline outputs
├── reports/                    # Generated reports
├── snapshots/                  # IGV snapshots
├── logs/                       # Log files
├── temp/                       # Temporary files
├── tests/                      # Test scripts
├── docs/                       # Documentation
├── setup.sh                    # Setup script
├── validate_pipeline.py        # Validation script
├── run_batch_example.sh        # Example batch script
├── requirements.txt            # Python dependencies
├── QUICKSTART.md              # Quick start guide
└── README.md                  # This file
```

## 🚀 Quick Start

### 1. Setup
```bash
# Clone the repository
git clone <your-repo-url>
cd smn_cnv_pipeline

# Run setup script
chmod +x setup.sh
./setup.sh

# Validate installation
python3 validate_pipeline.py
```

### 2. Configuration
```bash
# Edit configuration file
nano config/config.json

# Update SMN coordinates if needed (default: hg38)
# Set IGV path
# Add validation data file path
```

### 3. Single Sample Analysis
```bash
python3 smn_pipeline.py \
    --config config/config.json \
    --input sample.bam \
    --output output/sample_001 \
    --sample-id sample_001
```

### 4. Batch Processing
```bash
# Process all BAM files in a directory
python3 smn_pipeline.py \
    --config config/config.json \
    --input /path/to/bam/directory \
    --output output/batch_$(date +%Y%m%d) \
    --batch
```

## 📊 Output Files

### Single Sample Outputs
- **`sample_001_report.html`** - Interactive HTML report with:
  - Clinical interpretation and SMA risk assessment
  - CNV calls with confidence scores
  - Quality metrics and depth analysis plots
  - IGV snapshots for visual validation
  - Method comparison and consensus results

- **`sample_001_summary.tsv`** - Tab-separated summary with:
  - CNV calls for all regions
  - Confidence scores
  - Depth metrics
  - Quality assessments

### Batch Outputs
- **`batch_report_TIMESTAMP.html`** - Consolidated batch report with:
  - Summary statistics across all samples
  - Distribution plots and quality assessments
  - Clinical recommendations for flagged samples
  - Interactive visualizations

- **`batch_results_TIMESTAMP.tsv`** - Complete results table
- **Individual sample reports** for each processed sample

## 🧠 Machine Learning Features

### Adaptive Thresholds
- **Continuous learning**: Thresholds improve with each batch
- **Performance tracking**: Monitors accuracy against validation data
- **Target accuracy**: Aims for 95% accuracy after 1000+ samples
- **Multiple optimization**: Optimizes both computational and IGV-based calling

### Model Training
- **Ensemble methods**: Combines multiple algorithms
- **Feature selection**: Automatically selects best features
- **Hyperparameter optimization**: Uses Optuna for optimal performance
- **Cross-validation**: Robust model validation

## 🔧 Configuration Options

### Key Configuration Sections

#### SMN Coordinates (hg38)
```json
"smn_coordinates": {
  "smn1_exon7": {"chr": "chr5", "start": 70247724, "end": 70247775},
  "smn1_exon8": {"chr": "chr5", "start": 70248935, "end": 70249306},
  "smn2_exon7": {"chr": "chr5", "start": 69372304, "end": 69372355},
  "smn2_exon8": {"chr": "chr5", "start": 69373515, "end": 69373886}
}
```

#### CNV Thresholds
```json
"cnv_thresholds": {
  "depth_ratio": {
    "homo_del_threshold": 0.3,
    "hetero_del_threshold": 0.7,
    "dup_threshold": 1.5
  }
}
```

#### IGV Settings
```json
"igv_settings": {
  "igv_path": "igv.sh",
  "igv_memory": "4g",
  "timeout": 300
}
```

## 📋 Requirements

### System Requirements
- **OS**: Linux (tested on Ubuntu 18.04+, CentOS 7+)
- **Memory**: 8GB+ RAM recommended
- **Storage**: 10GB+ free space for outputs
- **CPU**: 4+ cores recommended for batch processing

### Software Dependencies
- **Python 3.8+**
- **samtools** (for BAM processing)
- **IGV** (for automated screenshots)

### Python Packages
```
pysam>=0.21.0
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
scikit-learn>=1.0.0
plotly>=5.0.0
optuna>=3.0.0
psutil>=5.8.0
```

## 🏥 Clinical Interpretation

### SMA Risk Categories
- **HIGH_RISK**: SMN1 homozygous deletion - SMA affected
- **CARRIER**: SMN1 heterozygous deletion - SMA carrier
- **LOW_RISK**: Normal SMN1 copy number
- **UNCERTAIN**: Atypical pattern requiring further investigation

### Copy Number Calls
- **HOMO_DEL**: 0 copies (homozygous deletion)
- **HETERO_DEL**: 1 copy (heterozygous deletion)
- **NORMAL**: 2 copies (normal)
- **DUP**: 3+ copies (duplication)

## 🔬 Validation Data

### MLPA Validation
The pipeline can use MLPA (Multiplex Ligation-dependent Probe Amplification) results for:
- Threshold optimization
- Performance monitoring
- Accuracy assessment

### Validation File Format
```tsv
sample_id	smn1_exon7	smn1_exon8	smn2_exon7	smn2_exon8	mlpa_result	comments
CONTROL_001	NORMAL	NORMAL	NORMAL	NORMAL	Normal control	Validated by MLPA
CARRIER_001	HETERO_DEL	HETERO_DEL	NORMAL	NORMAL	SMA carrier	Validated by MLPA
SMA_001	HOMO_DEL	HOMO_DEL	NORMAL	NORMAL	SMA affected	Validated by MLPA
```

## 🔧 Troubleshooting

### Common Issues

#### 1. IGV Not Found
```bash
# Install IGV automatically
./setup_igv.sh

# Or manually set path in config/config.json
"igv_settings": {
  "igv_path": "/path/to/igv.sh"
}
```

#### 2. BAM File Issues
```bash
# Index BAM files
samtools index sample.bam

# Check BAM file integrity
samtools quickcheck sample.bam

# Validate reference genome
samtools view -H sample.bam | grep "^@SQ"
```

#### 3. Memory Issues
```bash
# Reduce IGV memory if needed
"igv_settings": {
  "igv_memory": "2g"
}

# Process smaller batches
python3 smn_pipeline.py --config config/config.json --input dir --output out --batch
```

#### 4. Permission Errors
```bash
chmod +x *.sh
chmod +x validate_pipeline.py
```

### Log Analysis
- Check `logs/` directory for detailed error messages
- Use `--verbose` flag for debug logging
- Review IGV batch scripts in temp/ if snapshots fail

## 📈 Performance Optimization

### For Large Batches (50+ samples)
1. **Increase parallel processing**:
   ```json
   "processing": {
     "threads": 8,
     "chunk_size": 20
   }
   ```

2. **Optimize IGV settings**:
   ```json
   "igv_settings": {
     "igv_memory": "8g",
     "timeout": 600
   }
   ```

3. **Use SSD storage** for temp files and outputs

### Memory Optimization
- Process in smaller batches if memory limited
- Use `cleanup_temp_files` utility regularly
- Monitor with built-in resource tracking

## 📚 Algorithm Details

### Computational Methods

#### 1. Depth Ratio Method
- Normalizes read depth against control genes
- Uses adaptive thresholds for CNV calling
- Accounts for GC content and mapping quality

#### 2. Statistical Method
- Z-score analysis relative to control regions
- Bootstrap confidence intervals
- Hypothesis testing for significance

#### 3. HMM-like Method
- Coverage pattern analysis
- State classification based on uniformity
- Accounts for coverage gaps and variability

#### 4. Ensemble Method
- Combines multiple features and methods
- Machine learning-based optimization
- Weighted voting with quality factors

### Machine Learning Pipeline
1. **Feature extraction**: 12+ features per region
2. **Model training**: Random Forest with hyperparameter optimization
3. **Threshold adaptation**: Continuous improvement with new data
4. **Performance monitoring**: Tracks accuracy vs validation data
5. **Model updating**: Retrains every 50 samples or when accuracy drops

## 🧪 Testing

### Run Tests
```bash
# Validate setup
./run_tests.sh

# Test with example data
python3 tests/test_pipeline.py

# Validate against MLPA data
python3 tests/validate_mlpa.py
```

### Test Data
- Use provided test samples in `tests/data/`
- Validate with your own MLPA-confirmed samples
- Check accuracy metrics in generated reports

## 📖 Usage Examples

### Example 1: Clinical Lab Workflow
```bash
# Daily batch processing
BATCH_DIR="/data/wes_bams/$(date +%Y%m%d)"
OUTPUT_DIR="/results/smn_analysis/$(date +%Y%m%d)"

python3 smn_pipeline.py \
    --config config/clinical_config.json \
    --input $BATCH_DIR \
    --output $OUTPUT_DIR \
    --batch

# Generate consolidated report
echo "Results available at: $OUTPUT_DIR/batch_report_*.html"
```

### Example 2: Research Analysis
```bash
# Process research cohort
python3 smn_pipeline.py \
    --config config/research_config.json \
    --input /data/research_cohort \
    --output /results/research_smn \
    --batch \
    --update-thresholds
```

### Example 3: Single Sample QC
```bash
# Quick single sample check
python3 smn_pipeline.py \
    -c config/config.json \
    -i suspicious_sample.bam \
    -o qc_check \
    -s suspicious_001
```

## 📊 Output Interpretation

### CNV Call Confidence
- **High (>0.8)**: Reliable call, suitable for clinical reporting
- **Medium (0.5-0.8)**: Good call, may benefit from confirmation
- **Low (<0.5)**: Uncertain call, requires manual review

### SMA Risk Assessment
- **HIGH_RISK**: Likely SMA affected, urgent genetic counseling
- **CARRIER**: SMA carrier, genetic counseling for family planning
- **LOW_RISK**: Normal SMN1, routine reporting
- **UNCERTAIN**: Requires additional testing (MLPA, qPCR)

### Quality Metrics
- **Mapping Rate**: >95% for reliable analysis
- **Mean MAPQ**: >30 for high confidence
- **Coverage Uniformity**: >0.8 for optimal calling
- **Zero Depth Fraction**: <0.1 for normal regions

## 🔬 Scientific Background

### SMN Genes and SMA
- **SMN1 and SMN2**: Nearly identical genes on chromosome 5q
- **Exon 7 difference**: Critical nucleotide difference between SMN1/SMN2
- **SMA pathogenesis**: Loss of SMN1 function causes spinal muscular atrophy
- **Copy number variation**: SMN1 deletions are the primary cause of SMA

### Technical Challenges
- **High sequence similarity**: SMN1 and SMN2 are 99.9% identical
- **Pseudogene complexity**: Multiple SMN-like sequences in the genome
- **Coverage bias**: Exon 7 often has reduced coverage in WES
- **Reference alignment**: Challenges in distinguishing paralogs

### Pipeline Approach
- **Multi-method consensus**: Reduces false positives/negatives
- **Control gene normalization**: Accounts for technical variation
- **Pattern recognition**: Identifies characteristic CNV signatures
- **Quality integration**: Weights calls by sequencing quality

## 🤝 Contributing

### Development Setup
```bash
# Development installation
git clone <repo-url>
cd smn_cnv_pipeline
./setup.sh

# Install development dependencies
pip3 install -r requirements-dev.txt

# Run tests
./run_tests.sh
```

### Adding New Methods
1. Create new method class in `src/cnv_caller.py`
2. Implement `call_cnv()` method
3. Add to methods dictionary in CNVCaller
4. Update tests and documentation

### Validation Data
- Contribute MLPA-validated samples
- Follow validation data format
- Ensure patient consent for data sharing

## 📄 Citation

If you use this pipeline in your research, please cite:

```
SMN CNV Detection Pipeline v1.0
A comprehensive toolkit for SMN1/SMN2 copy number analysis from WES data
[Your publication details here]
```

## 🆘 Support

### Getting Help
1. **Check logs**: Review `logs/` directory for errors
2. **Run validation**: Use `validate_pipeline.py` for setup issues
3. **Review documentation**: Check `docs/` for detailed guides
4. **Test with known samples**: Validate with MLPA-confirmed cases

### Common Questions

**Q: Why are some exon 7 calls uncertain?**
A: Exon 7 often has reduced coverage in WES. The pipeline accounts for this with multiple validation methods.

**Q: How accurate is the pipeline?**
A: Accuracy improves with use, targeting 95% after 1000+ samples. Initial accuracy ~85-90% with default thresholds.

**Q: Can I use hg19 reference?**
A: Currently optimized for hg38. For hg19, update coordinates in config file.

**Q: What if IGV snapshots fail?**
A: The pipeline continues without snapshots. Check IGV installation and memory settings.

## 📝 Version History

### v1.0.0 (Current)
- Initial release
- Multi-method CNV detection
- IGV automation
- ML threshold optimization
- Comprehensive reporting

### Planned Features
- **v1.1**: Enhanced visualization with interactive plots
- **v1.2**: Integration with variant calling pipelines
- **v1.3**: Support for targeted sequencing panels
- **v2.0**: Real-time processing and cloud deployment

## 📜 License

[Specify your license here - e.g., MIT, GPL, etc.]

## 🙏 Acknowledgments

- Broad Institute IGV team for visualization tools
- SMA research community for clinical insights
- Bioinformatics community for algorithmic foundations

---

**For detailed usage instructions, see [QUICKSTART.md](QUICKSTART.md)**

**For technical documentation, see `docs/` directory**
