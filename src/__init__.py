"""
SMN CNV Pipeline Package
========================

A comprehensive bioinformatics pipeline for detecting copy number variations (CNVs) 
in SMN1 and SMN2 genes from whole exome sequencing (WES) data.

Modules:
--------
- depth_analyzer: Analyzes read depth and quality metrics for SMN regions
- cnv_caller: Multi-method CNV calling with ensemble approaches
- igv_automation: Automated IGV screenshot generation
- ml_threshold: Machine learning-based threshold optimization
- report_generator: HTML and TSV report generation
- utils: Utility functions and helpers

Version: 1.0.0
Author: SMN Pipeline Development Team
"""

__version__ = "1.0.0"
__author__ = "SMN Pipeline Development Team"

# Import main classes for easy access
try:
    from .depth_analyzer import DepthAnalyzer
    from .cnv_caller import CNVCaller
    from .igv_automation import IGVAutomation
    from .ml_threshold import MLThresholdOptimizer
    from .report_generator import ReportGenerator
    from .utils import setup_logging, validate_inputs
    
    __all__ = [
        'DepthAnalyzer',
        'CNVCaller', 
        'IGVAutomation',
        'MLThresholdOptimizer',
        'ReportGenerator',
        'setup_logging',
        'validate_inputs'
    ]
    
except ImportError as e:
    # Graceful handling of import errors during development/testing
    import warnings
    warnings.warn(f"Some modules could not be imported: {e}", ImportWarning)
    __all__ = []

# Module metadata
PIPELINE_INFO = {
    'name': 'SMN CNV Pipeline',
    'version': __version__,
    'description': 'CNV detection pipeline for SMN1/SMN2 genes',
    'author': __author__,
    'license': 'MIT',
    'python_requires': '>=3.8',
    'dependencies': [
    'dependencies': [
        'pysam>=0.21.0',
        'numpy>=1.21.0',
        'pandas>=1.5.0',
        'scipy>=1.9.0',
        'scikit-learn>=1.1.0',
        'plotly>=5.10.0',
        'optuna>=3.0.0',
        'psutil>=5.8.0'
    ],
    'supported_platforms': ['Linux'],
    'supported_references': ['hg38'],
    'input_formats': ['BAM'],
    'output_formats': ['HTML', 'TSV', 'PNG']
}

def get_pipeline_info():
    """Get pipeline information dictionary"""
    return PIPELINE_INFO.copy()

def check_dependencies():
    """Check if all required dependencies are available"""
    missing_deps = []
    required_packages = [
        'pysam', 'numpy', 'pandas', 'scipy', 'sklearn',
        'plotly', 'optuna', 'psutil'
    ]
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_deps.append(package)
    
    return {
        'all_available': len(missing_deps) == 0,
        'missing': missing_deps,
        'available': [p for p in required_packages if p not in missing_deps]
    }
