"""
Test suite for SMN CNV Pipeline
===============================

This package contains unit tests, integration tests, and validation tests
for the SMN CNV detection pipeline.

Test Categories:
---------------
- Unit tests: Test individual components and functions
- Integration tests: Test component interactions and workflows
- Validation tests: Test against known samples and MLPA data
- Performance tests: Test pipeline performance and resource usage

Usage:
------
Run all tests:
    python -m pytest tests/

Run specific test category:
    python -m pytest tests/unit/
    python -m pytest tests/integration/
    python -m pytest tests/validation/

Run with coverage:
    python -m pytest tests/ --cov=src --cov-report=html
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

# Add src to path for testing
TEST_DIR = Path(__file__).parent
PROJECT_ROOT = TEST_DIR.parent
SRC_DIR = PROJECT_ROOT / 'src'

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# Test configuration
TEST_CONFIG = {
    'test_data_dir': TEST_DIR / 'data',
    'temp_dir': None,  # Will be set per test
    'mock_config': {
        'pipeline_version': '1.0.0-test',
        'output_dir': 'test_output',
        'model_dir': 'test_models',
        'smn_coordinates': {
            'smn1_exon7': {'chr': 'chr5', 'start': 70247724, 'end': 70247775},
            'smn1_exon8': {'chr': 'chr5', 'start': 70248935, 'end': 70249306},
            'smn2_exon7': {'chr': 'chr5', 'start': 69372304, 'end': 69372355},
            'smn2_exon8': {'chr': 'chr5', 'start': 69373515, 'end': 69373886}
        }
    }
}

class TestSetupMixin:
    """Mixin class for common test setup and teardown"""
    
    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp(prefix='smn_test_')
        self.test_config = TEST_CONFIG['mock_config'].copy()
        self.test_config['output_dir'] = os.path.join(self.temp_dir, 'output')
        self.test_config['model_dir'] = os.path.join(self.temp_dir, 'models')
        
        # Create test directories
        os.makedirs(self.test_config['output_dir'], exist_ok=True)
        os.makedirs(self.test_config['model_dir'], exist_ok=True)
    
    def tearDown(self):
        """Clean up test environment"""
        if hasattr(self, 'temp_dir') and os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir, ignore_errors=True)

def create_mock_depth_metrics():
    """Create mock depth metrics for testing"""
    return {
        'smn1_exon7': {
            'mean_depth': 15.0, 'median_depth': 14.0, 'std_depth': 3.0,
            'coverage_uniformity': 0.8, 'zero_depth_positions': 2, 'low_depth_positions': 5,
            'region_length': 52, 'mean_mapq': 35.0, 'mean_baseq': 30.0, 'gc_content': 0.6,
            'min_depth': 5, 'max_depth': 25
        },
        'smn1_exon8': {
            'mean_depth': 30.0, 'median_depth': 29.0, 'std_depth': 5.0,
            'coverage_uniformity': 0.9, 'zero_depth_positions': 0, 'low_depth_positions': 2,
            'region_length': 371, 'mean_mapq': 40.0, 'mean_baseq': 32.0, 'gc_content': 0.5,
            'min_depth': 15, 'max_depth': 45
        },
        'smn2_exon7': {
            'mean_depth': 25.0, 'median_depth': 24.0, 'std_depth': 4.0,
            'coverage_uniformity': 0.85, 'zero_depth_positions': 1, 'low_depth_positions': 3,
            'region_length': 52, 'mean_mapq': 38.0, 'mean_baseq': 31.0, 'gc_content': 0.55,
            'min_depth': 10, 'max_depth': 35
        },
        'smn2_exon8': {
            'mean_depth': 28.0, 'median_depth': 27.0, 'std_depth': 6.0,
            'coverage_uniformity': 0.82, 'zero_depth_positions': 0, 'low_depth_positions': 1,
            'region_length': 371, 'mean_mapq': 36.0, 'mean_baseq': 30.0, 'gc_content': 0.52,
            'min_depth': 12, 'max_depth': 42
        },
        'normalized_depths': {
            'smn1_exon7': {'normalized_depth': 0.5, 'zscore': -2.5, 'depth_ratio': 0.5},
            'smn1_exon8': {'normalized_depth': 1.0, 'zscore': 0.0, 'depth_ratio': 1.0},
            'smn2_exon7': {'normalized_depth': 0.83, 'zscore': -1.0, 'depth_ratio': 0.83},
            'smn2_exon8': {'normalized_depth': 0.93, 'zscore': -0.3, 'depth_ratio': 0.93}
        },
        'overall': {
            'total_reads_sampled': 10000,
            'mapping_rate': 0.96,
            'mean_mapq': 37.0,
            'gc_content': 0.54
        },
        'controls': {
            'actb_exon1': {'mean_depth': 30.0},
            'gapdh_exon1': {'mean_depth': 32.0},
            'tbp_exon1': {'mean_depth': 28.0}
        }
    }

def create_mock_cnv_results():
    """Create mock CNV results for testing"""
    return {
        'smn1_exon7': {
            'call': 'HETERO_DEL',
            'confidence': 0.85,
            'methods': {
                'depth_ratio': {'call': 'HETERO_DEL', 'confidence': 0.87},
                'statistical': {'call': 'HETERO_DEL', 'confidence': 0.82},
                'hmm_like': {'call': 'HETERO_DEL', 'confidence': 0.79},
                'ensemble': {'call': 'HETERO_DEL', 'confidence': 0.86}
            }
        },
        'smn1_exon8': {
            'call': 'NORMAL',
            'confidence': 0.92,
            'methods': {
                'depth_ratio': {'call': 'NORMAL', 'confidence': 0.94},
                'statistical': {'call': 'NORMAL', 'confidence': 0.91},
                'hmm_like': {'call': 'NORMAL', 'confidence': 0.89},
                'ensemble': {'call': 'NORMAL', 'confidence': 0.93}
            }
        },
        'smn2_exon7': {
            'call': 'NORMAL',
            'confidence': 0.88,
            'methods': {
                'depth_ratio': {'call': 'NORMAL', 'confidence': 0.89},
                'statistical': {'call': 'NORMAL', 'confidence': 0.87},
                'hmm_like': {'call': 'NORMAL', 'confidence': 0.86},
                'ensemble': {'call': 'NORMAL', 'confidence': 0.90}
            }
        },
        'smn2_exon8': {
            'call': 'NORMAL',
            'confidence': 0.90,
            'methods': {
                'depth_ratio': {'call': 'NORMAL', 'confidence': 0.91},
                'statistical': {'call': 'NORMAL', 'confidence': 0.89},
                'hmm_like': {'call': 'NORMAL', 'confidence': 0.88},
                'ensemble': {'call': 'NORMAL', 'confidence': 0.92}
            }
        },
        'sma_risk': 'CARRIER'
    }

def create_mock_sample_result():
    """Create complete mock sample result for testing"""
    return {
        'sample_id': 'TEST_SAMPLE_001',
        'bam_file': 'test_sample.bam',
        'depth_metrics': create_mock_depth_metrics(),
        'cnv_calls': create_mock_cnv_results(),
        'igv_snapshots': {
            'smn1_exon7': 'snapshots/TEST_SAMPLE_001_smn1_exon7.png',
            'smn1_exon8': 'snapshots/TEST_SAMPLE_001_smn1_exon8.png',
            'smn2_exon7': 'snapshots/TEST_SAMPLE_001_smn2_exon7.png',
            'smn2_exon8': 'snapshots/TEST_SAMPLE_001_smn2_exon8.png'
        },
        'processing_time': '2024-01-01T12:00:00',
        'pipeline_version': '1.0.0-test'
    }

def create_mock_validation_data():
    """Create mock MLPA validation data"""
    return [
        {
            'sample_id': 'VAL_NORMAL_001',
            'mlpa_result': {
                'smn1_exon7': 'NORMAL',
                'smn1_exon8': 'NORMAL',
                'smn2_exon7': 'NORMAL',
                'smn2_exon8': 'NORMAL'
            }
        },
        {
            'sample_id': 'VAL_CARRIER_001',
            'mlpa_result': {
                'smn1_exon7': 'HETERO_DEL',
                'smn1_exon8': 'HETERO_DEL',
                'smn2_exon7': 'NORMAL',
                'smn2_exon8': 'NORMAL'
            }
        },
        {
            'sample_id': 'VAL_SMA_001',
            'mlpa_result': {
                'smn1_exon7': 'HOMO_DEL',
                'smn1_exon8': 'HOMO_DEL',
                'smn2_exon7': 'NORMAL',
                'smn2_exon8': 'NORMAL'
            }
        }
    ]

# Test utilities
def skip_if_no_dependencies():
    """Decorator to skip tests if dependencies are not available"""
    def decorator(test_func):
        import functools
        @functools.wraps(test_func)
        def wrapper(*args, **kwargs):
            try:
                import pysam, numpy, pandas, scipy, sklearn
                return test_func(*args, **kwargs)
            except ImportError as e:
                import unittest
                raise unittest.SkipTest(f"Dependencies not available: {e}")
        return wrapper
    return decorator

def requires_test_data():
    """Decorator to skip tests if test data is not available"""
    def decorator(test_func):
        import functools
        @functools.wraps(test_func)
        def wrapper(*args, **kwargs):
            test_data_dir = TEST_CONFIG['test_data_dir']
            if not test_data_dir.exists():
                import unittest
                raise unittest.SkipTest("Test data directory not found")
            return test_func(*args, **kwargs)
        return wrapper
    return decorator

__all__ = [
    'TEST_CONFIG',
    'TestSetupMixin',
    'create_mock_depth_metrics',
    'create_mock_cnv_results',
    'create_mock_sample_result',
    'create_mock_validation_data',
    'skip_if_no_dependencies',
    'requires_test_data'
]
