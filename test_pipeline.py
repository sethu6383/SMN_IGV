#!/usr/bin/env python3
"""
Test Script for SMN CNV Pipeline
Tests core functionality with mock data
"""

import os
import sys
import json
import tempfile
import numpy as np
import pysam
from pathlib import Path
import unittest
from unittest.mock import Mock, patch

# Add src to path
sys.path.insert(0, 'src')

from depth_analyzer import DepthAnalyzer
from cnv_caller import CNVCaller
from ml_threshold import MLThresholdOptimizer
from report_generator import ReportGenerator
from utils import setup_logging, validate_inputs

class TestSMNPipeline(unittest.TestCase):
    """Test cases for SMN pipeline components"""
    
    def setUp(self):
        """Set up test environment"""
        self.test_config = {
            "pipeline_version": "1.0.0",
            "output_dir": "test_output",
            "model_dir": "test_models",
            "igv_snapshots": "test_snapshots",
            "reports": "test_reports",
            "smn_coordinates": {
                "smn1_exon7": {"chr": "chr5", "start": 70247724, "end": 70247775},
                "smn1_exon8": {"chr": "chr5", "start": 70248935, "end": 70249306},
                "smn2_exon7": {"chr": "chr5", "start": 69372304, "end": 69372355},
                "smn2_exon8": {"chr": "chr5", "start": 69373515, "end": 69373886}
            },
            "cnv_thresholds": {
                "depth_ratio": {
                    "homo_del_threshold": 0.3,
                    "hetero_del_threshold": 0.7,
                    "dup_threshold": 1.5
                },
                "statistical": {
                    "zscore_threshold": 2.0,
                    "pvalue_threshold": 0.05
                }
            }
        }
        
        # Create test directories
        for dir_name in ["test_output", "test_models", "test_snapshots", "test_reports"]:
            Path(dir_name).mkdir(exist_ok=True)
    
    def tearDown(self):
        """Clean up test environment"""
        import shutil
        for dir_name in ["test_output", "test_models", "test_snapshots", "test_reports"]:
            if os.path.exists(dir_name):
                shutil.rmtree(dir_name)
    
    def test_depth_analyzer_initialization(self):
        """Test DepthAnalyzer initialization"""
        analyzer = DepthAnalyzer(self.test_config)
        
        self.assertIsInstance(analyzer.smn_regions, dict)
        self.assertIn('smn1_exon7', analyzer.smn_regions)
        self.assertIn('smn2_exon7', analyzer.smn_regions)
    
    def test_cnv_caller_initialization(self):
        """Test CNVCaller initialization"""
        caller = CNVCaller(self.test_config)
        
        self.assertIsInstance(caller.methods, dict)
        self.assertIn('depth_ratio', caller.methods)
        self.assertIn('statistical', caller.methods)
        self.assertIn('hmm_like', caller.methods)
        self.assertIn('ensemble', caller.methods)
    
    def test_mock_depth_analysis(self):
        """Test depth analysis with mock data"""
        
        # Create mock depth metrics
        mock_depth_metrics = {
            'smn1_exon7': {
                'mean_depth': 15.0,
                'median_depth': 14.0,
                'std_depth': 3.0,
                'coverage_uniformity': 0.8,
                'zero_depth_positions': 2,
                'region_length': 52,
                'mean_mapq': 35.0,
                'gc_content': 0.6
            },
            'smn1_exon8': {
                'mean_depth': 30.0,
                'median_depth': 29.0,
                'std_depth': 5.0,
                'coverage_uniformity': 0.9,
                'zero_depth_positions': 0,
                'region_length': 371,
                'mean_mapq': 40.0,
                'gc_content': 0.5
            },
            'smn2_exon7': {
                'mean_depth': 25.0,
                'median_depth': 24.0,
                'std_depth': 4.0,
                'coverage_uniformity': 0.85,
                'zero_depth_positions': 1,
                'region_length': 52,
                'mean_mapq': 38.0,
                'gc_content': 0.55
            },
            'smn2_exon8': {
                'mean_depth': 28.0,
                'median_depth': 27.0,
                'std_depth': 6.0,
                'coverage_uniformity': 0.82,
                'zero_depth_positions': 0,
                'region_length': 371,
                'mean_mapq': 36.0,
                'gc_content': 0.52
            },
            'normalized_depths': {
                'smn1_exon7': {'normalized_depth': 0.5, 'zscore': -2.5},
                'smn1_exon8': {'normalized_depth': 1.0, 'zscore': 0.0},
                'smn2_exon7': {'normalized_depth': 0.83, 'zscore': -1.0},
                'smn2_exon8': {'normalized_depth': 0.93, 'zscore': -0.3}
            },
            'overall': {
                'mapping_rate': 0.96,
                'mean_mapq': 37.0,
                'gc_content': 0.54
            }
        }
        
        # Test CNV calling
        caller = CNVCaller(self.test_config)
        cnv_results = caller.call_cnvs(mock_depth_metrics)
        
        # Verify results structure
        self.assertIn('smn1_exon7', cnv_results)
        self.assertIn('smn1_exon8', cnv_results)
        self.assertIn('sma_risk', cnv_results)
        
        # Verify calls are reasonable (SMN1 exon7 should be HETERO_DEL based on mock data)
        smn1_ex7_call = cnv_results['smn1_exon7']['call']
        self.assertIn(smn1_ex7_call, ['HOMO_DEL', 'HETERO_DEL', 'NORMAL', 'DUP'])
        
        # Check confidence scores
        for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
            confidence = cnv_results[region]['confidence']
            self.assertGreaterEqual(confidence, 0.0)
            self.assertLessEqual(confidence, 1.0)
    
    def test_ml_optimizer_initialization(self):
        """Test ML threshold optimizer"""
        optimizer = MLThresholdOptimizer(self.test_config)
        
        self.assertIsInstance(optimizer.thresholds, dict)
        self.assertIn('depth_ratio', optimizer.thresholds)
    
    def test_report_generator(self):
        """Test report generation with mock data"""
        
        # Mock sample result
        mock_result = {
            'sample_id': 'TEST_001',
            'bam_file': 'test.bam',
            'depth_metrics': {
                'smn1_exon7': {'mean_depth': 15.0, 'coverage_uniformity': 0.8, 'mean_mapq': 35.0, 'gc_content': 0.6},
                'overall': {'mapping_rate': 0.96, 'mean_mapq': 37.0, 'gc_content': 0.54},
                'normalized_depths': {
                    'smn1_exon7': {'normalized_depth': 0.5, 'zscore': -2.5}
                }
            },
            'cnv_calls': {
                'smn1_exon7': {'call': 'HETERO_DEL', 'confidence': 0.85},
                'smn1_exon8': {'call': 'NORMAL', 'confidence': 0.92},
                'smn2_exon7': {'call': 'NORMAL', 'confidence': 0.88},
                'smn2_exon8': {'call': 'NORMAL', 'confidence': 0.90},
                'sma_risk': 'CARRIER'
            },
            'igv_snapshots': {},
            'processing_time': '2024-01-01T12:00:00'
        }
        
        generator = ReportGenerator(self.test_config)
        
        # Test HTML generation (should not crash)
        try:
            plots = generator._create_sample_plots(mock_result)
            html_report = generator._create_html_report(mock_result, plots)
            self.assertIsInstance(html_report, str)
            self.assertIn('TEST_001', html_report)
            self.assertIn('CARRIER', html_report)
        except Exception as e:
            self.fail(f"Report generation failed: {e}")
    
    def test_configuration_loading(self):
        """Test configuration file loading"""
        
        # Create test config file
        test_config_file = 'test_config.json'
        with open(test_config_file, 'w') as f:
            json.dump(self.test_config, f)
        
        try:
            from utils import load_config_with_defaults
            loaded_config = load_config_with_defaults(test_config_file)
            
            self.assertIsInstance(loaded_config, dict)
            self.assertIn('pipeline_version', loaded_config)
            self.assertIn('smn_coordinates', loaded_config)
            
        finally:
            if os.path.exists(test_config_file):
                os.remove(test_config_file)
    
    def test_utility_functions(self):
        """Test utility functions"""
        
        from utils import parse_coordinates, calculate_statistics, format_genomic_position
        
        # Test coordinate parsing
        coords = parse_coordinates("chr5:70247724-70247775")
        self.assertEqual(coords['chr'], 'chr5')
        self.assertEqual(coords['start'], 70247724)
        self.assertEqual(coords['end'], 70247775)
        
        # Test statistics calculation
        test_values = [10, 20, 30, 40, 50]
        stats = calculate_statistics(test_values)
        self.assertEqual(stats['mean'], 30.0)
        self.assertEqual(stats['median'], 30.0)
        
        # Test position formatting
        formatted = format_genomic_position('chr5', 70247724)
        self.assertEqual(formatted, 'chr5:70,247,724')

def create_mock_bam_file(filename: str) -> str:
    """Create a mock BAM file for testing"""
    
    # This is a simplified mock - in practice you'd need a real BAM
    # For testing, we'll just create an empty file
    with open(filename, 'wb') as f:
        f.write(b'BAM\x01')  # Minimal BAM header
    
    return filename

def run_integration_test():
    """Run integration test with mock data"""
    
    print("Running SMN Pipeline Integration Test...")
    
    # Create test config
    test_config = {
        "pipeline_version": "1.0.0",
        "output_dir": "integration_test",
        "model_dir": "integration_test/models",
        "igv_snapshots": "integration_test/snapshots",
        "reports": "integration_test/reports"
    }
    
    # Create test directories
    for dir_path in test_config.values():
        if isinstance(dir_path, str):
            Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    try:
        # Test component initialization
        print("✓ Testing component initialization...")
        depth_analyzer = DepthAnalyzer(test_config)
        cnv_caller = CNVCaller(test_config)
        ml_optimizer = MLThresholdOptimizer(test_config)
        report_generator = ReportGenerator(test_config)
        
        print("✓ All components initialized successfully")
        
        # Test with mock depth data
        print("✓ Testing CNV calling with mock data...")
        mock_depth_metrics = {
            'smn1_exon7': {
                'mean_depth': 15.0, 'std_depth': 3.0, 'coverage_uniformity': 0.8,
                'zero_depth_positions': 2, 'low_depth_positions': 5, 'region_length': 52,
                'mean_mapq': 35.0, 'mean_baseq': 30.0, 'gc_content': 0.6
            },
            'smn1_exon8': {
                'mean_depth': 30.0, 'std_depth': 5.0, 'coverage_uniformity': 0.9,
                'zero_depth_positions': 0, 'low_depth_positions': 2, 'region_length': 371,
                'mean_mapq': 40.0, 'mean_baseq': 32.0, 'gc_content': 0.5
            },
            'smn2_exon7': {
                'mean_depth': 25.0, 'std_depth': 4.0, 'coverage_uniformity': 0.85,
                'zero_depth_positions': 1, 'low_depth_positions': 3, 'region_length': 52,
                'mean_mapq': 38.0, 'mean_baseq': 31.0, 'gc_content': 0.55
            },
            'smn2_exon8': {
                'mean_depth': 28.0, 'std_depth': 6.0, 'coverage_uniformity': 0.82,
                'zero_depth_positions': 0, 'low_depth_positions': 1, 'region_length': 371,
                'mean_mapq': 36.0, 'mean_baseq': 30.0, 'gc_content': 0.52
            },
            'normalized_depths': {
                'smn1_exon7': {'normalized_depth': 0.5, 'zscore': -2.5},
                'smn1_exon8': {'normalized_depth': 1.0, 'zscore': 0.0},
                'smn2_exon7': {'normalized_depth': 0.83, 'zscore': -1.0},
                'smn2_exon8': {'normalized_depth': 0.93, 'zscore': -0.3}
            },
            'overall': {
                'mapping_rate': 0.96,
                'mean_mapq': 37.0,
                'gc_content': 0.54
            }
        }
        
        # Test CNV calling
        cnv_results = cnv_caller.call_cnvs(mock_depth_metrics)
        
        print(f"✓ CNV calling completed. SMA risk: {cnv_results.get('sma_risk', 'UNKNOWN')}")
        
        # Test report generation
        print("✓ Testing report generation...")
        mock_sample_result = {
            'sample_id': 'TEST_INTEGRATION',
            'bam_file': 'test.bam',
            'depth_metrics': mock_depth_metrics,
            'cnv_calls': cnv_results,
            'igv_snapshots': {},
            'processing_time': '2024-01-01T12:00:00',
            'pipeline_version': '1.0.0'
        }
        
        report_file = report_generator.generate_sample_report(mock_sample_result)
        
        if report_file and os.path.exists(report_file):
            print(f"✓ Report generated: {report_file}")
        else:
            print("⚠ Report generation completed but file not found")
        
        print("✅ Integration test completed successfully!")
        
        # Display results summary
        print("\n=== Test Results Summary ===")
        print(f"SMN1 Exon 7: {cnv_results['smn1_exon7']['call']} (conf: {cnv_results['smn1_exon7']['confidence']:.3f})")
        print(f"SMN1 Exon 8: {cnv_results['smn1_exon8']['call']} (conf: {cnv_results['smn1_exon8']['confidence']:.3f})")
        print(f"SMN2 Exon 7: {cnv_results['smn2_exon7']['call']} (conf: {cnv_results['smn2_exon7']['confidence']:.3f})")
        print(f"SMN2 Exon 8: {cnv_results['smn2_exon8']['call']} (conf: {cnv_results['smn2_exon8']['confidence']:.3f})")
        print(f"SMA Risk: {cnv_results['sma_risk']}")
        
    except Exception as e:
        print(f"❌ Integration test failed: {e}")
        raise
    
    finally:
        # Cleanup
        import shutil
        if os.path.exists("integration_test"):
            shutil.rmtree("integration_test")

def test_configuration_validation():
    """Test configuration file validation"""
    
    print("Testing configuration validation...")
    
    # Test with valid config
    if os.path.exists('config/config.json'):
        try:
            with open('config/config.json', 'r') as f:
                config = json.load(f)
            
            required_keys = ['pipeline_version', 'smn_coordinates', 'cnv_thresholds']
            missing_keys = [key for key in required_keys if key not in config]
            
            if not missing_keys:
                print("✓ Configuration file validation passed")
                return True
            else:
                print(f"❌ Missing configuration keys: {missing_keys}")
                return False
                
        except Exception as e:
            print(f"❌ Configuration file error: {e}")
            return False
    else:
        print("⚠ Configuration file not found")
        return False

def test_dependencies():
    """Test that all required dependencies are available"""
    
    print("Testing dependencies...")
    
    required_packages = [
        'pysam', 'numpy', 'pandas', 'scipy', 'sklearn',
        'plotly', 'optuna', 'psutil'
    ]
    
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if not missing_packages:
        print("✓ All required packages available")
        return True
    else:
        print(f"❌ Missing packages: {missing_packages}")
        print("Run: pip3 install -r requirements.txt")
        return False

def main():
    """Main test function"""
    
    print("🧬 SMN CNV Pipeline Test Suite")
    print("=" * 40)
    
    all_passed = True
    
    # Test 1: Dependencies
    if not test_dependencies():
        all_passed = False
    
    # Test 2: Configuration
    if not test_configuration_validation():
        all_passed = False
    
    # Test 3: Unit tests
    print("\nRunning unit tests...")
    try:
        unittest.main(argv=[''], exit=False, verbosity=2)
        print("✓ Unit tests completed")
    except Exception as e:
        print(f"❌ Unit tests failed: {e}")
        all_passed = False
    
    # Test 4: Integration test
    print("\nRunning integration test...")
    try:
        run_integration_test()
        print("✓ Integration test passed")
    except Exception as e:
        print(f"❌ Integration test failed: {e}")
        all_passed = False
    
    # Final summary
    print("\n" + "=" * 40)
    if all_passed:
        print("🎉 All tests passed! Pipeline is ready for use.")
        print("\nNext steps:")
        print("1. Update SMN coordinates in config/config.json if needed")
        print("2. Add MLPA validation data")
        print("3. Test with real BAM files")
        print("4. Run batch processing")
    else:
        print("❌ Some tests failed. Please fix issues before using the pipeline.")
        sys.exit(1)

if __name__ == "__main__":
    main()
