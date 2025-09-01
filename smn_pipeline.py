#!/usr/bin/env python3
"""
SMN1/SMN2 CNV Detection Pipeline
Main pipeline script for detecting copy number variations in SMN genes
"""

import os
import sys
import json
import argparse
import logging
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import subprocess
import pysam
from datetime import datetime
import pickle

# Import custom modules
from src.depth_analyzer import DepthAnalyzer
from src.igv_automation import IGVAutomation
from src.ml_threshold import MLThresholdOptimizer
from src.report_generator import ReportGenerator
from src.cnv_caller import CNVCaller
from src.utils import setup_logging, validate_inputs

class SMNPipeline:
    """Main SMN CNV detection pipeline class"""
    
    def __init__(self, config_file: str):
        """Initialize pipeline with configuration"""
        with open(config_file, 'r') as f:
            self.config = json.load(f)
        
        self.setup_directories()
        self.logger = setup_logging(self.config['output_dir'])
        
        # Initialize components
        self.depth_analyzer = DepthAnalyzer(self.config)
        self.igv_automation = IGVAutomation(self.config)
        self.ml_optimizer = MLThresholdOptimizer(self.config)
        self.report_generator = ReportGenerator(self.config)
        self.cnv_caller = CNVCaller(self.config)
        
        self.results = []
    
    def setup_directories(self):
        """Create necessary output directories"""
        dirs = ['output_dir', 'igv_snapshots', 'reports', 'temp_files']
        for dir_key in dirs:
            if dir_key in self.config:
                Path(self.config[dir_key]).mkdir(parents=True, exist_ok=True)
    
    def process_single_sample(self, bam_file: str, sample_id: str) -> Dict:
        """Process a single BAM file for SMN CNV detection"""
        
        self.logger.info(f"Processing sample: {sample_id}")
        
        try:
            # Step 1: Extract depth and quality metrics
            depth_metrics = self.depth_analyzer.analyze_sample(bam_file, sample_id)
            
            # Step 2: Call CNVs using multiple methods
            cnv_results = self.cnv_caller.call_cnvs(depth_metrics)
            
            # Step 3: Generate IGV snapshots
            igv_snapshots = self.igv_automation.generate_snapshots(bam_file, sample_id)
            
            # Step 4: Compile results
            sample_result = {
                'sample_id': sample_id,
                'bam_file': bam_file,
                'depth_metrics': depth_metrics,
                'cnv_calls': cnv_results,
                'igv_snapshots': igv_snapshots,
                'processing_time': datetime.now().isoformat(),
                'pipeline_version': self.config['pipeline_version']
            }
            
            # Step 5: Generate individual sample report
            self.report_generator.generate_sample_report(sample_result)
            
            self.logger.info(f"Successfully processed sample: {sample_id}")
            return sample_result
            
        except Exception as e:
            self.logger.error(f"Error processing sample {sample_id}: {str(e)}")
            return {
                'sample_id': sample_id,
                'error': str(e),
                'status': 'failed'
            }
    
    def process_batch(self, bam_files: List[str], sample_ids: List[str] = None) -> List[Dict]:
        """Process multiple BAM files in batch"""
        
        if sample_ids is None:
            sample_ids = [Path(bam).stem for bam in bam_files]
        
        self.logger.info(f"Starting batch processing of {len(bam_files)} samples")
        
        batch_results = []
        for bam_file, sample_id in zip(bam_files, sample_ids):
            result = self.process_single_sample(bam_file, sample_id)
            batch_results.append(result)
            self.results.append(result)
        
        # Update ML thresholds after batch processing
        self.update_ml_thresholds(batch_results)
        
        # Generate consolidated batch report
        self.report_generator.generate_batch_report(batch_results)
        
        # Save results to TSV
        self.save_batch_tsv(batch_results)
        
        self.logger.info("Batch processing completed")
        return batch_results
    
    def update_ml_thresholds(self, batch_results: List[Dict]):
        """Update ML thresholds based on new results"""
        try:
            # Extract features for threshold optimization
            valid_results = [r for r in batch_results if 'error' not in r]
            
            if len(valid_results) > 0:
                self.ml_optimizer.update_thresholds(valid_results)
                self.logger.info("ML thresholds updated successfully")
            else:
                self.logger.warning("No valid results for threshold update")
                
        except Exception as e:
            self.logger.error(f"Error updating ML thresholds: {str(e)}")
    
    def save_batch_tsv(self, batch_results: List[Dict]):
        """Save batch results to TSV file"""
        
        # Flatten results for TSV format
        tsv_data = []
        for result in batch_results:
            if 'error' not in result:
                row = {
                    'sample_id': result['sample_id'],
                    'smn1_exon7_cnv': result['cnv_calls']['smn1_exon7']['call'],
                    'smn1_exon7_confidence': result['cnv_calls']['smn1_exon7']['confidence'],
                    'smn1_exon8_cnv': result['cnv_calls']['smn1_exon8']['call'],
                    'smn1_exon8_confidence': result['cnv_calls']['smn1_exon8']['confidence'],
                    'smn2_exon7_cnv': result['cnv_calls']['smn2_exon7']['call'],
                    'smn2_exon7_confidence': result['cnv_calls']['smn2_exon7']['confidence'],
                    'smn2_exon8_cnv': result['cnv_calls']['smn2_exon8']['call'],
                    'smn2_exon8_confidence': result['cnv_calls']['smn2_exon8']['confidence'],
                    'sma_risk_prediction': result['cnv_calls']['sma_risk'],
                    'depth_smn1_exon7': result['depth_metrics']['smn1_exon7']['mean_depth'],
                    'depth_smn1_exon8': result['depth_metrics']['smn1_exon8']['mean_depth'],
                    'depth_smn2_exon7': result['depth_metrics']['smn2_exon7']['mean_depth'],
                    'depth_smn2_exon8': result['depth_metrics']['smn2_exon8']['mean_depth'],
                    'mapping_quality': result['depth_metrics']['overall']['mean_mapq'],
                    'gc_content': result['depth_metrics']['overall']['gc_content'],
                    'processing_time': result['processing_time']
                }
                tsv_data.append(row)
            else:
                # Add failed samples
                row = {
                    'sample_id': result['sample_id'],
                    'error': result.get('error', 'Unknown error')
                }
                tsv_data.append(row)
        
        # Save to TSV
        df = pd.DataFrame(tsv_data)
        output_file = os.path.join(self.config['output_dir'], f"batch_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv")
        df.to_csv(output_file, sep='\t', index=False)
        
        self.logger.info(f"Batch results saved to: {output_file}")
        return output_file

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='SMN CNV Detection Pipeline')
    parser.add_argument('-c', '--config', required=True, help='Configuration file (JSON)')
    parser.add_argument('-i', '--input', required=True, help='Input BAM file or directory containing BAM files')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-s', '--sample-id', help='Sample ID (for single file processing)')
    parser.add_argument('--batch', action='store_true', help='Process all BAM files in directory')
    parser.add_argument('--update-thresholds', action='store_true', help='Force threshold update')
    
    args = parser.parse_args()
    
    # Load and update config with command line arguments
    with open(args.config, 'r') as f:
        config = json.load(f)
    
    config['output_dir'] = args.output
    
    # Initialize pipeline
    pipeline = SMNPipeline(args.config)
    
    if args.batch:
        # Process directory of BAM files
        bam_dir = Path(args.input)
        bam_files = list(bam_dir.glob("*.bam"))
        
        if not bam_files:
            print(f"No BAM files found in {args.input}")
            sys.exit(1)
        
        pipeline.process_batch([str(f) for f in bam_files])
        
    else:
        # Process single BAM file
        if not os.path.exists(args.input):
            print(f"Input file not found: {args.input}")
            sys.exit(1)
        
        sample_id = args.sample_id or Path(args.input).stem
        pipeline.process_single_sample(args.input, sample_id)
    
    print("Pipeline execution completed successfully!")

if __name__ == "__main__":
    main()
