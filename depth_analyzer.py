#!/usr/bin/env python3
"""
Depth Analysis Module for SMN CNV Pipeline
Analyzes read depth and quality metrics for SMN1/SMN2 regions
"""

import pysam
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
import logging
from scipy import stats
from collections import defaultdict

class DepthAnalyzer:
    """Analyzes read depth and quality metrics for SMN regions"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # SMN gene coordinates (hg38) - Update these with your specific coordinates
        self.smn_regions = {
            'smn1_exon7': {
                'chr': 'chr5',
                'start': 70247724,
                'end': 70247775,
                'strand': '-'
            },
            'smn1_exon8': {
                'chr': 'chr5', 
                'start': 70248935,
                'end': 70249306,
                'strand': '-'
            },
            'smn2_exon7': {
                'chr': 'chr5',
                'start': 69372304,
                'end': 69372355,
                'strand': '-'
            },
            'smn2_exon8': {
                'chr': 'chr5',
                'start': 69373515,
                'end': 69373886,
                'strand': '-'
            }
        }
        
        # Control regions for normalization
        self.control_regions = self._get_control_regions()
    
    def _get_control_regions(self) -> Dict:
        """Define control regions for depth normalization"""
        return {
            'actb_exon1': {'chr': 'chr7', 'start': 5527147, 'end': 5527252},
            'gapdh_exon1': {'chr': 'chr12', 'start': 6534405, 'end': 6534516},
            'tbp_exon1': {'chr': 'chr6', 'start': 170554001, 'end': 170554099}
        }
    
    def analyze_sample(self, bam_file: str, sample_id: str) -> Dict:
        """Analyze depth metrics for a single sample"""
        
        self.logger.info(f"Analyzing depth for sample: {sample_id}")
        
        try:
            with pysam.AlignmentFile(bam_file, 'rb') as bam:
                
                # Analyze SMN regions
                smn_metrics = {}
                for region_name, coords in self.smn_regions.items():
                    metrics = self._analyze_region(bam, region_name, coords)
                    smn_metrics[region_name] = metrics
                
                # Analyze control regions
                control_metrics = {}
                for region_name, coords in self.control_regions.items():
                    metrics = self._analyze_region(bam, region_name, coords)
                    control_metrics[region_name] = metrics
                
                # Calculate normalized depths
                normalized_depths = self._calculate_normalized_depths(smn_metrics, control_metrics)
                
                # Calculate overall sample metrics
                overall_metrics = self._calculate_overall_metrics(bam)
                
                return {
                    'sample_id': sample_id,
                    **smn_metrics,
                    'controls': control_metrics,
                    'normalized_depths': normalized_depths,
                    'overall': overall_metrics
                }
                
        except Exception as e:
            self.logger.error(f"Error analyzing sample {sample_id}: {str(e)}")
            raise
    
    def _analyze_region(self, bam: pysam.AlignmentFile, region_name: str, coords: Dict) -> Dict:
        """Analyze a specific genomic region"""
        
        chr_name = coords['chr']
        start = coords['start']
        end = coords['end']
        
        # Extract reads in region
        depths = []
        mapping_qualities = []
        base_qualities = []
        gc_content_values = []
        
        # Get pileup columns
        for pileup_column in bam.pileup(chr_name, start, end, truncate=True):
            pos = pileup_column.pos
            if start <= pos < end:
                
                depth = pileup_column.n
                depths.append(depth)
                
                # Calculate mapping quality for reads at this position
                mapqs = []
                baseqs = []
                
                for pileup_read in pileup_column.pileups:
                    if not pileup_read.is_del and not pileup_read.is_refskip:
                        mapqs.append(pileup_read.alignment.mapping_quality)
                        if pileup_read.query_position is not None:
                            baseqs.append(pileup_read.alignment.query_qualities[pileup_read.query_position])
                
                if mapqs:
                    mapping_qualities.append(np.mean(mapqs))
                if baseqs:
                    base_qualities.append(np.mean(baseqs))
        
        # Calculate GC content
        try:
            gc_content = self._calculate_gc_content(bam, chr_name, start, end)
        except:
            gc_content = 0.5  # Default GC content
        
        # Calculate statistics
        if depths:
            depth_stats = {
                'mean_depth': np.mean(depths),
                'median_depth': np.median(depths),
                'std_depth': np.std(depths),
                'min_depth': np.min(depths),
                'max_depth': np.max(depths),
                'zero_depth_positions': sum(1 for d in depths if d == 0),
                'low_depth_positions': sum(1 for d in depths if d < 10)  # Configurable threshold
            }
        else:
            depth_stats = {
                'mean_depth': 0,
                'median_depth': 0,
                'std_depth': 0,
                'min_depth': 0,
                'max_depth': 0,
                'zero_depth_positions': end - start,
                'low_depth_positions': end - start
            }
        
        quality_stats = {
            'mean_mapq': np.mean(mapping_qualities) if mapping_qualities else 0,
            'mean_baseq': np.mean(base_qualities) if base_qualities else 0,
            'gc_content': gc_content
        }
        
        # Additional metrics
        region_length = end - start
        coverage_uniformity = 1 - (depth_stats['std_depth'] / (depth_stats['mean_depth'] + 1e-8))
        
        return {
            **depth_stats,
            **quality_stats,
            'region_length': region_length,
            'coverage_uniformity': coverage_uniformity,
            'coordinates': f"{chr_name}:{start}-{end}"
        }
    
    def _calculate_gc_content(self, bam: pysam.AlignmentFile, chr_name: str, start: int, end: int) -> float:
        """Calculate GC content for a region (simplified version)"""
        # This is a placeholder - in practice, you'd want to use a reference genome
        # For now, return a default value
        return 0.5
    
    def _calculate_normalized_depths(self, smn_metrics: Dict, control_metrics: Dict) -> Dict:
        """Calculate normalized depths using control regions"""
        
        # Calculate median depth across control regions
        control_depths = [metrics['mean_depth'] for metrics in control_metrics.values()]
        median_control_depth = np.median(control_depths) if control_depths else 30.0  # Default expected depth
        
        normalized = {}
        for region_name, metrics in smn_metrics.items():
            normalized[region_name] = {
                'raw_depth': metrics['mean_depth'],
                'normalized_depth': metrics['mean_depth'] / (median_control_depth + 1e-8),
                'depth_ratio': metrics['mean_depth'] / (median_control_depth + 1e-8),
                'zscore': self._calculate_zscore(metrics['mean_depth'], control_depths)
            }
        
        return normalized
    
    def _calculate_zscore(self, value: float, control_values: List[float]) -> float:
        """Calculate z-score relative to control values"""
        if not control_values:
            return 0.0
        
        mean_control = np.mean(control_values)
        std_control = np.std(control_values)
        
        if std_control == 0:
            return 0.0
        
        return (value - mean_control) / std_control
    
    def _calculate_overall_metrics(self, bam: pysam.AlignmentFile) -> Dict:
        """Calculate overall sample metrics"""
        
        total_reads = 0
        total_mapped = 0
        mapq_sum = 0
        
        # Sample a subset of reads for efficiency
        sample_size = 10000
        read_count = 0
        
        for read in bam.fetch():
            if read_count >= sample_size:
                break
                
            total_reads += 1
            if not read.is_unmapped:
                total_mapped += 1
                mapq_sum += read.mapping_quality
                
            read_count += 1
        
        mapping_rate = total_mapped / (total_reads + 1e-8)
        mean_mapq = mapq_sum / (total_mapped + 1e-8)
        
        return {
            'total_reads_sampled': total_reads,
            'mapping_rate': mapping_rate,
            'mean_mapq': mean_mapq,
            'gc_content': 0.5  # Placeholder
        }
    
    def calculate_depth_ratios(self, sample_metrics: Dict) -> Dict:
        """Calculate various depth ratios for CNV calling"""
        
        ratios = {}
        
        # SMN1 vs SMN2 ratios
        smn1_ex7_depth = sample_metrics['smn1_exon7']['mean_depth']
        smn1_ex8_depth = sample_metrics['smn1_exon8']['mean_depth']
        smn2_ex7_depth = sample_metrics['smn2_exon7']['mean_depth']
        smn2_ex8_depth = sample_metrics['smn2_exon8']['mean_depth']
        
        ratios['smn1_smn2_ex7'] = smn1_ex7_depth / (smn2_ex7_depth + 1e-8)
        ratios['smn1_smn2_ex8'] = smn1_ex8_depth / (smn2_ex8_depth + 1e-8)
        ratios['ex7_ex8_smn1'] = smn1_ex7_depth / (smn1_ex8_depth + 1e-8)
        ratios['ex7_ex8_smn2'] = smn2_ex7_depth / (smn2_ex8_depth + 1e-8)
        
        # Total SMN depth
        total_smn_depth = smn1_ex7_depth + smn1_ex8_depth + smn2_ex7_depth + smn2_ex8_depth
        ratios['total_smn_depth'] = total_smn_depth
        
        return ratios
