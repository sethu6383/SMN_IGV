#!/usr/bin/env python3
"""
CNV Caller Module for SMN Pipeline
Implements multiple computational methods for CNV detection
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
import logging
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import IsolationForest
import pickle
import os

class CNVCaller:
    """Main CNV calling class with multiple methods"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Load thresholds if they exist
        self.thresholds = self._load_thresholds()
        
        # Initialize methods
        self.methods = {
            'depth_ratio': DepthRatioMethod(config),
            'statistical': StatisticalMethod(config),
            'hmm_like': HMMlikeMethod(config),
            'ensemble': EnsembleMethod(config)
        }
    
    def _load_thresholds(self) -> Dict:
        """Load adaptive thresholds from file"""
        threshold_file = os.path.join(self.config.get('model_dir', 'models'), 'thresholds.pkl')
        
        if os.path.exists(threshold_file):
            try:
                with open(threshold_file, 'rb') as f:
                    return pickle.load(f)
            except Exception as e:
                self.logger.warning(f"Could not load thresholds: {e}")
        
        # Default thresholds
        return {
            'depth_ratio': {
                'homo_del_threshold': 0.3,
                'hetero_del_threshold': 0.7,
                'dup_threshold': 1.5
            },
            'statistical': {
                'zscore_threshold': 2.0,
                'pvalue_threshold': 0.05
            },
            'ensemble': {
                'confidence_threshold': 0.8
            }
        }
    
    def call_cnvs(self, depth_metrics: Dict) -> Dict:
        """Main CNV calling function"""
        
        self.logger.info("Calling CNVs using multiple methods")
        
        results = {}
        
        # Call CNVs for each SMN region
        for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
            region_results = {}
            
            # Apply each method
            for method_name, method in self.methods.items():
                try:
                    call_result = method.call_cnv(region, depth_metrics, self.thresholds)
                    region_results[method_name] = call_result
                except Exception as e:
                    self.logger.error(f"Error in {method_name} for {region}: {e}")
                    region_results[method_name] = {
                        'call': 'UNKNOWN',
                        'confidence': 0.0,
                        'error': str(e)
                    }
            
            # Consensus calling
            consensus = self._make_consensus_call(region_results)
            region_results['consensus'] = consensus
            
            results[region] = region_results
        
        # Add SMA risk prediction
        results['sma_risk'] = self._predict_sma_risk(results)
        
        # Flatten for main output
        flattened_results = {}
        for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
            flattened_results[region] = {
                'call': results[region]['consensus']['call'],
                'confidence': results[region]['consensus']['confidence'],
                'methods': results[region]
            }
        
        flattened_results['sma_risk'] = results['sma_risk']
        
        return flattened_results
    
    def _make_consensus_call(self, method_results: Dict) -> Dict:
        """Make consensus call from multiple methods"""
        
        calls = []
        confidences = []
        
        for method_name, result in method_results.items():
            if 'error' not in result:
                calls.append(result['call'])
                confidences.append(result['confidence'])
        
        if not calls:
            return {
                'call': 'UNKNOWN',
                'confidence': 0.0,
                'method': 'consensus'
            }
        
        # Simple majority voting
        from collections import Counter
        call_counts = Counter(calls)
        consensus_call = call_counts.most_common(1)[0][0]
        
        # Average confidence for consensus call
        consensus_indices = [i for i, call in enumerate(calls) if call == consensus_call]
        consensus_confidence = np.mean([confidences[i] for i in consensus_indices])
        
        return {
            'call': consensus_call,
            'confidence': consensus_confidence,
            'method': 'consensus',
            'vote_distribution': dict(call_counts)
        }
    
    def _predict_sma_risk(self, cnv_results: Dict) -> str:
        """Predict SMA risk based on SMN1/SMN2 status"""
        
        smn1_ex7_call = cnv_results['smn1_exon7']['consensus']['call']
        smn1_ex8_call = cnv_results['smn1_exon8']['consensus']['call']
        
        # SMA risk logic
        if smn1_ex7_call == 'HOMO_DEL' and smn1_ex8_call == 'HOMO_DEL':
            return 'HIGH_RISK'
        elif smn1_ex7_call == 'HETERO_DEL' or smn1_ex8_call == 'HETERO_DEL':
            return 'CARRIER'
        elif smn1_ex7_call == 'NORMAL' and smn1_ex8_call == 'NORMAL':
            return 'LOW_RISK'
        else:
            return 'UNCERTAIN'

class DepthRatioMethod:
    """Depth ratio-based CNV calling"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def call_cnv(self, region: str, depth_metrics: Dict, thresholds: Dict) -> Dict:
        """Call CNV based on normalized depth ratios"""
        
        region_metrics = depth_metrics[region]
        normalized_depth = depth_metrics['normalized_depths'][region]['normalized_depth']
        
        # Get thresholds
        thresh = thresholds['depth_ratio']
        homo_del_thresh = thresh['homo_del_threshold']
        hetero_del_thresh = thresh['hetero_del_threshold']
        dup_thresh = thresh['dup_threshold']
        
        # Make call
        if normalized_depth < homo_del_thresh:
            call = 'HOMO_DEL'
            confidence = 1.0 - (normalized_depth / homo_del_thresh)
        elif normalized_depth < hetero_del_thresh:
            call = 'HETERO_DEL'
            confidence = 1.0 - abs(normalized_depth - 0.5) / 0.2
        elif normalized_depth > dup_thresh:
            call = 'DUP'
            confidence = min(1.0, (normalized_depth - dup_thresh) / dup_thresh)
        else:
            call = 'NORMAL'
            confidence = 1.0 - abs(normalized_depth - 1.0) / 0.3
        
        confidence = max(0.0, min(1.0, confidence))
        
        return {
            'call': call,
            'confidence': confidence,
            'normalized_depth': normalized_depth,
            'method': 'depth_ratio'
        }

class StatisticalMethod:
    """Statistical-based CNV calling using z-scores and hypothesis testing"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def call_cnv(self, region: str, depth_metrics: Dict, thresholds: Dict) -> Dict:
        """Call CNV based on statistical analysis"""
        
        region_metrics = depth_metrics[region]
        zscore = depth_metrics['normalized_depths'][region]['zscore']
        
        # Get thresholds
        thresh = thresholds['statistical']
        zscore_thresh = thresh['zscore_threshold']
        
        # Statistical test (simplified)
        pvalue = 2 * (1 - stats.norm.cdf(abs(zscore)))
        
        # Make call based on z-score
        if zscore < -zscore_thresh:
            if zscore < -3:
                call = 'HOMO_DEL'
                confidence = min(1.0, abs(zscore) / 3.0)
            else:
                call = 'HETERO_DEL'
                confidence = min(1.0, abs(zscore) / zscore_thresh)
        elif zscore > zscore_thresh:
            call = 'DUP'
            confidence = min(1.0, zscore / zscore_thresh)
        else:
            call = 'NORMAL'
            confidence = 1.0 - (abs(zscore) / zscore_thresh)
        
        confidence = max(0.0, min(1.0, confidence))
        
        return {
            'call': call,
            'confidence': confidence,
            'zscore': zscore,
            'pvalue': pvalue,
            'method': 'statistical'
        }

class HMMlikeMethod:
    """HMM-like method for CNV calling using coverage patterns"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def call_cnv(self, region: str, depth_metrics: Dict, thresholds: Dict) -> Dict:
        """Call CNV using HMM-like approach"""
        
        region_metrics = depth_metrics[region]
        
        # Get depth statistics
        mean_depth = region_metrics['mean_depth']
        std_depth = region_metrics['std_depth']
        coverage_uniformity = region_metrics['coverage_uniformity']
        zero_depth_positions = region_metrics['zero_depth_positions']
        region_length = region_metrics['region_length']
        
        # Calculate features
        cv = std_depth / (mean_depth + 1e-8)  # Coefficient of variation
        zero_fraction = zero_depth_positions / region_length
        
        # Simple state classification
        if zero_fraction > 0.8:  # High fraction of zero coverage
            call = 'HOMO_DEL'
            confidence = zero_fraction
        elif zero_fraction > 0.3 or cv > 2.0:  # Patchy coverage
            call = 'HETERO_DEL'
            confidence = max(zero_fraction, cv / 3.0)
        elif cv < 0.3 and coverage_uniformity > 0.8:  # Very uniform high coverage
            normalized_depth = depth_metrics['normalized_depths'][region]['normalized_depth']
            if normalized_depth > 1.5:
                call = 'DUP'
                confidence = min(1.0, normalized_depth - 1.0)
            else:
                call = 'NORMAL'
                confidence = coverage_uniformity
        else:
            call = 'NORMAL'
            confidence = 1.0 - cv
        
        confidence = max(0.0, min(1.0, confidence))
        
        return {
            'call': call,
            'confidence': confidence,
            'cv': cv,
            'zero_fraction': zero_fraction,
            'coverage_uniformity': coverage_uniformity,
            'method': 'hmm_like'
        }

class EnsembleMethod:
    """Ensemble method combining multiple features"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Try to load trained model
        self.model = self._load_model()
    
    def _load_model(self):
        """Load pre-trained ensemble model"""
        model_file = os.path.join(self.config.get('model_dir', 'models'), 'ensemble_model.pkl')
        
        if os.path.exists(model_file):
            try:
                with open(model_file, 'rb') as f:
                    return pickle.load(f)
            except Exception as e:
                self.logger.warning(f"Could not load ensemble model: {e}")
        
        return None
    
    def call_cnv(self, region: str, depth_metrics: Dict, thresholds: Dict) -> Dict:
        """Call CNV using ensemble approach"""
        
        # Extract features
        features = self._extract_features(region, depth_metrics)
        
        if self.model is None:
            # Fallback to rule-based ensemble
            return self._rule_based_ensemble(features, thresholds)
        else:
            # Use trained model
            return self._model_based_ensemble(features)
    
    def _extract_features(self, region: str, depth_metrics: Dict) -> np.ndarray:
        """Extract features for ensemble method"""
        
        region_metrics = depth_metrics[region]
        normalized = depth_metrics['normalized_depths'][region]
        overall = depth_metrics['overall']
        
        features = [
            normalized['normalized_depth'],
            normalized['zscore'],
            region_metrics['coverage_uniformity'],
            region_metrics['zero_depth_positions'] / region_metrics['region_length'],
            region_metrics['std_depth'] / (region_metrics['mean_depth'] + 1e-8),
            region_metrics['mean_mapq'],
            overall['mapping_rate'],
            region_metrics['gc_content']
        ]
        
        return np.array(features)
    
    def _rule_based_ensemble(self, features: np.ndarray, thresholds: Dict) -> Dict:
        """Rule-based ensemble when no model is available"""
        
        normalized_depth = features[0]
        zscore = features[1]
        coverage_uniformity = features[2]
        zero_fraction = features[3]
        cv = features[4]
        
        # Scoring system
        scores = {
            'HOMO_DEL': 0,
            'HETERO_DEL': 0,
            'NORMAL': 0,
            'DUP': 0
        }
        
        # Depth-based scoring
        if normalized_depth < 0.3:
            scores['HOMO_DEL'] += 3
        elif normalized_depth < 0.7:
            scores['HETERO_DEL'] += 2
        elif normalized_depth > 1.5:
            scores['DUP'] += 2
        else:
            scores['NORMAL'] += 1
        
        # Z-score based scoring
        if zscore < -2:
            scores['HOMO_DEL'] += 2
        elif zscore < -1:
            scores['HETERO_DEL'] += 1
        elif zscore > 2:
            scores['DUP'] += 1
        else:
            scores['NORMAL'] += 1
        
        # Coverage pattern scoring
        if zero_fraction > 0.5:
            scores['HOMO_DEL'] += 2
        elif zero_fraction > 0.2:
            scores['HETERO_DEL'] += 1
        
        if cv > 1.5:
            scores['HETERO_DEL'] += 1
        
        if coverage_uniformity > 0.8:
            scores['NORMAL'] += 1
            scores['DUP'] += 1
        
        # Make call
        max_score = max(scores.values())
        best_calls = [call for call, score in scores.items() if score == max_score]
        
        # If tie, prefer more conservative call
        call_priority = ['NORMAL', 'HETERO_DEL', 'DUP', 'HOMO_DEL']
        for preferred_call in call_priority:
            if preferred_call in best_calls:
                final_call = preferred_call
                break
        else:
            final_call = best_calls[0]
        
        confidence = max_score / 6.0  # Normalize by max possible score
        confidence = max(0.0, min(1.0, confidence))
        
        return {
            'call': final_call,
            'confidence': confidence,
            'scores': scores,
            'method': 'ensemble_rule_based'
        }
    
    def _model_based_ensemble(self, features: np.ndarray) -> Dict:
        """Model-based ensemble prediction"""
        
        try:
            prediction = self.model.predict([features])[0]
            confidence = max(self.model.predict_proba([features])[0])
            
            return {
                'call': prediction,
                'confidence': confidence,
                'method': 'ensemble_model_based'
            }
        except Exception as e:
            self.logger.error(f"Error in model prediction: {e}")
            # Fallback to simple depth ratio
            normalized_depth = features[0]
            if normalized_depth < 0.3:
                call = 'HOMO_DEL'
            elif normalized_depth < 0.7:
                call = 'HETERO_DEL'
            elif normalized_depth > 1.5:
                call = 'DUP'
            else:
                call = 'NORMAL'
            
            return {
                'call': call,
                'confidence': 0.5,
                'method': 'ensemble_fallback'
            }
