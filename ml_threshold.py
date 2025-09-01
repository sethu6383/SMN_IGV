#!/usr/bin/env python3
"""
ML Threshold Optimizer for SMN Pipeline
Adapts thresholds based on performance data
"""

import numpy as np
import pandas as pd
import pickle
import os
import json
import logging
from typing import Dict, List, Tuple, Optional
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.calibration import CalibratedClassifierCV
import optuna
from datetime import datetime

class MLThresholdOptimizer:
    """ML-based threshold optimization system"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Model directory
        self.model_dir = Path(config.get('model_dir', 'models'))
        self.model_dir.mkdir(parents=True, exist_ok=True)
        
        # Files for persistent storage
        self.threshold_file = self.model_dir / 'thresholds.pkl'
        self.performance_file = self.model_dir / 'performance_history.json'
        self.model_file = self.model_dir / 'ensemble_model.pkl'
        self.scaler_file = self.model_dir / 'feature_scaler.pkl'
        
        # Load existing data
        self.thresholds = self._load_thresholds()
        self.performance_history = self._load_performance_history()
        self.trained_models = self._load_trained_models()
        
        # Training data accumulator
        self.training_data = []
        self.validation_data = []
        
        # Load MLPA validated samples if available
        self._load_validation_data()
    
    def _load_thresholds(self) -> Dict:
        """Load current thresholds"""
        if self.threshold_file.exists():
            try:
                with open(self.threshold_file, 'rb') as f:
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
            },
            'last_updated': datetime.now().isoformat(),
            'samples_processed': 0,
            'performance_metrics': {}
        }
    
    def _load_performance_history(self) -> List[Dict]:
        """Load performance history"""
        if self.performance_file.exists():
            try:
                with open(self.performance_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                self.logger.warning(f"Could not load performance history: {e}")
        return []
    
    def _load_trained_models(self) -> Dict:
        """Load trained models"""
        models = {}
        
        if self.model_file.exists():
            try:
                with open(self.model_file, 'rb') as f:
                    models['classifier'] = pickle.load(f)
            except Exception as e:
                self.logger.warning(f"Could not load trained model: {e}")
        
        if self.scaler_file.exists():
            try:
                with open(self.scaler_file, 'rb') as f:
                    models['scaler'] = pickle.load(f)
            except Exception as e:
                self.logger.warning(f"Could not load scaler: {e}")
        
        return models
    
    def _load_validation_data(self):
        """Load MLPA validated samples"""
        validation_file = self.config.get('validation_data_file')
        
        if validation_file and os.path.exists(validation_file):
            try:
                df = pd.read_csv(validation_file, sep='\t')
                self.validation_data = df.to_dict('records')
                self.logger.info(f"Loaded {len(self.validation_data)} validation samples")
            except Exception as e:
                self.logger.error(f"Error loading validation data: {e}")
    
    def update_thresholds(self, batch_results: List[Dict]):
        """Update thresholds based on new batch results"""
        
        self.logger.info("Updating ML thresholds based on new results")
        
        # Add to training data
        for result in batch_results:
            if 'error' not in result:
                self.training_data.append(result)
        
        # Update sample count
        self.thresholds['samples_processed'] += len(batch_results)
        
        # If we have enough samples, retrain models
        if len(self.training_data) >= 20:  # Minimum samples for retraining
            self._retrain_models()
        
        # Optimize thresholds using current data
        self._optimize_thresholds()
        
        # Save updated thresholds
        self._save_thresholds()
        
        # Record performance
        self._record_performance()
    
    def _retrain_models(self):
        """Retrain ML models with accumulated data"""
        
        try:
            # Prepare training data
            features, labels = self._prepare_training_data()
            
            if len(features) < 10:  # Not enough data
                self.logger.warning("Insufficient data for model retraining")
                return
            
            # Split data
            from sklearn.model_selection import train_test_split
            X_train, X_test, y_train, y_test = train_test_split(
                features, labels, test_size=0.2, random_state=42, stratify=labels
            )
            
            # Scale features
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)
            
            # Train model with hyperparameter optimization
            best_model = self._optimize_model_hyperparameters(X_train_scaled, y_train)
            
            # Evaluate model
            y_pred = best_model.predict(X_test_scaled)
            accuracy = accuracy_score(y_test, y_pred)
            
            self.logger.info(f"Model retrained with accuracy: {accuracy:.3f}")
            
            # Save models if performance is good
            if accuracy > 0.7:  # Minimum acceptable accuracy
                self.trained_models['classifier'] = best_model
                self.trained_models['scaler'] = scaler
                self._save_trained_models()
                
                # Clear training data to prevent memory issues
                self.training_data = self.training_data[-100:]  # Keep last 100 samples
            
        except Exception as e:
            self.logger.error(f"Error retraining models: {e}")
    
    def _prepare_training_data(self) -> Tuple[np.ndarray, np.ndarray]:
        """Prepare training data from accumulated results"""
        
        features = []
        labels = []
        
        # Combine training data with validation data
        all_data = self.training_data + self.validation_data
        
        for sample in all_data:
            try:
                # Extract features for each region
                for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
                    if region in sample.get('depth_metrics', {}):
                        sample_features = self._extract_sample_features(sample, region)
                        sample_label = self._get_sample_label(sample, region)
                        
                        if sample_features is not None and sample_label is not None:
                            features.append(sample_features)
                            labels.append(sample_label)
                            
            except Exception as e:
                self.logger.warning(f"Error processing sample for training: {e}")
                continue
        
        return np.array(features), np.array(labels)
    
    def _extract_sample_features(self, sample: Dict, region: str) -> Optional[np.ndarray]:
        """Extract features for ML model"""
        
        try:
            depth_metrics = sample['depth_metrics']
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
                region_metrics['mean_baseq'],
                overall['mapping_rate'],
                region_metrics['gc_content'],
                region_metrics['min_depth'],
                region_metrics['max_depth']
            ]
            
            return np.array(features)
            
        except Exception as e:
            return None
    
    def _get_sample_label(self, sample: Dict, region: str) -> Optional[str]:
        """Get ground truth label for sample"""
        
        # Check if this is validation data with known labels
        if 'mlpa_result' in sample:
            return sample['mlpa_result'].get(region)
        
        # Check if we have consensus call (for training on our own calls)
        if 'cnv_calls' in sample and region in sample['cnv_calls']:
            return sample['cnv_calls'][region]['call']
        
        return None
    
    def _optimize_model_hyperparameters(self, X_train: np.ndarray, y_train: np.ndarray):
        """Optimize model hyperparameters using Optuna"""
        
        def objective(trial):
            # Suggest hyperparameters
            n_estimators = trial.suggest_int('n_estimators', 50, 300)
            max_depth = trial.suggest_int('max_depth', 3, 20)
            min_samples_split = trial.suggest_int('min_samples_split', 2, 20)
            min_samples_leaf = trial.suggest_int('min_samples_leaf', 1, 10)
            
            # Create and train model
            model = RandomForestClassifier(
                n_estimators=n_estimators,
                max_depth=max_depth,
                min_samples_split=min_samples_split,
                min_samples_leaf=min_samples_leaf,
                random_state=42,
                n_jobs=-1
            )
            
            # Cross-validation score
            from sklearn.model_selection import cross_val_score
            scores = cross_val_score(model, X_train, y_train, cv=3, scoring='accuracy')
            return scores.mean()
        
        # Optimize hyperparameters
        study = optuna.create_study(direction='maximize')
        study.optimize(objective, n_trials=50, timeout=300)  # 5 min timeout
        
        # Train final model with best parameters
        best_params = study.best_params
        final_model = RandomForestClassifier(**best_params, random_state=42, n_jobs=-1)
        final_model.fit(X_train, y_train)
        
        # Calibrate probabilities
        calibrated_model = CalibratedClassifierCV(final_model, cv=3)
        calibrated_model.fit(X_train, y_train)
        
        return calibrated_model
    
    def _optimize_thresholds(self):
        """Optimize thresholds using Optuna"""
        
        if len(self.validation_data) < 10:  # Need validation data for optimization
            self.logger.warning("Insufficient validation data for threshold optimization")
            return
        
        def objective(trial):
            # Suggest threshold values
            homo_del_threshold = trial.suggest_float('homo_del_threshold', 0.1, 0.5)
            hetero_del_threshold = trial.suggest_float('hetero_del_threshold', 0.5, 0.9)
            dup_threshold = trial.suggest_float('dup_threshold', 1.2, 2.0)
            zscore_threshold = trial.suggest_float('zscore_threshold', 1.5, 3.0)
            confidence_threshold = trial.suggest_float('confidence_threshold', 0.6, 0.95)
            
            # Test thresholds on validation data
            test_thresholds = {
                'depth_ratio': {
                    'homo_del_threshold': homo_del_threshold,
                    'hetero_del_threshold': hetero_del_threshold,
                    'dup_threshold': dup_threshold
                },
                'statistical': {
                    'zscore_threshold': zscore_threshold,
                    'pvalue_threshold': 0.05
                },
                'ensemble': {
                    'confidence_threshold': confidence_threshold
                }
            }
            
            accuracy = self._evaluate_thresholds(test_thresholds)
            return accuracy
        
        # Optimize thresholds
        study = optuna.create_study(direction='maximize')
        study.optimize(objective, n_trials=100, timeout=600)  # 10 min timeout
        
        # Update thresholds with best values
        best_params = study.best_params
        self.thresholds['depth_ratio'].update({
            'homo_del_threshold': best_params['homo_del_threshold'],
            'hetero_del_threshold': best_params['hetero_del_threshold'],
            'dup_threshold': best_params['dup_threshold']
        })
        self.thresholds['statistical']['zscore_threshold'] = best_params['zscore_threshold']
        self.thresholds['ensemble']['confidence_threshold'] = best_params['confidence_threshold']
        
        self.logger.info(f"Thresholds optimized with accuracy: {study.best_value:.3f}")
    
    def _evaluate_thresholds(self, test_thresholds: Dict) -> float:
        """Evaluate threshold performance on validation data"""
        
        correct_predictions = 0
        total_predictions = 0
        
        for sample in self.validation_data:
            try:
                for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
                    if region in sample.get('depth_metrics', {}) and region in sample.get('mlpa_result', {}):
                        
                        # Make prediction using test thresholds
                        prediction = self._make_threshold_prediction(sample, region, test_thresholds)
                        ground_truth = sample['mlpa_result'][region]
                        
                        if prediction == ground_truth:
                            correct_predictions += 1
                        total_predictions += 1
                        
            except Exception as e:
                continue
        
        if total_predictions == 0:
            return 0.0
        
        return correct_predictions / total_predictions
    
    def _make_threshold_prediction(self, sample: Dict, region: str, thresholds: Dict) -> str:
        """Make CNV prediction using given thresholds"""
        
        try:
            normalized_depth = sample['depth_metrics']['normalized_depths'][region]['normalized_depth']
            zscore = sample['depth_metrics']['normalized_depths'][region]['zscore']
            
            # Depth ratio method
            depth_thresholds = thresholds['depth_ratio']
            if normalized_depth < depth_thresholds['homo_del_threshold']:
                depth_call = 'HOMO_DEL'
            elif normalized_depth < depth_thresholds['hetero_del_threshold']:
                depth_call = 'HETERO_DEL'
            elif normalized_depth > depth_thresholds['dup_threshold']:
                depth_call = 'DUP'
            else:
                depth_call = 'NORMAL'
            
            # Statistical method
            zscore_threshold = thresholds['statistical']['zscore_threshold']
            if zscore < -zscore_threshold:
                if zscore < -3:
                    stat_call = 'HOMO_DEL'
                else:
                    stat_call = 'HETERO_DEL'
            elif zscore > zscore_threshold:
                stat_call = 'DUP'
            else:
                stat_call = 'NORMAL'
            
            # Simple consensus (could be improved)
            if depth_call == stat_call:
                return depth_call
            else:
                # If disagreement, prefer depth ratio for SMN
                return depth_call
                
        except Exception as e:
            return 'NORMAL'  # Default fallback
    
    def _save_thresholds(self):
        """Save current thresholds to file"""
        try:
            self.thresholds['last_updated'] = datetime.now().isoformat()
            with open(self.threshold_file, 'wb') as f:
                pickle.dump(self.thresholds, f)
            self.logger.info("Thresholds saved successfully")
        except Exception as e:
            self.logger.error(f"Error saving thresholds: {e}")
    
    def _save_trained_models(self):
        """Save trained models to file"""
        try:
            if 'classifier' in self.trained_models:
                with open(self.model_file, 'wb') as f:
                    pickle.dump(self.trained_models['classifier'], f)
            
            if 'scaler' in self.trained_models:
                with open(self.scaler_file, 'wb') as f:
                    pickle.dump(self.trained_models['scaler'], f)
                    
            self.logger.info("Models saved successfully")
        except Exception as e:
            self.logger.error(f"Error saving models: {e}")
    
    def _record_performance(self):
        """Record current performance metrics"""
        
        if not self.validation_data:
            return
        
        # Evaluate current performance
        accuracy = self._evaluate_thresholds(self.thresholds)
        
        performance_record = {
            'timestamp': datetime.now().isoformat(),
            'samples_processed': self.thresholds['samples_processed'],
            'accuracy': accuracy,
            'thresholds': self.thresholds.copy()
        }
        
        self.performance_history.append(performance_record)
        
        # Keep only last 100 records
        self.performance_history = self.performance_history[-100:]
        
        # Save performance history
        try:
            with open(self.performance_file, 'w') as f:
                json.dump(self.performance_history, f, indent=2)
        except Exception as e:
            self.logger.error(f"Error saving performance history: {e}")
        
        self.logger.info(f"Current accuracy: {accuracy:.3f}")
    
    def get_performance_summary(self) -> Dict:
        """Get performance summary"""
        
        if not self.performance_history:
            return {'status': 'no_data'}
        
        recent_performance = self.performance_history[-10:]  # Last 10 records
        accuracies = [p['accuracy'] for p in recent_performance]
        
        return {
            'current_accuracy': accuracies[-1] if accuracies else 0.0,
            'mean_accuracy': np.mean(accuracies),
            'std_accuracy': np.std(accuracies),
            'samples_processed': self.thresholds['samples_processed'],
            'last_updated': self.thresholds.get('last_updated'),
            'target_accuracy': 0.95,
            'accuracy_trend': 'improving' if len(accuracies) > 1 and accuracies[-1] > accuracies[0] else 'stable'
        }
    
    def should_retrain(self) -> bool:
        """Determine if models should be retrained"""
        
        # Retrain every 50 samples or if accuracy drops
        samples_since_update = self.thresholds['samples_processed']
        
        if samples_since_update > 0 and samples_since_update % 50 == 0:
            return True
        
        # Check if accuracy is declining
        if len(self.performance_history) > 5:
            recent_acc = [p['accuracy'] for p in self.performance_history[-5:]]
            if len(recent_acc) > 2 and recent_acc[-1] < recent_acc[0] - 0.05:
                return True
        
        return False
    
    def export_model_performance(self) -> str:
        """Export detailed model performance report"""
        
        report_file = self.model_dir / f'performance_report_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
        
        report = {
            'summary': self.get_performance_summary(),
            'thresholds': self.thresholds,
            'performance_history': self.performance_history,
            'validation_samples': len(self.validation_data),
            'training_samples': len(self.training_data)
        }
        
        try:
            with open(report_file, 'w') as f:
                json.dump(report, f, indent=2)
            
            self.logger.info(f"Performance report exported: {report_file}")
            return str(report_file)
            
        except Exception as e:
            self.logger.error(f"Error exporting performance report: {e}")
            return None
