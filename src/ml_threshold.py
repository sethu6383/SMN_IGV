#!/usr/bin/env python3
"""
Complete ML Threshold Optimizer for SMN Pipeline
"""

import numpy as np
import pandas as pd
import pickle
import os
import json
import logging
from typing import Dict, List, Tuple, Optional
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.calibration import CalibratedClassifierCV
import optuna
from datetime import datetime
from pathlib import Path

class MLThresholdOptimizer:
    """Complete ML-based threshold optimization system"""
    
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
        
        # Default thresholds from config
        default_thresholds = self.config.get('cnv_thresholds', {})
        default_thresholds.update({
            'last_updated': datetime.now().isoformat(),
            'samples_processed': 0,
            'performance_metrics': {}
        })
        return default_thresholds
    
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
                    self.logger.info("Loaded trained classifier model")
            except Exception as e:
                self.logger.warning(f"Could not load trained model: {e}")
        
        if self.scaler_file.exists():
            try:
                with open(self.scaler_file, 'rb') as f:
                    models['scaler'] = pickle.load(f)
                    self.logger.info("Loaded feature scaler")
            except Exception as e:
                self.logger.warning(f"Could not load scaler: {e}")
        
        return models
    
    def _load_validation_data(self):
        """Load MLPA validated samples"""
        validation_file = self.config.get('validation_data_file')
        
        if validation_file and os.path.exists(validation_file):
            try:
                df = pd.read_csv(validation_file, sep='\t')
                
                # Convert to list of dictionaries
                self.validation_data = []
                for _, row in df.iterrows():
                    sample_dict = {
                        'sample_id': row['sample_id'],
                        'mlpa_result': {
                            'smn1_exon7': row['smn1_exon7'],
                            'smn1_exon8': row['smn1_exon8'],
                            'smn2_exon7': row['smn2_exon7'],
                            'smn2_exon8': row['smn2_exon8']
                        }
                    }
                    self.validation_data.append(sample_dict)
                
                self.logger.info(f"Loaded {len(self.validation_data)} validation samples")
            except Exception as e:
                self.logger.error(f"Error loading validation data: {e}")
    
    def update_thresholds(self, batch_results: List[Dict]):
        """Update thresholds based on new batch results"""
        
        self.logger.info("Updating ML thresholds based on new results")
        
        # Add to training data
        valid_results = [r for r in batch_results if 'error' not in r]
        self.training_data.extend(valid_results)
        
        # Update sample count
        self.thresholds['samples_processed'] += len(valid_results)
        
        # If we have enough samples, retrain models
        min_samples = self.config.get('ml_optimization', {}).get('min_samples_for_training', 20)
        if len(self.training_data) >= min_samples:
            self._retrain_models()
        
        # Optimize thresholds using current data
        if self.validation_data:
            self._optimize_thresholds()
        
        # Save updated thresholds
        self._save_thresholds()
        
        # Record performance
        self._record_performance()
    
    def _retrain_models(self):
        """Retrain ML models with accumulated data"""
        
        try:
            self.logger.info("Retraining ML models...")
            
            # Prepare training data
            features, labels = self._prepare_training_data()
            
            if len(features) < 10:
                self.logger.warning("Insufficient data for model retraining")
                return
            
            # Check class distribution
            unique_labels, counts = np.unique(labels, return_counts=True)
            self.logger.info(f"Training data distribution: {dict(zip(unique_labels, counts))}")
            
            # Split data
            X_train, X_test, y_train, y_test = train_test_split(
                features, labels, test_size=0.2, random_state=42, 
                stratify=labels if len(unique_labels) > 1 else None
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
            if accuracy > 0.7:
                self.trained_models['classifier'] = best_model
                self.trained_models['scaler'] = scaler
                self._save_trained_models()
                
                # Keep reasonable amount of training data
                if len(self.training_data) > 200:
                    self.training_data = self.training_data[-200:]
                
                self.logger.info("Models saved successfully")
            else:
                self.logger.warning(f"Model accuracy too low ({accuracy:.3f}), not saving")
            
        except Exception as e:
            self.logger.error(f"Error retraining models: {e}")
    
    def _prepare_training_data(self) -> Tuple[np.ndarray, np.ndarray]:
        """Prepare training data from accumulated results and validation data"""
        
        features = []
        labels = []
        
        # Process validation data (ground truth)
        for sample in self.validation_data:
            try:
                # For validation data, we need to simulate depth metrics
                # In practice, you'd have actual depth metrics from processing these samples
                sample_id = sample['sample_id']
                
                # Skip if we don't have processed this validation sample yet
                matching_result = None
                for training_sample in self.training_data:
                    if training_sample.get('sample_id') == sample_id:
                        matching_result = training_sample
                        break
                
                if matching_result:
                    for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
                        sample_features = self._extract_sample_features(matching_result, region)
                        sample_label = sample['mlpa_result'].get(region)
                        
                        if sample_features is not None and sample_label is not None:
                            features.append(sample_features)
                            labels.append(sample_label)
                            
            except Exception as e:
                self.logger.warning(f"Error processing validation sample {sample.get('sample_id', 'unknown')}: {e}")
                continue
        
        # Process training data (our own calls - use high confidence calls only)
        for sample in self.training_data:
            try:
                for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
                    if (region in sample.get('cnv_calls', {}) and 
                        sample['cnv_calls'][region].get('confidence', 0) > 0.9):  # High confidence only
                        
                        sample_features = self._extract_sample_features(sample, region)
                        sample_label = sample['cnv_calls'][region]['call']
                        
                        if sample_features is not None and sample_label is not None:
                            features.append(sample_features)
                            labels.append(sample_label)
                            
            except Exception as e:
                self.logger.warning(f"Error processing training sample: {e}")
                continue
        
        if len(features) == 0:
            raise ValueError("No valid training data available")
        
        return np.array(features), np.array(labels)
    
    def _extract_sample_features(self, sample: Dict, region: str) -> Optional[np.ndarray]:
        """Extract features for ML model"""
        
        try:
            depth_metrics = sample.get('depth_metrics', {})
            if region not in depth_metrics:
                return None
                
            region_metrics = depth_metrics[region]
            normalized = depth_metrics.get('normalized_depths', {}).get(region, {})
            overall = depth_metrics.get('overall', {})
            
            features = [
                normalized.get('normalized_depth', 1.0),
                normalized.get('zscore', 0.0),
                region_metrics.get('coverage_uniformity', 0.5),
                region_metrics.get('zero_depth_positions', 0) / max(region_metrics.get('region_length', 1), 1),
                region_metrics.get('std_depth', 0) / max(region_metrics.get('mean_depth', 1), 1e-8),
                region_metrics.get('mean_mapq', 30) / 60.0,
                region_metrics.get('mean_baseq', 30) / 40.0,
                overall.get('mapping_rate', 0.95),
                region_metrics.get('gc_content', 0.5),
                region_metrics.get('min_depth', 0),
                region_metrics.get('max_depth', 0) / max(region_metrics.get('mean_depth', 1), 1e-8),
                region_metrics.get('low_depth_positions', 0) / max(region_metrics.get('region_length', 1), 1)
            ]
            
            return np.array(features, dtype=float)
            
        except Exception as e:
            self.logger.warning(f"Error extracting features for {region}: {e}")
            return None
    
    def _optimize_model_hyperparameters(self, X_train: np.ndarray, y_train: np.ndarray):
        """Optimize model hyperparameters using Optuna"""
        
        def objective(trial):
            n_estimators = trial.suggest_int('n_estimators', 50, 300)
            max_depth = trial.suggest_int('max_depth', 3, 20)
            min_samples_split = trial.suggest_int('min_samples_split', 2, 20)
            min_samples_leaf = trial.suggest_int('min_samples_leaf', 1, 10)
            max_features = trial.suggest_categorical('max_features', ['sqrt', 'log2', None])
            
            model = RandomForestClassifier(
                n_estimators=n_estimators,
                max_depth=max_depth,
                min_samples_split=min_samples_split,
                min_samples_leaf=min_samples_leaf,
                max_features=max_features,
                random_state=42,
                n_jobs=-1,
                class_weight='balanced'
            )
            
            # Cross-validation score
            cv_folds = min(5, len(np.unique(y_train)))
            scores = cross_val_score(model, X_train, y_train, cv=cv_folds, scoring='accuracy')
            return scores.mean()
        
        # Optimize hyperparameters
        study = optuna.create_study(direction='maximize')
        study.optimize(objective, n_trials=50, timeout=300)
        
        # Train final model with best parameters
        best_params = study.best_params
        final_model = RandomForestClassifier(**best_params, random_state=42, n_jobs=-1, class_weight='balanced')
        final_model.fit(X_train, y_train)
        
        # Calibrate probabilities for better confidence estimates
        calibrated_model = CalibratedClassifierCV(final_model, cv=3)
        calibrated_model.fit(X_train, y_train)
        
        self.logger.info(f"Best hyperparameters: {best_params}")
        return calibrated_model
    
    def _optimize_thresholds(self):
        """Optimize thresholds using validation data"""
        
        if len(self.validation_data) < 5:
            self.logger.warning("Insufficient validation data for threshold optimization")
            return
        
        def objective(trial):
            # Suggest threshold values
            homo_del_threshold = trial.suggest_float('homo_del_threshold', 0.05, 0.6)
            hetero_del_threshold = trial.suggest_float('hetero_del_threshold', 0.4, 0.9)
            dup_threshold = trial.suggest_float('dup_threshold', 1.1, 3.0)
            zscore_threshold = trial.suggest_float('zscore_threshold', 1.0, 4.0)
            confidence_threshold = trial.suggest_float('confidence_threshold', 0.5, 0.95)
            
            # Ensure logical ordering
            if hetero_del_threshold <= homo_del_threshold:
                return 0.0
            if dup_threshold <= hetero_del_threshold:
                return 0.0
            
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
            
            accuracy = self._evaluate_thresholds_on_validation(test_thresholds)
            return accuracy
        
        try:
            # Optimize thresholds
            study = optuna.create_study(direction='maximize')
            study.optimize(objective, n_trials=100, timeout=600)
            
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
            
        except Exception as e:
            self.logger.error(f"Error optimizing thresholds: {e}")
    
    def _evaluate_thresholds_on_validation(self, test_thresholds: Dict) -> float:
        """Evaluate threshold performance on validation data"""
        
        correct_predictions = 0
        total_predictions = 0
        
        # We need to simulate the pipeline's prediction logic here
        for val_sample in self.validation_data:
            # Find corresponding processed sample
            sample_id = val_sample['sample_id']
            matching_result = None
            
            for training_sample in self.training_data:
                if training_sample.get('sample_id') == sample_id:
                    matching_result = training_sample
                    break
            
            if not matching_result:
                continue
            
            for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
                if region in val_sample['mlpa_result']:
                    try:
                        prediction = self._make_threshold_prediction(matching_result, region, test_thresholds)
                        ground_truth = val_sample['mlpa_result'][region]
                        
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
            depth_metrics = sample.get('depth_metrics', {})
            if region not in depth_metrics:
                return 'NORMAL'
            
            normalized_depth = depth_metrics.get('normalized_depths', {}).get(region, {}).get('normalized_depth', 1.0)
            zscore = depth_metrics.get('normalized_depths', {}).get(region, {}).get('zscore', 0.0)
            
            # Apply thresholds
            depth_thresholds = thresholds['depth_ratio']
            stat_thresholds = thresholds['statistical']
            
            # Depth ratio prediction
            if normalized_depth < depth_thresholds['homo_del_threshold']:
                depth_call = 'HOMO_DEL'
            elif normalized_depth < depth_thresholds['hetero_del_threshold']:
                depth_call = 'HETERO_DEL'
            elif normalized_depth > depth_thresholds['dup_threshold']:
                depth_call = 'DUP'
            else:
                depth_call = 'NORMAL'
            
            # Statistical prediction
            zscore_threshold = stat_thresholds['zscore_threshold']
            if zscore < -zscore_threshold:
                if zscore < -3:
                    stat_call = 'HOMO_DEL'
                else:
                    stat_call = 'HETERO_DEL'
            elif zscore > zscore_threshold:
                stat_call = 'DUP'
            else:
                stat_call = 'NORMAL'
            
            # Simple consensus
            if depth_call == stat_call:
                return depth_call
            else:
                # For SMN genes, prefer depth ratio method
                return depth_call
                
        except Exception as e:
            return 'NORMAL'
    
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
        
        current_accuracy = 0.0
        if self.validation_data:
            current_accuracy = self._evaluate_thresholds_on_validation(self.thresholds)
        
        performance_record = {
            'timestamp': datetime.now().isoformat(),
            'samples_processed': self.thresholds['samples_processed'],
            'accuracy': current_accuracy,
            'thresholds': self.thresholds.copy(),
            'validation_samples': len(self.validation_data),
            'training_samples': len(self.training_data)
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
        
        self.logger.info(f"Current accuracy: {current_accuracy:.3f}")
        
        # Update thresholds performance metrics
        self.thresholds['performance_metrics'] = {
            'current_accuracy': current_accuracy,
            'last_evaluation': datetime.now().isoformat()
        }
    
    def get_performance_summary(self) -> Dict:
        """Get comprehensive performance summary"""
        
        if not self.performance_history:
            return {'status': 'no_data'}
        
        recent_performance = self.performance_history[-10:]
        accuracies = [p['accuracy'] for p in recent_performance]
        
        # Calculate trends
        accuracy_trend = 'stable'
        if len(accuracies) > 1:
            recent_avg = np.mean(accuracies[-3:]) if len(accuracies) >= 3 else accuracies[-1]
            older_avg = np.mean(accuracies[:3]) if len(accuracies) >= 6 else accuracies[0]
            
            if recent_avg > older_avg + 0.02:
                accuracy_trend = 'improving'
            elif recent_avg < older_avg - 0.02:
                accuracy_trend = 'declining'
        
        return {
            'current_accuracy': accuracies[-1] if accuracies else 0.0,
            'mean_accuracy': np.mean(accuracies),
            'std_accuracy': np.std(accuracies),
            'min_accuracy': np.min(accuracies),
            'max_accuracy': np.max(accuracies),
            'samples_processed': self.thresholds['samples_processed'],
            'last_updated': self.thresholds.get('last_updated'),
            'target_accuracy': self.config.get('ml_optimization', {}).get('target_accuracy', 0.95),
            'accuracy_trend': accuracy_trend,
            'validation_samples': len(self.validation_data),
            'training_samples': len(self.training_data),
            'model_available': 'classifier' in self.trained_models,
            'target_reached': accuracies[-1] >= 0.95 if accuracies else False
        }
    
    def should_retrain(self) -> bool:
        """Determine if models should be retrained"""
        
        retrain_interval = self.config.get('ml_optimization', {}).get('retrain_interval', 50)
        samples_processed = self.thresholds['samples_processed']
        
        # Retrain at intervals
        if samples_processed > 0 and samples_processed % retrain_interval == 0:
            return True
        
        # Retrain if accuracy is declining
        if len(self.performance_history) >= 5:
            recent_accuracies = [p['accuracy'] for p in self.performance_history[-5:]]
            if len(recent_accuracies) >= 3:
                recent_avg = np.mean(recent_accuracies[-2:])
                older_avg = np.mean(recent_accuracies[:2])
                
                if recent_avg < older_avg - 0.05:  # 5% drop
                    self.logger.info("Accuracy declining, triggering retrain")
                    return True
        
        # Retrain if we have significantly more data
        if len(self.training_data) > self.thresholds.get('last_training_size', 0) + 100:
            return True
        
        return False
    
    def export_model_performance(self) -> str:
        """Export detailed model performance report"""
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = self.model_dir / f'performance_report_{timestamp}.json'
        
        # Compile comprehensive report
        report = {
            'metadata': {
                'generated': datetime.now().isoformat(),
                'pipeline_version': self.config.get('pipeline_version', '1.0.0'),
                'samples_processed': self.thresholds['samples_processed']
            },
            'performance_summary': self.get_performance_summary(),
            'current_thresholds': self.thresholds,
            'performance_history': self.performance_history,
            'model_info': {
                'model_available': 'classifier' in self.trained_models,
                'scaler_available': 'scaler' in self.trained_models,
                'validation_samples': len(self.validation_data),
                'training_samples': len(self.training_data)
            }
        }
        
        # Add feature importance if model is available
        if 'classifier' in self.trained_models:
            try:
                model = self.trained_models['classifier']
                if hasattr(model, 'feature_importances_'):
                    feature_names = [
                        'normalized_depth', 'zscore', 'coverage_uniformity', 'zero_fraction',
                        'cv', 'mapq_norm', 'baseq_norm', 'mapping_rate', 'gc_content',
                        'min_depth', 'max_mean_ratio', 'low_depth_fraction'
                    ]
                    
                    importance_dict = dict(zip(feature_names, model.feature_importances_))
                    report['feature_importance'] = importance_dict
            except Exception as e:
                self.logger.warning(f"Could not extract feature importance: {e}")
        
        try:
            with open(report_file, 'w') as f:
                json.dump(report, f, indent=2)
            
            self.logger.info(f"Performance report exported: {report_file}")
            return str(report_file)
            
        except Exception as e:
            self.logger.error(f"Error exporting performance report: {e}")
            return None
    
    def reset_learning(self):
        """Reset ML learning (useful for testing or restarting)"""
        
        self.logger.info("Resetting ML learning state")
        
        # Clear training data
        self.training_data = []
        
        # Reset thresholds to defaults
        self.thresholds = self.config.get('cnv_thresholds', {})
        self.thresholds.update({
            'last_updated': datetime.now().isoformat(),
            'samples_processed': 0,
            'performance_metrics': {}
        })
        
        # Clear models
        self.trained_models = {}
        
        # Clear performance history
        self.performance_history = []
        
        # Save reset state
        self._save_thresholds()
        
        # Remove model files
        for model_file in [self.model_file, self.scaler_file, self.performance_file]:
            if model_file.exists():
                model_file.unlink()
        
        self.logger.info("ML learning state reset successfully")
    
    def get_learning_status(self) -> Dict:
        """Get current learning status"""
        
        performance = self.get_performance_summary()
        
        # Determine learning phase
        samples_processed = self.thresholds['samples_processed']
        if samples_processed < 20:
            phase = 'initialization'
        elif samples_processed < 100:
            phase = 'early_learning'
        elif samples_processed < 500:
            phase = 'active_learning'
        elif samples_processed < 1000:
            phase = 'refinement'
        else:
            phase = 'mature'
        
        # Determine if target is reached
        target_reached = performance.get('current_accuracy', 0) >= 0.95
        
        return {
            'phase': phase,
            'samples_processed': samples_processed,
            'target_reached': target_reached,
            'accuracy': performance.get('current_accuracy', 0),
            'trend': performance.get('accuracy_trend', 'unknown'),
            'next_milestone': self._get_next_milestone(samples_processed),
            'recommendations': self._get_learning_recommendations(performance, samples_processed)
        }
    
    def _get_next_milestone(self, samples_processed: int) -> Dict:
        """Get next learning milestone"""
        
        milestones = [20, 50, 100, 200, 500, 1000]
        
        for milestone in milestones:
            if samples_processed < milestone:
                return {
                    'samples': milestone,
                    'description': f"Reach {milestone} samples for improved accuracy"
                }
        
        return {
            'samples': samples_processed + 100,
            'description': "Continue processing for sustained high accuracy"
        }
    
    def _get_learning_recommendations(self, performance: Dict, samples_processed: int) -> List[str]:
        """Get learning recommendations based on current status"""
        
        recommendations = []
        current_accuracy = performance.get('current_accuracy', 0)
        trend = performance.get('accuracy_trend', 'stable')
        
        # Sample count recommendations
        if samples_processed < 20:
            recommendations.append("Process more samples to enable ML training (need 20+ samples)")
        elif samples_processed < 100:
            recommendations.append("Continue processing samples for better threshold optimization")
        
        # Accuracy recommendations
        if current_accuracy < 0.8:
            recommendations.append("Consider adding more MLPA validation data to improve accuracy")
            recommendations.append("Review failed samples for patterns in the reports")
        elif current_accuracy < 0.9:
            recommendations.append("Good progress! Continue processing samples to reach 95% target")
        elif current_accuracy >= 0.95:
            recommendations.append("Target accuracy reached! System is performing optimally")
        
        # Trend recommendations
        if trend == 'declining':
            recommendations.append("Performance declining - consider model retraining")
            recommendations.append("Check for systematic changes in sequencing protocol")
        elif trend == 'improving':
            recommendations.append("Performance improving - continue current approach")
        
        # Data recommendations
        validation_samples = performance.get('validation_samples', 0)
        if validation_samples < 20:
            recommendations.append("Add more MLPA validation samples for better optimization")
        
        # Model recommendations
        if not performance.get('model_available', False) and samples_processed > 50:
            recommendations.append("Sufficient data available - trigger model training")
        
        return recommendations
