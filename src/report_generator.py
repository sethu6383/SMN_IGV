#!/usr/bin/env python3
"""
Report Generator for SMN Pipeline
Generates HTML and TSV reports with clinical interpretations
"""

import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List
import logging
from datetime import datetime
import base64
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo

class ReportGenerator:
    """Generates comprehensive reports for SMN CNV analysis"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        self.reports_dir = Path(config.get('reports', 'reports'))
        self.reports_dir.mkdir(parents=True, exist_ok=True)
        
        # Clinical interpretation guidelines
        self.clinical_guidelines = self._load_clinical_guidelines()
    
    def _load_clinical_guidelines(self) -> Dict:
        """Load clinical interpretation guidelines"""
        return {
            'sma_risk_interpretation': {
                'HIGH_RISK': {
                    'description': 'High risk for SMA - Both SMN1 exons deleted',
                    'recommendation': 'Genetic counseling recommended. Consider confirmatory testing.',
                    'clinical_significance': 'Pathogenic',
                    'inheritance': 'Autosomal recessive'
                },
                'CARRIER': {
                    'description': 'SMA carrier - Heterozygous SMN1 deletion',
                    'recommendation': 'Genetic counseling for family planning. Partner testing recommended.',
                    'clinical_significance': 'Carrier status',
                    'inheritance': 'Autosomal recessive'
                },
                'LOW_RISK': {
                    'description': 'Low risk for SMA - Normal SMN1 copy number',
                    'recommendation': 'No additional testing typically required.',
                    'clinical_significance': 'Benign',
                    'inheritance': 'Not applicable'
                },
                'UNCERTAIN': {
                    'description': 'Uncertain significance - Atypical copy number pattern',
                    'recommendation': 'Consider orthogonal testing methods (MLPA, qPCR). Clinical correlation advised.',
                    'clinical_significance': 'Variant of uncertain significance',
                    'inheritance': 'Unknown'
                }
            },
            'copy_number_ranges': {
                'HOMO_DEL': '0 copies',
                'HETERO_DEL': '1 copy',
                'NORMAL': '2 copies',
                'DUP': '3+ copies'
            }
        }
    
    def generate_sample_report(self, sample_result: Dict) -> str:
        """Generate comprehensive report for a single sample"""
        
        sample_id = sample_result['sample_id']
        self.logger.info(f"Generating report for sample: {sample_id}")
        
        # Create sample-specific directory
        sample_dir = self.reports_dir / sample_id
        sample_dir.mkdir(exist_ok=True)
        
        # Generate plots
        plots = self._create_sample_plots(sample_result)
        
        # Create HTML report
        html_report = self._create_html_report(sample_result, plots)
        html_file = sample_dir / f"{sample_id}_report.html"
        
        with open(html_file, 'w') as f:
            f.write(html_report)
        
        # Generate TSV summary
        tsv_file = self._create_sample_tsv(sample_result, sample_dir)
        
        self.logger.info(f"Reports generated: {html_file}, {tsv_file}")
        return str(html_file)
    
    def generate_batch_report(self, batch_results: List[Dict]) -> str:
        """Generate consolidated batch report"""
        
        self.logger.info(f"Generating batch report for {len(batch_results)} samples")
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        batch_dir = self.reports_dir / f"batch_{timestamp}"
        batch_dir.mkdir(exist_ok=True)
        
        # Create batch plots
        batch_plots = self._create_batch_plots(batch_results)
        
        # Generate HTML report
        html_report = self._create_batch_html_report(batch_results, batch_plots)
        html_file = batch_dir / f"batch_report_{timestamp}.html"
        
        with open(html_file, 'w') as f:
            f.write(html_report)
        
        # Generate consolidated TSV
        tsv_file = self._create_batch_tsv(batch_results, batch_dir, timestamp)
        
        self.logger.info(f"Batch reports generated: {html_file}, {tsv_file}")
        return str(html_file)
    
    def _create_sample_plots(self, sample_result: Dict) -> Dict:
        """Create plots for single sample"""
        
        plots = {}
        
        if 'depth_metrics' not in sample_result:
            return plots
        
        depth_metrics = sample_result['depth_metrics']
        
        # 1. Depth comparison plot
        plots['depth_comparison'] = self._create_depth_comparison_plot(depth_metrics)
        
        # 2. Copy number plot
        plots['copy_number'] = self._create_copy_number_plot(sample_result)
        
        # 3. Quality metrics plot
        plots['quality_metrics'] = self._create_quality_metrics_plot(depth_metrics)
        
        # 4. CNV confidence plot
        plots['cnv_confidence'] = self._create_cnv_confidence_plot(sample_result)
        
        return plots
    
    def _create_depth_comparison_plot(self, depth_metrics: Dict) -> str:
        """Create depth comparison plot"""
        
        regions = ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']
        depths = [depth_metrics[region]['mean_depth'] for region in regions]
        normalized_depths = [depth_metrics['normalized_depths'][region]['normalized_depth'] for region in regions]
        
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=('Raw Depth', 'Normalized Depth'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        # Raw depths
        fig.add_trace(
            go.Bar(x=regions, y=depths, name='Raw Depth', marker_color='lightblue'),
            row=1, col=1
        )
        
        # Normalized depths
        colors = ['red' if d < 0.7 else 'green' if d < 1.3 else 'orange' for d in normalized_depths]
        fig.add_trace(
            go.Bar(x=regions, y=normalized_depths, name='Normalized Depth', marker_color=colors),
            row=1, col=2
        )
        
        # Add reference lines
        fig.add_hline(y=1.0, line_dash="dash", line_color="black", row=1, col=2)
        fig.add_hline(y=0.5, line_dash="dash", line_color="red", row=1, col=2)
        fig.add_hline(y=1.5, line_dash="dash", line_color="orange", row=1, col=2)
        
        fig.update_layout(
            title="SMN Gene Coverage Analysis",
            showlegend=False,
            height=400
        )
        
        return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    
    def _create_copy_number_plot(self, sample_result: Dict) -> str:
        """Create copy number visualization"""
        
        if 'cnv_calls' not in sample_result:
            return ""
        
        cnv_calls = sample_result['cnv_calls']
        regions = ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']
        
        # Map calls to copy numbers
        call_to_copies = {
            'HOMO_DEL': 0,
            'HETERO_DEL': 1,
            'NORMAL': 2,
            'DUP': 3
        }
        
        copy_numbers = []
        confidences = []
        colors = []
        
        for region in regions:
            call = cnv_calls[region]['call']
            confidence = cnv_calls[region]['confidence']
            
            copy_numbers.append(call_to_copies.get(call, 2))
            confidences.append(confidence)
            
            # Color based on call
            if call == 'HOMO_DEL':
                colors.append('red')
            elif call == 'HETERO_DEL':
                colors.append('orange')
            elif call == 'DUP':
                colors.append('purple')
            else:
                colors.append('green')
        
        fig = go.Figure()
        
        # Add bars with confidence as opacity
        for i, (region, copies, conf, color) in enumerate(zip(regions, copy_numbers, confidences, colors)):
            fig.add_trace(go.Bar(
                x=[region],
                y=[copies],
                marker_color=color,
                marker_opacity=conf,
                name=f"{region} (conf: {conf:.2f})",
                text=f"{copies} copies",
                textposition="middle center"
            ))
        
        fig.update_layout(
            title="Predicted Copy Numbers",
            yaxis_title="Copy Number",
            xaxis_title="SMN Regions",
            showlegend=False,
            height=400
        )
        
        return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    
    def _create_quality_metrics_plot(self, depth_metrics: Dict) -> str:
        """Create quality metrics visualization"""
        
        regions = ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']
        
        metrics_data = []
        for region in regions:
            metrics = depth_metrics[region]
            metrics_data.append({
                'Region': region,
                'Coverage Uniformity': metrics['coverage_uniformity'],
                'Mean MAPQ': metrics['mean_mapq'] / 60.0,  # Normalize to 0-1
                'Zero Depth Fraction': metrics['zero_depth_positions'] / metrics['region_length'],
                'GC Content': metrics['gc_content']
            })
        
        df = pd.DataFrame(metrics_data)
        
        # Create radar chart
        fig = go.Figure()
        
        for i, row in df.iterrows():
            fig.add_trace(go.Scatterpolar(
                r=[row['Coverage Uniformity'], row['Mean MAPQ'], 
                   1-row['Zero Depth Fraction'], row['GC Content']],
                theta=['Coverage Uniformity', 'Mapping Quality', 
                       'Coverage Completeness', 'GC Content'],
                fill='toself',
                name=row['Region']
            ))
        
        fig.update_layout(
            polar=dict(
                radialaxis=dict(visible=True, range=[0, 1])
            ),
            title="Quality Metrics by Region",
            height=500
        )
        
        return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    
    def _create_cnv_confidence_plot(self, sample_result: Dict) -> str:
        """Create CNV confidence visualization"""
        
        if 'cnv_calls' not in sample_result:
            return ""
        
        cnv_calls = sample_result['cnv_calls']
        regions = ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']
        
        # Get method-wise confidences
        methods = ['depth_ratio', 'statistical', 'hmm_like', 'ensemble']
        method_data = []
        
        for region in regions:
            if 'methods' in cnv_calls[region]:
                region_methods = cnv_calls[region]['methods']
                for method in methods:
                    if method in region_methods and 'confidence' in region_methods[method]:
                        method_data.append({
                            'Region': region,
                            'Method': method,
                            'Confidence': region_methods[method]['confidence'],
                            'Call': region_methods[method]['call']
                        })
        
        if not method_data:
            return ""
        
        df = pd.DataFrame(method_data)
        
        fig = px.bar(df, x='Region', y='Confidence', color='Method',
                     title='CNV Calling Confidence by Method',
                     barmode='group')
        
        fig.add_hline(y=0.8, line_dash="dash", line_color="red",
                      annotation_text="High Confidence Threshold")
        
        fig.update_layout(height=400)
        
        return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    
    def _create_batch_plots(self, batch_results: List[Dict]) -> Dict:
        """Create batch-level plots"""
        
        plots = {}
        
        # Filter successful results
        valid_results = [r for r in batch_results if 'error' not in r]
        
        if not valid_results:
            return plots
        
        # 1. Batch overview plot
        plots['batch_overview'] = self._create_batch_overview_plot(valid_results)
        
        # 2. SMA risk distribution
        plots['sma_distribution'] = self._create_sma_distribution_plot(valid_results)
        
        # 3. Quality distribution
        plots['quality_distribution'] = self._create_quality_distribution_plot(valid_results)
        
        # 4. Depth distribution
        plots['depth_distribution'] = self._create_depth_distribution_plot(valid_results)
        
        return plots
    
    def _create_batch_overview_plot(self, batch_results: List[Dict]) -> str:
        """Create batch overview plot"""
        
        # Collect CNV calls
        call_counts = {'NORMAL': 0, 'HETERO_DEL': 0, 'HOMO_DEL': 0, 'DUP': 0}
        
        for result in batch_results:
            cnv_calls = result.get('cnv_calls', {})
            for region in ['smn1_exon7', 'smn1_exon8']:  # Focus on SMN1
                if region in cnv_calls:
                    call = cnv_calls[region]['call']
                    call_counts[call] = call_counts.get(call, 0) + 1
        
        fig = go.Figure(data=[go.Pie(
            labels=list(call_counts.keys()),
            values=list(call_counts.values()),
            title="SMN1 CNV Distribution"
        )])
        
        fig.update_layout(title="Batch CNV Call Distribution", height=400)
        
        return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    
    def _create_sma_distribution_plot(self, batch_results: List[Dict]) -> str:
        """Create SMA risk distribution plot"""
        
        sma_risks = []
        for result in batch_results:
            sma_risk = result.get('cnv_calls', {}).get('sma_risk', 'UNCERTAIN')
            sma_risks.append(sma_risk)
        
        risk_counts = pd.Series(sma_risks).value_counts()
        
        colors = {
            'LOW_RISK': 'green',
            'CARRIER': 'orange', 
            'HIGH_RISK': 'red',
            'UNCERTAIN': 'gray'
        }
        
        fig = go.Figure(data=[go.Bar(
            x=risk_counts.index,
            y=risk_counts.values,
            marker_color=[colors.get(risk, 'blue') for risk in risk_counts.index]
        )])
        
        fig.update_layout(
            title="SMA Risk Distribution",
            xaxis_title="Risk Category",
            yaxis_title="Number of Samples",
            height=400
        )
        
        return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    
    def _create_quality_distribution_plot(self, batch_results: List[Dict]) -> str:
        """Create quality metrics distribution"""
        
        quality_data = []
        for result in batch_results:
            sample_id = result['sample_id']
            overall = result.get('depth_metrics', {}).get('overall', {})
            
            quality_data.append({
                'Sample': sample_id,
                'Mapping Rate': overall.get('mapping_rate', 0),
                'Mean MAPQ': overall.get('mean_mapq', 0)
            })
        
        df = pd.DataFrame(quality_data)
        
        fig = make_subplots(rows=1, cols=2, subplot_titles=('Mapping Rate', 'Mean MAPQ'))
        
        fig.add_trace(go.Histogram(x=df['Mapping Rate'], name='Mapping Rate'), row=1, col=1)
        fig.add_trace(go.Histogram(x=df['Mean MAPQ'], name='Mean MAPQ'), row=1, col=2)
        
        fig.update_layout(title="Quality Metrics Distribution", height=400, showlegend=False)
        
        return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    
    def _create_depth_distribution_plot(self, batch_results: List[Dict]) -> str:
        """Create depth distribution plot"""
        
        depth_data = []
        for result in batch_results:
            sample_id = result['sample_id']
            normalized_depths = result.get('depth_metrics', {}).get('normalized_depths', {})
            
            for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
                if region in normalized_depths:
                    depth_data.append({
                        'Sample': sample_id,
                        'Region': region,
                        'Normalized Depth': normalized_depths[region]['normalized_depth']
                    })
        
        df = pd.DataFrame(depth_data)
        
        fig = px.box(df, x='Region', y='Normalized Depth', 
                     title='Normalized Depth Distribution by Region')
        
        # Add reference lines
        fig.add_hline(y=1.0, line_dash="dash", line_color="black")
        fig.add_hline(y=0.5, line_dash="dash", line_color="red")
        fig.add_hline(y=1.5, line_dash="dash", line_color="orange")
        
        fig.update_layout(height=400)
        
        return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    
    def _create_html_report(self, sample_result: Dict, plots: Dict) -> str:
        """Create HTML report for single sample"""
        
        sample_id = sample_result['sample_id']
        
        # Get clinical interpretation
        clinical_interp = self._get_clinical_interpretation(sample_result)
        
        # Encode IGV snapshots
        encoded_snapshots = self._encode_snapshots(sample_result.get('igv_snapshots', {}))
        
        html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SMN CNV Report - {sample_id}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background-color: white; padding: 20px; border-radius: 10px; }}
        .header {{ background-color: #2c3e50; color: white; padding: 20px; border-radius: 5px; margin-bottom: 20px; }}
        .section {{ margin-bottom: 30px; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
        .cnv-call {{ padding: 10px; margin: 5px; border-radius: 5px; display: inline-block; font-weight: bold; }}
        .normal {{ background-color: #d4edda; color: #155724; }}
        .hetero-del {{ background-color: #fff3cd; color: #856404; }}
        .homo-del {{ background-color: #f8d7da; color: #721c24; }}
        .dup {{ background-color: #d1ecf1; color: #0c5460; }}
        .risk-high {{ background-color: #f8d7da; color: #721c24; padding: 15px; border-radius: 5px; }}
        .risk-carrier {{ background-color: #fff3cd; color: #856404; padding: 15px; border-radius: 5px; }}
        .risk-low {{ background-color: #d4edda; color: #155724; padding: 15px; border-radius: 5px; }}
        .snapshot {{ margin: 10px; text-align: center; }}
        .snapshot img {{ max-width: 100%; border: 1px solid #ddd; border-radius: 5px; }}
        .metrics-table {{ width: 100%; border-collapse: collapse; }}
        .metrics-table th, .metrics-table td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
        .metrics-table th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>SMN CNV Analysis Report</h1>
            <h2>Sample: {sample_id}</h2>
            <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="section">
            <h3>Clinical Summary</h3>
            {clinical_interp['summary_html']}
        </div>
        
        <div class="section">
            <h3>CNV Calls</h3>
            {self._create_cnv_calls_table(sample_result)}
        </div>
        
        <div class="section">
            <h3>Depth Analysis</h3>
            {plots.get('depth_comparison', '')}
        </div>
        
        <div class="section">
            <h3>Copy Number Predictions</h3>
            {plots.get('copy_number', '')}
        </div>
        
        <div class="section">
            <h3>Quality Metrics</h3>
            {plots.get('quality_metrics', '')}
            {self._create_quality_metrics_table(sample_result)}
        </div>
        
        <div class="section">
            <h3>Method Confidence</h3>
            {plots.get('cnv_confidence', '')}
        </div>
        
        <div class="section">
            <h3>IGV Snapshots</h3>
            {self._create_snapshots_html(encoded_snapshots)}
        </div>
        
        <div class="section">
            <h3>Technical Details</h3>
            {self._create_technical_details_table(sample_result)}
        </div>
    </div>
</body>
</html>
"""
        
        return html_template
    
    def _get_clinical_interpretation(self, sample_result: Dict) -> Dict:
        """Generate clinical interpretation"""
        
        sma_risk = sample_result.get('cnv_calls', {}).get('sma_risk', 'UNCERTAIN')
        guidelines = self.clinical_guidelines['sma_risk_interpretation'][sma_risk]
        
        # Create HTML summary
        risk_class = f"risk-{sma_risk.lower().replace('_', '-')}"
        
        summary_html = f"""
        <div class="{risk_class}">
            <h4>SMA Risk Assessment: {sma_risk}</h4>
            <p><strong>Description:</strong> {guidelines['description']}</p>
            <p><strong>Clinical Significance:</strong> {guidelines['clinical_significance']}</p>
            <p><strong>Recommendation:</strong> {guidelines['recommendation']}</p>
            <p><strong>Inheritance Pattern:</strong> {guidelines['inheritance']}</p>
        </div>
        """
        
        return {
            'risk_level': sma_risk,
            'guidelines': guidelines,
            'summary_html': summary_html
        }
    
    def _create_cnv_calls_table(self, sample_result: Dict) -> str:
        """Create CNV calls table HTML"""
        
        if 'cnv_calls' not in sample_result:
            return "<p>No CNV calls available</p>"
        
        cnv_calls = sample_result['cnv_calls']
        regions = ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']
        
        table_html = '<table class="metrics-table"><thead><tr><th>Region</th><th>CNV Call</th><th>Confidence</th><th>Copy Number</th></tr></thead><tbody>'
        
        call_to_copies = {'HOMO_DEL': '0', 'HETERO_DEL': '1', 'NORMAL': '2', 'DUP': '3+'}
        
        for region in regions:
            if region in cnv_calls:
                call = cnv_calls[region]['call']
                confidence = cnv_calls[region]['confidence']
                copies = call_to_copies.get(call, 'Unknown')
                
                call_class = call.lower().replace('_', '-')
                
                table_html += f"""
                <tr>
                    <td>{region.upper()}</td>
                    <td><span class="cnv-call {call_class}">{call}</span></td>
                    <td>{confidence:.3f}</td>
                    <td>{copies}</td>
                </tr>
                """
        
        table_html += '</tbody></table>'
        return table_html
    
    def _create_quality_metrics_table(self, sample_result: Dict) -> str:
        """Create quality metrics table"""
        
        if 'depth_metrics' not in sample_result:
            return ""
        
        overall = sample_result['depth_metrics'].get('overall', {})
        
        table_html = f"""
        <table class="metrics-table">
            <thead><tr><th>Metric</th><th>Value</th><th>Status</th></tr></thead>
            <tbody>
                <tr><td>Mapping Rate</td><td>{overall.get('mapping_rate', 0):.3f}</td><td>{'Good' if overall.get('mapping_rate', 0) > 0.95 else 'Poor'}</td></tr>
                <tr><td>Mean MAPQ</td><td>{overall.get('mean_mapq', 0):.1f}</td><td>{'Good' if overall.get('mean_mapq', 0) > 30 else 'Poor'}</td></tr>
                <tr><td>GC Content</td><td>{overall.get('gc_content', 0):.3f}</td><td>Normal</td></tr>
            </tbody>
        </table>
        """
        
        return table_html
    
    def _encode_snapshots(self, snapshots: Dict) -> Dict:
        """Encode IGV snapshots as base64 for HTML embedding"""
        
        encoded = {}
        
        for region_name, snapshot_file in snapshots.items():
            if snapshot_file and os.path.exists(snapshot_file):
                try:
                    with open(snapshot_file, 'rb') as f:
                        encoded_image = base64.b64encode(f.read()).decode('utf-8')
                    encoded[region_name] = encoded_image
                except Exception as e:
                    self.logger.warning(f"Could not encode snapshot {snapshot_file}: {e}")
                    encoded[region_name] = None
            else:
                encoded[region_name] = None
        
        return encoded
    
    def _create_snapshots_html(self, encoded_snapshots: Dict) -> str:
        """Create HTML for IGV snapshots"""
        
        if not encoded_snapshots:
            return "<p>No IGV snapshots available</p>"
        
        html = ""
        for region_name, encoded_image in encoded_snapshots.items():
            if encoded_image:
                html += f"""
                <div class="snapshot">
                    <h4>{region_name.upper().replace('_', ' ')}</h4>
                    <img src="data:image/png;base64,{encoded_image}" alt="{region_name} snapshot">
                </div>
                """
            else:
                html += f"""
                <div class="snapshot">
                    <h4>{region_name.upper().replace('_', ' ')}</h4>
                    <p>Snapshot not available</p>
                </div>
                """
        
        return html
    
    def _create_technical_details_table(self, sample_result: Dict) -> str:
        """Create technical details table"""
        
        table_html = f"""
        <table class="metrics-table">
            <thead><tr><th>Parameter</th><th>Value</th></tr></thead>
            <tbody>
                <tr><td>BAM File</td><td>{sample_result.get('bam_file', 'N/A')}</td></tr>
                <tr><td>Processing Time</td><td>{sample_result.get('processing_time', 'N/A')}</td></tr>
                <tr><td>Pipeline Version</td><td>{sample_result.get('pipeline_version', 'N/A')}</td></tr>
                <tr><td>Reference Genome</td><td>hg38</td></tr>
                <tr><td>Analysis Type</td><td>Whole Exome Sequencing (WES)</td></tr>
            </tbody>
        </table>
        """
        
        return table_html
    
    def _create_batch_html_report(self, batch_results: List[Dict], plots: Dict) -> str:
        """Create HTML report for batch processing"""
        
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        total_samples = len(batch_results)
        successful_samples = len([r for r in batch_results if 'error' not in r])
        
        html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SMN CNV Batch Report</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        .container {{ max-width: 1400px; margin: 0 auto; background-color: white; padding: 20px; border-radius: 10px; }}
        .header {{ background-color: #2c3e50; color: white; padding: 20px; border-radius: 5px; margin-bottom: 20px; }}
        .section {{ margin-bottom: 30px; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; }}
        .stat-card {{ background-color: #f8f9fa; padding: 15px; border-radius: 5px; text-align: center; }}
        .stat-number {{ font-size: 2em; font-weight: bold; color: #2c3e50; }}
        .stat-label {{ color: #6c757d; }}
        .results-table {{ width: 100%; border-collapse: collapse; font-size: 0.9em; }}
        .results-table th, .results-table td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
        .results-table th {{ background-color: #f2f2f2; position: sticky; top: 0; }}
        .results-container {{ max-height: 600px; overflow-y: auto; }}
        .cnv-call {{ padding: 3px 8px; border-radius: 3px; font-weight: bold; }}
        .normal {{ background-color: #d4edda; color: #155724; }}
        .hetero-del {{ background-color: #fff3cd; color: #856404; }}
        .homo-del {{ background-color: #f8d7da; color: #721c24; }}
        .dup {{ background-color: #d1ecf1; color: #0c5460; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>SMN CNV Batch Analysis Report</h1>
            <p>Analysis completed: {timestamp}</p>
        </div>
        
        <div class="section">
            <h3>Batch Summary</h3>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-number">{total_samples}</div>
                    <div class="stat-label">Total Samples</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{successful_samples}</div>
                    <div class="stat-label">Successful</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{total_samples - successful_samples}</div>
                    <div class="stat-label">Failed</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{self._count_high_risk(batch_results)}</div>
                    <div class="stat-label">High Risk</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{self._count_carriers(batch_results)}</div>
                    <div class="stat-label">Carriers</div>
                </div>
            </div>
        </div>
        
        <div class="section">
            <h3>Distribution Analysis</h3>
            {plots.get('batch_overview', '')}
            {plots.get('sma_distribution', '')}
        </div>
        
        <div class="section">
            <h3>Quality Assessment</h3>
            {plots.get('quality_distribution', '')}
            {plots.get('depth_distribution', '')}
        </div>
        
        <div class="section">
            <h3>Detailed Results</h3>
            <div class="results-container">
                {self._create_batch_results_table(batch_results)}
            </div>
        </div>
        
        <div class="section">
            <h3>Clinical Recommendations</h3>
            {self._create_batch_clinical_recommendations(batch_results)}
        </div>
    </div>
</body>
</html>
"""
        
        return html_template
    
    def _count_high_risk(self, batch_results: List[Dict]) -> int:
        """Count high-risk samples"""
        return sum(1 for r in batch_results if r.get('cnv_calls', {}).get('sma_risk') == 'HIGH_RISK')
    
    def _count_carriers(self, batch_results: List[Dict]) -> int:
        """Count carrier samples"""
        return sum(1 for r in batch_results if r.get('cnv_calls', {}).get('sma_risk') == 'CARRIER')
    
    def _create_batch_results_table(self, batch_results: List[Dict]) -> str:
        """Create detailed results table for batch"""
        
        table_html = """
        <table class="results-table">
            <thead>
                <tr>
                    <th>Sample ID</th>
                    <th>SMA Risk</th>
                    <th>SMN1 Ex7</th>
                    <th>SMN1 Ex8</th>
                    <th>SMN2 Ex7</th>
                    <th>SMN2 Ex8</th>
                    <th>Confidence</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>
        """
        
        for result in batch_results:
            sample_id = result['sample_id']
            
            if 'error' in result:
                table_html += f"""
                <tr>
                    <td>{sample_id}</td>
                    <td colspan="7" style="color: red;">Error: {result['error']}</td>
                </tr>
                """
            else:
                cnv_calls = result.get('cnv_calls', {})
                sma_risk = cnv_calls.get('sma_risk', 'UNCERTAIN')
                
                # Get calls and confidences
                regions = ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']
                calls = []
                confidences = []
                
                for region in regions:
                    if region in cnv_calls:
                        call = cnv_calls[region]['call']
                        conf = cnv_calls[region]['confidence']
                        call_class = call.lower().replace('_', '-')
                        calls.append(f'<span class="cnv-call {call_class}">{call}</span>')
                        confidences.append(conf)
                    else:
                        calls.append('N/A')
                        confidences.append(0.0)
                
                avg_confidence = np.mean([c for c in confidences if c > 0])
                
                table_html += f"""
                <tr>
                    <td>{sample_id}</td>
                    <td><span class="cnv-call {sma_risk.lower().replace('_', '-')}">{sma_risk}</span></td>
                    <td>{calls[0]}</td>
                    <td>{calls[1]}</td>
                    <td>{calls[2]}</td>
                    <td>{calls[3]}</td>
                    <td>{avg_confidence:.3f}</td>
                    <td>Success</td>
                </tr>
                """
        
        table_html += '</tbody></table>'
        return table_html
    
    def _create_batch_clinical_recommendations(self, batch_results: List[Dict]) -> str:
        """Create clinical recommendations for batch"""
        
        high_risk_samples = []
        carrier_samples = []
        uncertain_samples = []
        
        for result in batch_results:
            if 'error' not in result:
                sma_risk = result.get('cnv_calls', {}).get('sma_risk', 'UNCERTAIN')
                sample_id = result['sample_id']
                
                if sma_risk == 'HIGH_RISK':
                    high_risk_samples.append(sample_id)
                elif sma_risk == 'CARRIER':
                    carrier_samples.append(sample_id)
                elif sma_risk == 'UNCERTAIN':
                    uncertain_samples.append(sample_id)
        
        recommendations_html = ""
        
        if high_risk_samples:
            recommendations_html += f"""
            <div class="risk-high">
                <h4>High Risk Samples ({len(high_risk_samples)})</h4>
                <p><strong>Samples:</strong> {', '.join(high_risk_samples)}</p>
                <p><strong>Action Required:</strong> Immediate genetic counseling and confirmatory testing recommended.</p>
            </div>
            """
        
        if carrier_samples:
            recommendations_html += f"""
            <div class="risk-carrier">
                <h4>Carrier Samples ({len(carrier_samples)})</h4>
                <p><strong>Samples:</strong> {', '.join(carrier_samples)}</p>
                <p><strong>Action Required:</strong> Genetic counseling for family planning. Partner testing recommended.</p>
            </div>
            """
        
        if uncertain_samples:
            recommendations_html += f"""
            <div class="risk-uncertain">
                <h4>Uncertain Samples ({len(uncertain_samples)})</h4>
                <p><strong>Samples:</strong> {', '.join(uncertain_samples)}</p>
                <p><strong>Action Required:</strong> Consider orthogonal testing methods (MLPA, qPCR). Review IGV snapshots.</p>
            </div>
            """
        
        return recommendations_html or "<p>No specific recommendations for this batch.</p>"
    
    def _create_sample_tsv(self, sample_result: Dict, sample_dir: Path) -> str:
        """Create TSV summary for single sample"""
        
        sample_id = sample_result['sample_id']
        tsv_file = sample_dir / f"{sample_id}_summary.tsv"
        
        # Flatten data for TSV
        if 'cnv_calls' in sample_result:
            cnv_calls = sample_result['cnv_calls']
            
            tsv_data = [{
                'sample_id': sample_id,
                'sma_risk': cnv_calls.get('sma_risk', 'UNCERTAIN'),
                'smn1_exon7_call': cnv_calls.get('smn1_exon7', {}).get('call', 'N/A'),
                'smn1_exon7_confidence': cnv_calls.get('smn1_exon7', {}).get('confidence', 0),
                'smn1_exon8_call': cnv_calls.get('smn1_exon8', {}).get('call', 'N/A'),
                'smn1_exon8_confidence': cnv_calls.get('smn1_exon8', {}).get('confidence', 0),
                'smn2_exon7_call': cnv_calls.get('smn2_exon7', {}).get('call', 'N/A'),
                'smn2_exon7_confidence': cnv_calls.get('smn2_exon7', {}).get('confidence', 0),
                'smn2_exon8_call': cnv_calls.get('smn2_exon8', {}).get('call', 'N/A'),
                'smn2_exon8_confidence': cnv_calls.get('smn2_exon8', {}).get('confidence', 0),
                'processing_time': sample_result.get('processing_time', 'N/A')
            }]
        else:
            tsv_data = [{
                'sample_id': sample_id,
                'error': sample_result.get('error', 'Unknown error')
            }]
        
        df = pd.DataFrame(tsv_data)
        df.to_csv(tsv_file, sep='\t', index=False)
        
        return str(tsv_file)

class ClinicalReportGenerator:
    """Specialized clinical report generator with enhanced interpretations"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Clinical guidelines and references
        self.clinical_guidelines = self._load_enhanced_clinical_guidelines()
        
    def _load_enhanced_clinical_guidelines(self) -> Dict:
        """Load enhanced clinical interpretation guidelines"""
        return {
            'sma_types': {
                'type1': {
                    'description': 'Severe infantile SMA (Werdnig-Hoffmann disease)',
                    'onset': '0-6 months',
                    'motor_milestones': 'Never sits without support',
                    'life_expectancy': '<2 years without intervention',
                    'smn1_copies': 0,
                    'smn2_copies': '1-2'
                },
                'type2': {
                    'description': 'Intermediate SMA (Dubowitz disease)',
                    'onset': '6-18 months',
                    'motor_milestones': 'Sits without support, never walks',
                    'life_expectancy': 'Variable, often into adulthood',
                    'smn1_copies': 0,
                    'smn2_copies': '3'
                },
                'type3': {
                    'description': 'Mild SMA (Kugelberg-Welander disease)',
                    'onset': '>18 months',
                    'motor_milestones': 'Walks independently',
                    'life_expectancy': 'Normal or near-normal',
                    'smn1_copies': 0,
                    'smn2_copies': '4+'
                },
                'type4': {
                    'description': 'Adult-onset SMA',
                    'onset': '>30 years',
                    'motor_milestones': 'Normal development initially',
                    'life_expectancy': 'Normal',
                    'smn1_copies': 0,
                    'smn2_copies': '4+'
                }
            },
            'therapeutic_options': {
                'nusinersen': {
                    'name': 'Nusinersen (Spinraza)',
                    'mechanism': 'Antisense oligonucleotide',
                    'indication': 'All SMA types',
                    'administration': 'Intrathecal injection'
                },
                'onasemnogene': {
                    'name': 'Onasemnogene abeparvovec (Zolgensma)',
                    'mechanism': 'Gene replacement therapy',
                    'indication': 'SMA Type 1, <2 years old',
                    'administration': 'Single IV infusion'
                },
                'risdiplam': {
                    'name': 'Risdiplam (Evrysdi)',
                    'mechanism': 'SMN2 splicing modifier',
                    'indication': 'All SMA types, >2 months old',
                    'administration': 'Oral daily'
                }
            },
            'genetic_counseling_points': [
                'SMA is inherited in an autosomal recessive manner',
                'Carrier frequency is approximately 1 in 50-60 individuals',
                'Risk for SMA in offspring depends on partner carrier status',
                'Prenatal and preimplantation genetic testing available',
                'Family cascade testing recommended for at-risk relatives'
            ]
        }
    
    def generate_clinical_summary(self, cnv_results: Dict, patient_info: Dict = None) -> Dict:
        """Generate comprehensive clinical summary"""
        
        sma_risk = cnv_results.get('sma_risk', 'UNCERTAIN')
        
        # Determine likely SMA type based on SMN2 copies (if available)
        predicted_sma_type = self._predict_sma_type(cnv_results)
        
        # Generate clinical recommendations
        clinical_recommendations = self._generate_clinical_recommendations(sma_risk, predicted_sma_type)
        
        # Generate follow-up actions
        follow_up_actions = self._generate_follow_up_actions(sma_risk, cnv_results)
        
        return {
            'sma_risk_level': sma_risk,
            'predicted_sma_type': predicted_sma_type,
            'clinical_recommendations': clinical_recommendations,
            'follow_up_actions': follow_up_actions,
            'genetic_counseling_indicated': sma_risk in ['HIGH_RISK', 'CARRIER'],
            'urgent_referral_needed': sma_risk == 'HIGH_RISK',
            'therapeutic_options': self._get_therapeutic_options(sma_risk, predicted_sma_type),
            'family_testing_recommended': sma_risk in ['HIGH_RISK', 'CARRIER']
        }
    
    def _predict_sma_type(self, cnv_results: Dict) -> Optional[str]:
        """Predict likely SMA type based on copy number results"""
        
        smn1_ex7 = cnv_results.get('smn1_exon7', {}).get('call', 'NORMAL')
        smn1_ex8 = cnv_results.get('smn1_exon8', {}).get('call', 'NORMAL')
        smn2_ex7 = cnv_results.get('smn2_exon7', {}).get('call', 'NORMAL')
        smn2_ex8 = cnv_results.get('smn2_exon8', {}).get('call', 'NORMAL')
        
        # Only predict if SMN1 is deleted
        if not (smn1_ex7 in ['HOMO_DEL', 'HETERO_DEL'] or smn1_ex8 in ['HOMO_DEL', 'HETERO_DEL']):
            return None
        
        # Estimate SMN2 copies (simplified - would need more sophisticated analysis in practice)
        smn2_copies = 0
        if smn2_ex7 == 'NORMAL':
            smn2_copies += 2
        elif smn2_ex7 == 'DUP':
            smn2_copies += 3
        elif smn2_ex7 == 'HETERO_DEL':
            smn2_copies += 1
        
        if smn2_ex8 == 'NORMAL':
            smn2_copies += 2
        elif smn2_ex8 == 'DUP':
            smn2_copies += 3
        elif smn2_ex8 == 'HETERO_DEL':
            smn2_copies += 1
        
        # Average SMN2 copies
        estimated_smn2_copies = smn2_copies / 2
        
        # Predict SMA type based on SMN2 copies
        if estimated_smn2_copies <= 2:
            return 'type1'
        elif estimated_smn2_copies == 3:
            return 'type2'
        elif estimated_smn2_copies >= 4:
            return 'type3'
        
        return 'uncertain'
    
    def _generate_clinical_recommendations(self, sma_risk: str, predicted_sma_type: Optional[str]) -> List[str]:
        """Generate clinical recommendations based on results"""
        
        recommendations = []
        
        if sma_risk == 'HIGH_RISK':
            recommendations.extend([
                'URGENT: Refer to neurology/genetics specialist immediately',
                'Consider confirmatory testing with MLPA or qPCR',
                'Genetic counseling strongly recommended',
                'Discuss therapeutic options if diagnosis confirmed',
                'Family cascade testing recommended'
            ])
            
            if predicted_sma_type:
                sma_info = self.clinical_guidelines['sma_types'].get(predicted_sma_type, {})
                if predicted_sma_type == 'type1':
                    recommendations.extend([
                        'Early intervention critical - consider immediate therapy',
                        'Respiratory and nutritional support may be needed',
                        'Discuss prognosis and quality of life considerations'
                    ])
                elif predicted_sma_type in ['type2', 'type3']:
                    recommendations.extend([
                        'Discuss long-term management and prognosis',
                        'Physical therapy and supportive care planning',
                        'Monitor for disease progression'
                    ])
        
        elif sma_risk == 'CARRIER':
            recommendations.extend([
                'Genetic counseling recommended for family planning',
                'Partner testing recommended before conception',
                'Discuss reproductive options if partner is also carrier',
                'Family cascade testing may be appropriate',
                'No immediate medical intervention needed'
            ])
        
        elif sma_risk == 'LOW_RISK':
            recommendations.extend([
                'No further SMA-specific testing typically required',
                'Routine genetic counseling if family history present',
                'Results consistent with normal SMA risk'
            ])
        
        elif sma_risk == 'UNCERTAIN':
            recommendations.extend([
                'Consider orthogonal testing methods (MLPA, qPCR)',
                'Genetic counseling recommended for interpretation',
                'Clinical correlation advised',
                'May require repeat testing or additional analysis'
            ])
        
        return recommendations
    
    def _generate_follow_up_actions(self, sma_risk: str, cnv_results: Dict) -> List[Dict]:
        """Generate specific follow-up actions with timelines"""
        
        actions = []
        
        if sma_risk == 'HIGH_RISK':
            actions.extend([
                {
                    'action': 'Specialist referral',
                    'specialty': 'Neurology/Genetics',
                    'urgency': 'Urgent (within 1-2 weeks)',
                    'purpose': 'Confirm diagnosis and discuss management'
                },
                {
                    'action': 'Confirmatory testing',
                    'test': 'MLPA or quantitative PCR',
                    'urgency': 'Urgent (within 1 week)',
                    'purpose': 'Validate WES findings'
                },
                {
                    'action': 'Genetic counseling',
                    'urgency': 'Urgent (within 1-2 weeks)',
                    'purpose': 'Discuss implications and family planning'
                }
            ])
        
        elif sma_risk == 'CARRIER':
            actions.extend([
                {
                    'action': 'Genetic counseling',
                    'urgency': 'Routine (within 4-6 weeks)',
                    'purpose': 'Discuss carrier status and family planning'
                },
                {
                    'action': 'Partner testing',
                    'urgency': 'Before conception if planning pregnancy',
                    'purpose': 'Assess recurrence risk'
                },
                {
                    'action': 'Family screening',
                    'urgency': 'As appropriate',
                    'purpose': 'Identify other at-risk family members'
                }
            ])
        
        elif sma_risk == 'UNCERTAIN':
            actions.extend([
                {
                    'action': 'Repeat/confirmatory testing',
                    'test': 'MLPA or alternative method',
                    'urgency': 'Routine (within 2-4 weeks)',
                    'purpose': 'Clarify uncertain results'
                },
                {
                    'action': 'Genetic counseling',
                    'urgency': 'Routine (within 4 weeks)',
                    'purpose': 'Interpret results and discuss next steps'
                }
            ])
        
        # Check for quality issues that might affect interpretation
        quality_issues = self._assess_result_quality(cnv_results)
        if quality_issues:
            actions.append({
                'action': 'Technical review',
                'urgency': 'Routine',
                'purpose': f'Address quality concerns: {"; ".join(quality_issues)}'
            })
        
        return actions
    
    def _assess_result_quality(self, cnv_results: Dict) -> List[str]:
        """Assess quality issues that might affect clinical interpretation"""
        
        quality_issues = []
        
        for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
            if region in cnv_results:
                confidence = cnv_results[region].get('confidence', 0)
                
                if confidence < 0.7:
                    quality_issues.append(f"Low confidence for {region}")
                
                # Check for method disagreement if available
                methods = cnv_results[region].get('methods', {})
                if len(methods) > 1:
                    calls = [m.get('call') for m in methods.values() if 'call' in m]
                    if len(set(calls)) > 1:  # Disagreement between methods
                        quality_issues.append(f"Method disagreement for {region}")
        
        return quality_issues
    
    def _get_therapeutic_options(self, sma_risk: str, predicted_sma_type: Optional[str]) -> List[Dict]:
        """Get relevant therapeutic options"""
        
        if sma_risk != 'HIGH_RISK':
            return []
        
        therapeutic_options = []
        
        for therapy_id, therapy_info in self.clinical_guidelines['therapeutic_options'].items():
            # Add all therapies for high-risk patients with appropriate indications
            therapeutic_options.append({
                'name': therapy_info['name'],
                'mechanism': therapy_info['mechanism'],
                'indication': therapy_info['indication'],
                'administration': therapy_info['administration'],
                'note': 'Consult specialist for eligibility and timing'
            })
        
        return therapeutic_options

class QualityMetricsReporter:
    """Generates quality metrics reports for batches and individual samples"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def generate_batch_quality_report(self, batch_results: List[Dict]) -> Dict:
        """Generate comprehensive quality metrics for a batch"""
        
        valid_results = [r for r in batch_results if 'error' not in r]
        
        if not valid_results:
            return {'error': 'No valid results for quality analysis'}
        
        # Collect metrics
        depth_metrics = self._collect_depth_metrics(valid_results)
        confidence_metrics = self._collect_confidence_metrics(valid_results)
        coverage_metrics = self._collect_coverage_metrics(valid_results)
        quality_scores = self._calculate_quality_scores(valid_results)
        
        # Generate summary statistics
        quality_report = {
            'batch_summary': {
                'total_samples': len(batch_results),
                'valid_samples': len(valid_results),
                'failed_samples': len(batch_results) - len(valid_results),
                'success_rate': len(valid_results) / len(batch_results)
            },
            'depth_metrics': depth_metrics,
            'confidence_metrics': confidence_metrics,
            'coverage_metrics': coverage_metrics,
            'quality_scores': quality_scores,
            'recommendations': self._generate_quality_recommendations(valid_results)
        }
        
        return quality_report
    
    def _collect_depth_metrics(self, results: List[Dict]) -> Dict:
        """Collect depth metrics across samples"""
        
        depth_data = {region: [] for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']}
        normalized_depth_data = {region: [] for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']}
        
        for result in results:
            depth_metrics = result.get('depth_metrics', {})
            
            for region in depth_data.keys():
                if region in depth_metrics:
                    depth_data[region].append(depth_metrics[region].get('mean_depth', 0))
                    
                normalized_depths = depth_metrics.get('normalized_depths', {})
                if region in normalized_depths:
                    normalized_depth_data[region].append(normalized_depths[region].get('normalized_depth', 1.0))
        
        # Calculate statistics
        depth_stats = {}
        for region, depths in depth_data.items():
            if depths:
                depth_stats[region] = {
                    'mean': np.mean(depths),
                    'median': np.median(depths),
                    'std': np.std(depths),
                    'min': np.min(depths),
                    'max': np.max(depths),
                    'samples': len(depths)
                }
        
        normalized_depth_stats = {}
        for region, depths in normalized_depth_data.items():
            if depths:
                normalized_depth_stats[region] = {
                    'mean': np.mean(depths),
                    'median': np.median(depths),
                    'std': np.std(depths),
                    'cv': np.std(depths) / np.mean(depths) if np.mean(depths) > 0 else 0,
                    'samples': len(depths)
                }
        
        return {
            'raw_depths': depth_stats,
            'normalized_depths': normalized_depth_stats
        }
    
    def _collect_confidence_metrics(self, results: List[Dict]) -> Dict:
        """Collect confidence metrics across samples"""
        
        confidence_data = {region: [] for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']}
        
        for result in results:
            cnv_calls = result.get('cnv_calls', {})
            
            for region in confidence_data.keys():
                if region in cnv_calls:
                    confidence_data[region].append(cnv_calls[region].get('confidence', 0))
        
        # Calculate statistics
        confidence_stats = {}
        for region, confidences in confidence_data.items():
            if confidences:
                confidence_stats[region] = {
                    'mean': np.mean(confidences),
                    'median': np.median(confidences),
                    'min': np.min(confidences),
                    'max': np.max(confidences),
                    'high_confidence': sum(1 for c in confidences if c > 0.8) / len(confidences),
                    'low_confidence': sum(1 for c in confidences if c < 0.5) / len(confidences),
                    'samples': len(confidences)
                }
        
        return confidence_stats
    
    def _collect_coverage_metrics(self, results: List[Dict]) -> Dict:
        """Collect coverage quality metrics"""
        
        coverage_data = {
            'mapping_rates': [],
            'mean_mapqs': [],
            'coverage_uniformity': {region: [] for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']}
        }
        
        for result in results:
            depth_metrics = result.get('depth_metrics', {})
            overall = depth_metrics.get('overall', {})
            
            # Overall metrics
            coverage_data['mapping_rates'].append(overall.get('mapping_rate', 0))
            coverage_data['mean_mapqs'].append(overall.get('mean_mapq', 0))
            
            # Region-specific coverage uniformity
            for region in coverage_data['coverage_uniformity'].keys():
                if region in depth_metrics:
                    uniformity = depth_metrics[region].get('coverage_uniformity', 0)
                    coverage_data['coverage_uniformity'][region].append(uniformity)
        
        # Calculate statistics
        coverage_stats = {
            'mapping_rate': {
                'mean': np.mean(coverage_data['mapping_rates']),
                'min': np.min(coverage_data['mapping_rates']),
                'samples_below_95': sum(1 for r in coverage_data['mapping_rates'] if r < 0.95)
            },
            'mean_mapq': {
                'mean': np.mean(coverage_data['mean_mapqs']),
                'min': np.min(coverage_data['mean_mapqs']),
                'samples_below_30': sum(1 for m in coverage_data['mean_mapqs'] if m < 30)
            },
            'coverage_uniformity': {}
        }
        
        for region, uniformities in coverage_data['coverage_uniformity'].items():
            if uniformities:
                coverage_stats['coverage_uniformity'][region] = {
                    'mean': np.mean(uniformities),
                    'min': np.min(uniformities),
                    'samples_below_80': sum(1 for u in uniformities if u < 0.8) / len(uniformities)
                }
        
        return coverage_stats
    
    def _calculate_quality_scores(self, results: List[Dict]) -> Dict:
        """Calculate overall quality scores for samples"""
        
        quality_scores = []
        region_quality = {region: [] for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']}
        
        for result in results:
            sample_quality = self._calculate_sample_quality_score(result)
            quality_scores.append(sample_quality['overall'])
            
            for region, score in sample_quality['regions'].items():
                if region in region_quality:
                    region_quality[region].append(score)
        
        # Overall quality distribution
        quality_distribution = {
            'high_quality': sum(1 for q in quality_scores if q > 0.8) / len(quality_scores),
            'medium_quality': sum(1 for q in quality_scores if 0.6 <= q <= 0.8) / len(quality_scores),
            'low_quality': sum(1 for q in quality_scores if q < 0.6) / len(quality_scores),
            'mean_score': np.mean(quality_scores)
        }
        
        # Region quality scores
        region_quality_stats = {}
        for region, scores in region_quality.items():
            if scores:
                region_quality_stats[region] = {
                    'mean': np.mean(scores),
                    'min': np.min(scores),
                    'max': np.max(scores)
                }
        
        return {
            'distribution': quality_distribution,
            'region_quality': region_quality_stats
        }
    
    def _calculate_sample_quality_score(self, result: Dict) -> Dict:
        """Calculate quality score for individual sample"""
        
        depth_metrics = result.get('depth_metrics', {})
        cnv_calls = result.get('cnv_calls', {})
        overall = depth_metrics.get('overall', {})
        
        # Overall quality factors
        mapping_rate_score = min(1.0, overall.get('mapping_rate', 0) / 0.95)
        mapq_score = min(1.0, overall.get('mean_mapq', 0) / 40.0)
        
        region_scores = {}
        region_quality_factors = []
        
        for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
            if region in depth_metrics and region in cnv_calls:
                # Region-specific quality factors
                uniformity = depth_metrics[region].get('coverage_uniformity', 0)
                depth = depth_metrics[region].get('mean_depth', 0)
                confidence = cnv_calls[region].get('confidence', 0)
                
                depth_score = min(1.0, depth / 20.0)  # Target depth 20x
                uniformity_score = uniformity
                confidence_score = confidence
                
                region_score = np.mean([depth_score, uniformity_score, confidence_score, mapping_rate_score, mapq_score])
                region_scores[region] = region_score
                region_quality_factors.append(region_score)
        
        overall_score = np.mean(region_quality_factors) if region_quality_factors else 0
        
        return {
            'overall': overall_score,
            'regions': region_scores
        }
    
    def _generate_quality_recommendations(self, results: List[Dict]) -> List[str]:
        """Generate quality improvement recommendations"""
        
        recommendations = []
        
        # Analyze common issues
        low_mapping_samples = 0
        low_depth_samples = 0
        low_confidence_samples = 0
        
        for result in results:
            depth_metrics = result.get('depth_metrics', {})
            cnv_calls = result.get('cnv_calls', {})
            overall = depth_metrics.get('overall', {})
            
            if overall.get('mapping_rate', 1.0) < 0.95:
                low_mapping_samples += 1
            
            # Check for low depth in SMN regions
            smn_depths = []
            for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
                if region in depth_metrics:
                    smn_depths.append(depth_metrics[region].get('mean_depth', 0))
            
            if smn_depths and np.mean(smn_depths) < 15:
                low_depth_samples += 1
            
            # Check for low confidence calls
            confidences = []
            for region in ['smn1_exon7', 'smn1_exon8', 'smn2_exon7', 'smn2_exon8']:
                if region in cnv_calls:
                    confidences.append(cnv_calls[region].get('confidence', 0))
            
            if confidences and np.mean(confidences) < 0.7:
                low_confidence_samples += 1
        
        total_samples = len(results)
        
        # Generate recommendations
        if low_mapping_samples > total_samples * 0.1:
            recommendations.append(f"High mapping rate issues ({low_mapping_samples} samples) - check library prep and alignment parameters")
        
        if low_depth_samples > total_samples * 0.1:
            recommendations.append(f"Low depth in SMN regions ({low_depth_samples} samples) - consider targeted enrichment or increased sequencing depth")
        
        if low_confidence_samples > total_samples * 0.2:
            recommendations.append(f"Low confidence calls ({low_confidence_samples} samples) - review with orthogonal methods")
        
        if not recommendations:
            recommendations.append("Overall good quality metrics across the batch")
        
        return recommendations
    
    def _create_batch_tsv(self, batch_results: List[Dict], batch_dir: Path, timestamp: str) -> str:
        """Create consolidated TSV for batch"""
        
        tsv_file = batch_dir / f"batch_results_{timestamp}.tsv"
        
        # Collect all data
        all_data = []
        for result in batch_results:
            if 'error' not in result and 'cnv_calls' in result:
                cnv_calls = result['cnv_calls']
                depth_metrics = result.get('depth_metrics', {})
                
                row_data = {
                    'sample_id': result['sample_id'],
                    'sma_risk': cnv_calls.get('sma_risk', 'UNCERTAIN'),
                    'smn1_exon7_call': cnv_calls.get('smn1_exon7', {}).get('call', 'N/A'),
                    'smn1_exon7_confidence': cnv_calls.get('smn1_exon7', {}).get('confidence', 0),
                    'smn1_exon8_call': cnv_calls.get('smn1_exon8', {}).get('call', 'N/A'),
                    'smn1_exon8_confidence': cnv_calls.get('smn1_exon8', {}).get('confidence', 0),
                    'smn2_exon7_call': cnv_calls.get('smn2_exon7', {}).get('call', 'N/A'),
                    'smn2_exon7_confidence': cnv_calls.get('smn2_exon7', {}).get('confidence', 0),
                    'smn2_exon8_call': cnv_calls.get('smn2_exon8', {}).get('call', 'N/A'),
                    'smn2_exon8_confidence': cnv_calls.get('smn2_exon8', {}).get('confidence', 0)
                }
                
                # Add depth metrics if available
                if 'smn1_exon7' in depth_metrics:
                    row_data.update({
                        'smn1_exon7_depth': depth_metrics['smn1_exon7']['mean_depth'],
                        'smn1_exon8_depth': depth_metrics['smn1_exon8']['mean_depth'],
                        'smn2_exon7_depth': depth_metrics['smn2_exon7']['mean_depth'],
                        'smn2_exon8_depth': depth_metrics['smn2_exon8']['mean_depth'],
                        'mapping_rate': depth_metrics.get('overall', {}).get('mapping_rate', 0),
                        'mean_mapq': depth_metrics.get('overall', {}).get('mean_mapq', 0)
                    })
                
                all_data.append(row_data)
            else:
                all_data.append({
                    'sample_id': result['sample_id'],
                    'error': result.get('error', 'Processing failed')
                })
        
        df = pd.DataFrame(all_data)
        df.to_csv(tsv_file, sep='\t', index=False)
        
        return str(tsv_file)
