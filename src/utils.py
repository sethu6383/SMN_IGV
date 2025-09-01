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
