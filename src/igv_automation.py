#!/usr/bin/env python3
"""
IGV Automation Module for SMN Pipeline
Automates IGV snapshots for visual CNV assessment
"""

import os
import subprocess
import logging
import time
from typing import Dict, List
from pathlib import Path
import tempfile

class IGVAutomation:
    """Automates IGV for taking screenshots of SMN regions"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # IGV settings
        self.igv_path = config.get('igv_path', 'igv.sh')
        self.igv_memory = config.get('igv_memory', '4g')
        self.snapshot_dir = config.get('igv_snapshots', 'snapshots')
        
        # Create snapshot directory
        Path(self.snapshot_dir).mkdir(parents=True, exist_ok=True)
        
        # SMN regions for snapshots
        self.regions = {
            'smn1_exon7': 'chr5:70247700-70247800',
            'smn1_exon8': 'chr5:70248900-70249350',
            'smn2_exon7': 'chr5:69372280-69372380', 
            'smn2_exon8': 'chr5:69373490-69373920',
            'smn1_overview': 'chr5:70247000-70250000',
            'smn2_overview': 'chr5:69372000-69375000'
        }
    
    def generate_snapshots(self, bam_file: str, sample_id: str) -> Dict:
        """Generate IGV snapshots for all SMN regions"""
        
        self.logger.info(f"Generating IGV snapshots for sample: {sample_id}")
        
        snapshots = {}
        
        try:
            # Create IGV batch script
            batch_script = self._create_batch_script(bam_file, sample_id)
            
            # Run IGV in batch mode
            success = self._run_igv_batch(batch_script)
            
            if success:
                # Collect generated snapshots
                for region_name in self.regions.keys():
                    snapshot_file = os.path.join(
                        self.snapshot_dir, 
                        f"{sample_id}_{region_name}.png"
                    )
                    if os.path.exists(snapshot_file):
                        snapshots[region_name] = snapshot_file
                    else:
                        self.logger.warning(f"Snapshot not found: {snapshot_file}")
                        snapshots[region_name] = None
            else:
                self.logger.error(f"IGV batch execution failed for sample: {sample_id}")
            
            # Clean up batch script
            if os.path.exists(batch_script):
                os.remove(batch_script)
                
        except Exception as e:
            self.logger.error(f"Error generating snapshots for {sample_id}: {str(e)}")
        
        return snapshots
    
    def _create_batch_script(self, bam_file: str, sample_id: str) -> str:
        """Create IGV batch script for automated snapshot generation"""
        
        # Create temporary batch script
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bat', delete=False) as f:
            batch_script = f.name
            
            # IGV batch commands
            commands = [
                f"new",  # New session
                f"genome hg38",  # Set genome
                f"load {bam_file}",  # Load BAM file
                "",  # Empty line for safety
                "snapshotDirectory " + self.snapshot_dir,  # Set snapshot directory
                "maxPanelHeight 1000",  # Set panel height
                ""
            ]
            
            # Add commands for each region
            for region_name, region_coord in self.regions.items():
                commands.extend([
                    f"goto {region_coord}",  # Navigate to region
                    "sort position",  # Sort by position
                    "collapse",  # Collapse tracks
                    f"snapshot {sample_id}_{region_name}.png",  # Take snapshot
                    ""
                ])
            
            commands.append("exit")  # Exit IGV
            
            # Write commands to file
            f.write("\n".join(commands))
        
        self.logger.debug(f"Created IGV batch script: {batch_script}")
        return batch_script
    
    def _run_igv_batch(self, batch_script: str) -> bool:
        """Execute IGV in batch mode"""
        
        try:
            # Construct IGV command
            cmd = [
                self.igv_path,
                '-b', batch_script,
                '--memory', self.igv_memory
            ]
            
            self.logger.debug(f"Running IGV command: {' '.join(cmd)}")
            
            # Run IGV
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if process.returncode == 0:
                self.logger.info("IGV batch execution completed successfully")
                return True
            else:
                self.logger.error(f"IGV batch execution failed: {process.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            self.logger.error("IGV batch execution timed out")
            return False
        except Exception as e:
            self.logger.error(f"Error running IGV batch: {str(e)}")
            return False
    
    def generate_custom_snapshot(self, bam_file: str, region: str, output_file: str) -> bool:
        """Generate a custom snapshot for a specific region"""
        
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bat', delete=False) as f:
                batch_script = f.name
                
                commands = [
                    "new",
                    "genome hg38",
                    f"load {bam_file}",
                    f"snapshotDirectory {os.path.dirname(output_file)}",
                    f"goto {region}",
                    "sort position",
                    "collapse", 
                    f"snapshot {os.path.basename(output_file)}",
                    "exit"
                ]
                
                f.write("\n".join(commands))
            
            # Run IGV
            success = self._run_igv_batch(batch_script)
            
            # Clean up
            os.remove(batch_script)
            
            return success and os.path.exists(output_file)
            
        except Exception as e:
            self.logger.error(f"Error generating custom snapshot: {str(e)}")
            return False
    
    def validate_snapshots(self, snapshots: Dict) -> Dict:
        """Validate generated snapshots and extract basic metrics"""
        
        validated = {}
        
        for region_name, snapshot_file in snapshots.items():
            if snapshot_file and os.path.exists(snapshot_file):
                # Get file size as basic validation
                file_size = os.path.getsize(snapshot_file)
                
                validated[region_name] = {
                    'file_path': snapshot_file,
                    'file_size': file_size,
                    'exists': True,
                    'valid': file_size > 1000  # Minimum expected size
                }
            else:
                validated[region_name] = {
                    'file_path': snapshot_file,
                    'exists': False,
                    'valid': False
                }
        
        return validated
    
    def create_igv_session(self, bam_file: str, sample_id: str) -> str:
        """Create an IGV session file for manual review"""
        
        session_file = os.path.join(self.snapshot_dir, f"{sample_id}_session.xml")
        
        session_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg38" hasGeneTrack="true" hasSequenceTrack="true" version="8">
    <Resources>
        <Resource path="{bam_file}"/>
    </Resources>
    <Panel height="400" name="Panel1" width="1200">
        <Track altColor="0,0,178" autoScale="true" color="175,175,175" 
               colorScale="ContinuousColorScale;0.0;308.0;255,255,255;175,175,175" 
               displayMode="COLLAPSED" featureVisibilityWindow="-1" 
               fontSize="10" id="{bam_file}" name="{sample_id}" 
               showReference="false" sortable="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" 
                      maximum="308.0" minimum="0.0" type="LINEAR"/>
        </Track>
    </Panel>
    <Panel height="60" name="FeaturePanel" width="1200">
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" 
               displayMode="COLLAPSED" featureVisibilityWindow="-1" 
               fontSize="10" height="35" id="Reference sequence" 
               name="Reference sequence" showReference="false" 
               sortable="false" visible="true"/>
    </Panel>
    <Panel height="60" name="Panel2" width="1200">
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" 
               displayMode="COLLAPSED" featureVisibilityWindow="-1" 
               fontSize="10" height="35" id="Gene" name="Gene" 
               showReference="false" sortable="false" visible="true"/>
    </Panel>
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>"""
        
        try:
            with open(session_file, 'w') as f:
                f.write(session_xml)
            
            self.logger.info(f"Created IGV session file: {session_file}")
            return session_file
            
        except Exception as e:
            self.logger.error(f"Error creating IGV session: {str(e)}")
            return None
    
    def batch_snapshot_analysis(self, snapshots: Dict) -> Dict:
        """Analyze snapshots for visual CNV indicators"""
        
        analysis_results = {}
        
        for region_name, snapshot_info in snapshots.items():
            if isinstance(snapshot_info, dict) and snapshot_info.get('valid', False):
                # Basic image analysis could be added here
                # For now, just record that the snapshot is available
                analysis_results[region_name] = {
                    'snapshot_available': True,
                    'file_size': snapshot_info['file_size'],
                    'visual_quality': 'good' if snapshot_info['file_size'] > 5000 else 'poor'
                }
            else:
                analysis_results[region_name] = {
                    'snapshot_available': False,
                    'visual_quality': 'unavailable'
                }
        
        return analysis_results
