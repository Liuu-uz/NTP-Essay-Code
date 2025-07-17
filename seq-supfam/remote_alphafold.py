#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simplified Remote AlphaFold Predictor
Run ProteiniX prediction on remote server, convert CIF to PDB, and download results
"""

import os
import sys
import subprocess
import time

class SimplifiedAlphaFoldPredictor:
    def __init__(self, sequence_id):
        self.sequence_id = sequence_id
        self.json_filename = f"{sequence_id}.json"
        
        # Fixed configuration
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.upload_dir = "~/student/students_webserver/zhijing/input_jsons"
        self.work_dir = "~/student/students_webserver/zhijing"
        
        print(f"✅ Initialized predictor:")
        print(f"   JSON file: {self.json_filename}")
        print(f"   Sequence ID: {self.sequence_id}")
    
    def run_prediction_and_conversion(self):
        """Run prediction and conversion interactively via SSH"""
        output_dir = f"{self.sequence_id}_output"
        
        print(f"\n🚀 Starting AlphaFold prediction and conversion...")
        print("="*60)
        print(f"📋 JSON file: {self.upload_dir}/{self.json_filename}")
        print(f"📋 Sequence ID: {self.sequence_id}")
        print(f"📋 Output directory: {output_dir}")
        print("="*60)
        print("🔄 SSH session will open. Commands to run:")
        print()
        print(f"   cd {self.work_dir}")
        print(f"   conda activate protoneix")
        print(f"   protenix predict --input {self.upload_dir}/{self.json_filename} --out_dir ./{output_dir} --seeds 101")
        print()
        print("   # After prediction completes, run conversion:")
        print(f"   python convert_cif_to_pdb.py {output_dir}/{self.sequence_id}/seed_101/predictions/{self.sequence_id}_seed_101_sample_0.cif predicted_structures/{self.sequence_id}_sample_0.pdb")
        print()
        print("⚠️  After both steps complete, type 'exit' to continue with download")
        print("="*60)
        
        input("Press Enter to start SSH session...")
        
        # Start interactive SSH session
        ssh_cmd = f"ssh -t {self.remote_server} 'cd {self.work_dir}; bash'"
        
        print(f"\n🔗 Connecting to remote server...")
        subprocess.run(ssh_cmd, shell=True)
        
        print(f"\n✅ SSH session ended")
        return output_dir
    
    def download_pdb_results(self):
        """Download PDB files from predicted_structures directory"""
        print(f"\n📥 Downloading PDB results...")
        
        # Create local results directory
        local_results_dir = f"alphafold_results"
        if not os.path.exists(local_results_dir):
            os.makedirs(local_results_dir)
            print(f"📁 Created results directory: {local_results_dir}")
        
        # Download PDB files from predicted_structures
        print(f"   Downloading PDB files from predicted_structures/...")
        pdb_pattern = f"predicted_structures/{self.sequence_id}_sample_*.pdb"
        remote_pdb_path = f"{self.remote_server}:{self.work_dir}/{pdb_pattern}"
        
        cmd = f"scp {remote_pdb_path} {local_results_dir}/"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        downloaded_files = []
        
        if result.returncode == 0:
            print(f"   ✅ Downloaded PDB files successfully")
            downloaded_files.append(f"{local_results_dir}/{self.sequence_id}_sample_0.pdb")
            
            # Check downloaded file
            pdb_file = f"{local_results_dir}/{self.sequence_id}_sample_0.pdb"
            if os.path.exists(pdb_file):
                file_size = os.path.getsize(pdb_file)
                print(f"   📄 {pdb_file} ({file_size} bytes)")
                
                # Show first few lines
                with open(pdb_file, 'r') as f:
                    first_line = f.readline().strip()
                    print(f"   📝 First line: {first_line}")
        else:
            print(f"   ❌ Failed to download PDB files")
            print(f"   Error: {result.stderr}")
            
            # Try downloading entire predicted_structures directory
            print(f"   Trying to download entire predicted_structures directory...")
            remote_dir = f"{self.remote_server}:{self.work_dir}/predicted_structures"
            local_dir = os.path.join(local_results_dir, "predicted_structures")
            
            cmd = f"scp -r {remote_dir} {local_results_dir}/"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"   ✅ Downloaded predicted_structures directory")
                downloaded_files.append(local_dir)
        
        return downloaded_files
    
    def process_sequence(self):
        """Main function to process sequence"""
        print(f"🧬 AlphaFold Predictor for {self.sequence_id}")
        print("="*50)
        
        try:
            # 1. Run prediction and conversion via SSH
            output_dir = self.run_prediction_and_conversion()
            
            # 2. Download PDB results
            downloaded_files = self.download_pdb_results()
            
            # 3. Summary
            print(f"\n{'='*60}")
            print("🎉 Process completed!")
            
            if downloaded_files:
                print(f"📋 Downloaded files:")
                for file_path in downloaded_files:
                    print(f"   📄 {file_path}")
                print(f"   🧬 PDB files ready for molecular visualization")
            else:
                print("⚠️  No PDB files downloaded")
                print("💡 Please check if prediction and conversion completed successfully")
            
            print("="*60)
            return True
            
        except Exception as e:
            print(f"❌ Error during processing: {e}")
            return False

def main():
    if len(sys.argv) != 2:
        print("Usage:")
        print("  python simplified_alphafold.py SEQUENCE_ID")
        print()
        print("Example:")
        print("  python simplified_alphafold.py sequence_1")
        return
    
    sequence_id = sys.argv[1]
    
    try:
        predictor = SimplifiedAlphaFoldPredictor(sequence_id)
        predictor.process_sequence()
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    main()