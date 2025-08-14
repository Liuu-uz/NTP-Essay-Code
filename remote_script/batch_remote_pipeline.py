#!/usr/bin/env python3
"""
Batch Remote Processing Script - Run all AlphaFold predictions and subsequent processing steps on remote server
Please place this script in the remote server's ~/student/students_webserver/zhijing/ directory
"""

import os
import sys
import subprocess
import argparse
import time
import glob
import json
from pathlib import Path

class BatchRemotePipeline:
    def __init__(self):
        self.base_dir = os.path.expanduser("~/student/students_webserver/zhijing")
        self.input_dir = os.path.join(self.base_dir, "input_jsons")
        self.output_dir = os.path.join(self.base_dir, "predicted_structures")
        self.results_dir = os.path.join(self.base_dir, "batch_results")
        
        # Create result directories
        os.makedirs(self.results_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)
        
    def get_json_files(self):
        """Get all JSON files"""
        if not os.path.exists(self.input_dir):
            print(f"ERROR: Input directory does not exist: {self.input_dir}")
            return []
        
        json_files = glob.glob(os.path.join(self.input_dir, "*.json"))
        
        if not json_files:
            print(f"ERROR: No JSON files found in input directory: {self.input_dir}")
            return []
        
        # Return only filenames, not full paths
        json_filenames = [os.path.basename(f) for f in json_files]
        print(f"Found {len(json_filenames)} JSON files")
        return json_filenames
    
    def run_single_alphafold_prediction(self, json_filename):
        """Run single AlphaFold prediction"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"Starting AlphaFold prediction: {json_filename}")
        
        # Ensure we're in the correct directory
        os.chdir(self.base_dir)
        
        # Build command
        cmd = [
            "protenix", "predict",
            "--input", f"{self.input_dir}/{json_filename}",
            "--out_dir", f"./{output_name}_output",
            "--seeds", "101",
            "--use_msa_server"
        ]
        
        print(f"Executing command: {' '.join(cmd)}")
        print("=" * 50)
        
        # Use real-time output mode
        try:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                bufsize=1
            )
            
            # Display output in real-time
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    print(f"  {output.strip()}")
            
            # Wait for process completion
            return_code = process.poll()
            
            if return_code != 0:
                print(f"ERROR: AlphaFold prediction failed {json_filename} (return code: {return_code})")
                return False
            
            print("=" * 50)
            print(f"SUCCESS: AlphaFold prediction completed: {json_filename}")
            return True
            
        except Exception as e:
            print(f"ERROR: AlphaFold prediction error {json_filename}: {str(e)}")
            return False
    
    def convert_single_cif_to_pdb(self, json_filename):
        """Convert single CIF file to PDB file"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"Converting CIF to PDB: {json_filename}")
        
        # Search for CIF file
        actual_output_dir = f"{self.base_dir}/{output_name}_output"
        cif_pattern = f"{actual_output_dir}/{output_name}/seed_101/predictions/{output_name}_seed_101_sample_0.cif"
        
        if not os.path.exists(cif_pattern):
            # Try searching with wildcards
            cif_files = glob.glob(f"{actual_output_dir}/{output_name}/seed_101/predictions/*sample_0.cif")
            if cif_files:
                cif_pattern = cif_files[0]
                print(f"  Found CIF file: {cif_pattern}")
            else:
                print(f"ERROR: CIF file not found: {json_filename}")
                return False
        else:
            print(f"  Found CIF file: {cif_pattern}")
        
        # Generate output PDB filename
        pdb_file = f"{self.output_dir}/{output_name}_sample_0.pdb"
        
        # Run conversion script
        cmd = ["python", "convert_cif_to_pdb.py", cif_pattern, pdb_file]
        print(f"  Executing: {' '.join(cmd)}")
        
        try:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True
            )
            
            stdout, _ = process.communicate()
            
            if process.returncode != 0:
                print(f"ERROR: CIF to PDB conversion failed {json_filename}:")
                print(f"  Output: {stdout}")
                return False
            
            print(f"  Conversion output: {stdout.strip()}")
            print(f"SUCCESS: CIF to PDB conversion completed: {json_filename}")
            return True
            
        except Exception as e:
            print(f"ERROR: CIF to PDB conversion error {json_filename}: {str(e)}")
            return False
    
    def convert_single_pdb_to_dat(self, json_filename):
        """Convert single PDB file to DAT file"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"Converting PDB to DAT: {json_filename}")
        
        # Check if PDB file exists
        pdb_file = f"{self.output_dir}/{output_name}_sample_0.pdb"
        if not os.path.exists(pdb_file):
            print(f"ERROR: PDB file does not exist: {json_filename}")
            return False
        
        # Temporarily move other PDB files, process only current file
        temp_dir = f"{self.base_dir}/temp_pdb_backup"
        os.makedirs(temp_dir, exist_ok=True)
        
        # Backup other PDB files
        other_pdb_files = []
        for f in glob.glob(f"{self.output_dir}/*.pdb"):
            if f != pdb_file:
                other_pdb_files.append(f)
                backup_name = os.path.basename(f)
                os.rename(f, f"{temp_dir}/{backup_name}")
        
        # Run conversion script (only process current PDB file)
        cmd = ["python", "import_pdbs_to_dat.py"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Restore other PDB files
        for backup_file in glob.glob(f"{temp_dir}/*.pdb"):
            original_name = os.path.basename(backup_file)
            os.rename(backup_file, f"{self.output_dir}/{original_name}")
        
        if result.returncode != 0:
            print(f"ERROR: PDB to DAT conversion failed {json_filename}: {result.stderr}")
            return False
        
        print(f"SUCCESS: PDB to DAT conversion completed: {json_filename}")
        return True
    
    def run_single_dali_comparison(self, json_filename):
        """Run DALI structure comparison for single sequence"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"Starting DALI structure comparison: {json_filename}")
        
        # Ensure we're in the correct directory
        os.chdir(self.base_dir)
        
        # Create sequence-specific dali_results directory
        sequence_dali_dir = f"{self.results_dir}/{output_name}_dali_results"
        os.makedirs(sequence_dali_dir, exist_ok=True)
        
        # Find generated DAT files
        dat_files_in_query = []
        if os.path.exists("query_structures_DAT"):
            dat_files_in_query = [f for f in os.listdir("query_structures_DAT") if f.endswith('.dat')]
        
        if not dat_files_in_query:
            print(f"ERROR: No query structure DAT files found: {json_filename}")
            return False
        
        # Use the latest generated DAT file as query structure
        query_name = os.path.splitext(dat_files_in_query[-1])[0]  # Use the last file
        print(f"Using query structure: {query_name}")
        
        # Check reference DAT directory
        dat_dir = "DAT"
        if not os.path.exists(dat_dir):
            print(f"ERROR: DAT directory does not exist: {dat_dir}")
            return False
        
        # Get all reference DAT files
        ref_dat_files = [f for f in os.listdir(dat_dir) if f.endswith('.dat')]
        if not ref_dat_files:
            print(f"ERROR: No .dat files found in reference DAT directory")
            return False
        
        print(f"Found {len(ref_dat_files)} reference structures")
        
        # Run DALI comparison for each reference DAT file
        successful_comparisons = 0
        for dat_file in ref_dat_files:
            ref_id = os.path.splitext(dat_file)[0]
            
            cmd = [
                "/home/wenhao/6tx0/software/dali/DaliLite.v5/bin/dali.pl",
                "--cd1", query_name,
                "--cd2", ref_id,
                "--dat1", "query_structures_DAT",
                "--dat2", "DAT",
                "--outfmt", "summary",
                "--clean"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                continue
            
            # Move result files to sequence-specific directory
            default_output = f"{query_name}.txt"
            target_output = f"{sequence_dali_dir}/{query_name}_vs_{ref_id}.txt"
            
            if os.path.exists(default_output):
                try:
                    os.rename(default_output, target_output)
                    successful_comparisons += 1
                except Exception as e:
                    print(f"WARNING: Failed to move output file: {e}")
        
        print(f"SUCCESS: DALI comparison completed {json_filename}: {successful_comparisons} successful comparisons")
        return successful_comparisons > 0
    
    def extract_single_results(self, json_filename):
        """Extract Z-score results for single sequence"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"Extracting results: {json_filename}")
        
        sequence_dali_dir = f"{self.results_dir}/{output_name}_dali_results"
        
        if not os.path.exists(sequence_dali_dir):
            print(f"ERROR: DALI results directory does not exist: {json_filename}")
            return False
        
        # Modified extract_zscores.py logic, specifically for this sequence
        zscores = {}
        
        def parse_z_score(file_path):
            with open(file_path) as f:
                for line in f:
                    if line.strip().startswith("1:"):
                        parts = line.strip().split()
                        try:
                            z = float(parts[2])  # Third field is Z-score
                            return z
                        except:
                            return None
            return None
        
        # Iterate through this sequence's results directory
        for fname in os.listdir(sequence_dali_dir):
            if fname.endswith(".txt"):
                fpath = os.path.join(sequence_dali_dir, fname)
                z = parse_z_score(fpath)
                if z is not None:
                    zscores[fname] = z
        
        # Sort output
        sorted_results = sorted(zscores.items(), key=lambda x: -x[1])
        
        # Save to sequence-specific CSV file
        csv_file = f"{self.results_dir}/{output_name}_zscores.csv"
        
        if sorted_results:
            print(f"Found {len(sorted_results)} Z-score results")
            
            best = sorted_results[0]
            print(f"Best match: {best[0]} with Z = {best[1]:.2f}")
            
            # Save as CSV
            with open(csv_file, "w") as f:
                f.write("filename,Z-score\n")
                for fname, z in sorted_results:
                    f.write(f"{fname},{z:.2f}\n")
            
            print(f"SUCCESS: Z-scores saved to: {csv_file}")
            return True
        else:
            print(f"ERROR: No valid Z-score results found: {json_filename}")
            # Create empty CSV file to indicate processing completed but no results
            with open(csv_file, "w") as f:
                f.write("filename,Z-score\n")
                f.write("No matches found,0.0\n")
            return False
    
    def cleanup_intermediate_files(self, json_filename):
        """Clean up intermediate files to save space"""
        output_name = os.path.splitext(json_filename)[0]
        
        # Delete AlphaFold output directory (keep PDB files)
        alphafold_output_dir = f"{self.base_dir}/{output_name}_output"
        if os.path.exists(alphafold_output_dir):
            try:
                subprocess.run(["rm", "-rf", alphafold_output_dir], check=True)
                print(f"Cleaned AlphaFold output directory: {output_name}")
            except:
                print(f"WARNING: Cleanup failed: {alphafold_output_dir}")
        
        # Clean temporary files in query_structures_DAT
        if os.path.exists("query_structures_DAT"):
            for dat_file in glob.glob("query_structures_DAT/*.dat"):
                try:
                    os.remove(dat_file)
                except:
                    pass
    
    def process_single_sequence(self, json_filename):
        """Process complete pipeline for single sequence"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"\n{'='*60}")
        print(f"Processing sequence: {json_filename}")
        print(f"Output name: {output_name}")
        print(f"{'='*60}")
        
        try:
            # Step 1: AlphaFold prediction
            if not self.run_single_alphafold_prediction(json_filename):
                return False
            
            # Step 2: Convert CIF to PDB
            if not self.convert_single_cif_to_pdb(json_filename):
                return False
            
            # Step 3: Convert PDB to DAT
            if not self.convert_single_pdb_to_dat(json_filename):
                return False
            
            # Step 4: DALI comparison
            if not self.run_single_dali_comparison(json_filename):
                return False
            
            # Step 5: Extract results
            if not self.extract_single_results(json_filename):
                print(f"WARNING: Z-score extraction failed, but continuing processing")
            
            # Step 6: Clean up intermediate files
            self.cleanup_intermediate_files(json_filename)
            
            print(f"SUCCESS: Sequence processing completed: {json_filename}")
            return True
            
        except Exception as e:
            print(f"ERROR: Sequence processing failed {json_filename}: {str(e)}")
            return False
    
    def run_batch_pipeline(self):
        """Run batch processing pipeline"""
        print(f"Starting batch remote processing pipeline...")
        print(f"Working directory: {self.base_dir}")
        print(f"Input directory: {self.input_dir}")
        print(f"Results directory: {self.results_dir}")
        print("="*60)
        
        # Get all JSON files
        json_files = self.get_json_files()
        if not json_files:
            print("ERROR: No JSON files found")
            return False
        
        # Statistics variables
        total_files = len(json_files)
        successful_processes = 0
        failed_processes = 0
        
        print(f"Starting to process {total_files} sequences...")
        
        # Process sequences one by one
        for i, json_filename in enumerate(json_files, 1):
            print(f"\n[{i}/{total_files}] Currently processing: {json_filename}")
            
            start_time = time.time()
            
            if self.process_single_sequence(json_filename):
                successful_processes += 1
                elapsed = time.time() - start_time
                print(f"SUCCESS: Sequence completed ({elapsed:.1f}s): {json_filename}")
            else:
                failed_processes += 1
                elapsed = time.time() - start_time
                print(f"ERROR: Sequence failed ({elapsed:.1f}s): {json_filename}")
            
            # Show progress
            remaining = total_files - i
            if remaining > 0:
                print(f"Progress: {i}/{total_files} completed, {remaining} remaining")
        
        # Final statistics
        print("\n" + "="*60)
        print("Batch processing statistics:")
        print(f"   Total sequences: {total_files}")
        print(f"   Successfully processed: {successful_processes}")
        print(f"   Failed count: {failed_processes}")
        print(f"   Success rate: {successful_processes/total_files*100:.1f}%")
        
        # Show result files
        csv_files = glob.glob(f"{self.results_dir}/*_zscores.csv")
        print(f"\nGenerated result files: {len(csv_files)}")
        print(f"Results directory: {self.results_dir}")
        
        if csv_files:
            print("SUCCESS: Batch processing completed!")
            print(f"Next steps:")
            print(f"   1. Check results: ls -la {self.results_dir}/")
            print(f"   2. Download results: scp -r webserver@coulomb.phys.ucl.ac.uk:{self.results_dir}/ ./")
            return True
        else:
            print("ERROR: No result files generated")
            return False

def main():
    parser = argparse.ArgumentParser(
        description='Batch remote AlphaFold prediction and structure comparison pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  python batch_remote_pipeline.py                 # Batch process all JSON files
  python batch_remote_pipeline.py --check         # Check environment
  python batch_remote_pipeline.py --results-only  # Extract existing results only
        """
    )
    
    parser.add_argument('--check', '-c', action='store_true', 
                       help='Only check environment, do not run pipeline')
    parser.add_argument('--results-only', '-r', action='store_true', 
                       help='Only extract results, skip prediction steps')
    
    args = parser.parse_args()
    
    pipeline = BatchRemotePipeline()
    
    if args.check:
        print("Environment check mode")
        # Basic environment check
        required_scripts = [
            "convert_cif_to_pdb.py",
            "import_pdbs_to_dat.py"
        ]
        
        all_good = True
        for script in required_scripts:
            if os.path.exists(script):
                print(f"SUCCESS: Found script: {script}")
            else:
                print(f"ERROR: Missing script: {script}")
                all_good = False
        
        if os.path.exists(pipeline.input_dir):
            json_count = len(glob.glob(os.path.join(pipeline.input_dir, "*.json")))
            print(f"SUCCESS: Input directory exists, contains {json_count} JSON files")
        else:
            print(f"ERROR: Input directory does not exist: {pipeline.input_dir}")
            all_good = False
        
        if all_good:
            print("SUCCESS: Environment check passed, can run batch pipeline")
        else:
            print("ERROR: Environment check failed")
        return
    
    if args.results_only:
        print("Results extraction only mode")
        # Here you can add logic for extracting results only
        print("Feature under development...")
        return
    
    # Run batch pipeline
    success = pipeline.run_batch_pipeline()
    
    if success:
        print("\nSUCCESS: Batch pipeline executed successfully!")
        sys.exit(0)
    else:
        print("\nERROR: Batch pipeline execution failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()