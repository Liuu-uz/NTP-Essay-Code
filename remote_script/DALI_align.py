#!/usr/bin/env python3
"""
Standalone script to convert PDB to DAT and run DALI comparison
This version preserves HTML intermediate files instead of looking for PDB outputs
Usage: python pdb_to_dali.py <pdb_filename>
Example: python pdb_to_dali.py sequence_1_sample_0.pdb
"""

import os
import sys
import subprocess
import argparse
import glob
import shutil
import re
from pathlib import Path

class PDBToDaliProcessor:
    def __init__(self, pdb_filename, max_refs=None):
        self.pdb_filename = pdb_filename
        self.pdb_basename = os.path.splitext(pdb_filename)[0]
        self.max_refs = max_refs
        self.base_dir = os.path.expanduser("~/student/students_webserver/zhijing")
        self.predicted_structures_dir = os.path.join(self.base_dir, "predicted_structures")
        self.results_dir = os.path.join(self.base_dir, "single_dali_results")
        
        # Create results directory
        os.makedirs(self.results_dir, exist_ok=True)
        
        # Ensure we're working in the correct directory
        os.chdir(self.base_dir)
        
        print(f"🔧 Initialized PDBToDaliProcessor:")
        print(f"   PDB file: {self.pdb_filename}")
        print(f"   PDB basename: {self.pdb_basename}")
        print(f"   Base directory: {self.base_dir}")
        print(f"   Results directory: {self.results_dir}")
        
    def convert_pdb_to_dat(self):
        """Convert PDB file to DAT file using existing script"""
        print(f"🔄 Converting PDB file to DAT file...")
        
        # Check if PDB file exists in predicted_structures directory
        pdb_path = os.path.join(self.predicted_structures_dir, self.pdb_filename)
        if not os.path.exists(pdb_path):
            print(f"❌ PDB file not found: {pdb_path}")
            return False
        
        print(f"Found PDB file: {pdb_path}")
        
        # Backup other PDB files temporarily to ensure only our target file is processed
        temp_dir = f"{self.base_dir}/temp_pdb_backup"
        os.makedirs(temp_dir, exist_ok=True)
        
        # Backup other PDB files
        other_pdb_files = []
        for f in glob.glob(f"{self.predicted_structures_dir}/*.pdb"):
            if f != pdb_path:
                other_pdb_files.append(f)
                backup_name = os.path.basename(f)
                shutil.move(f, f"{temp_dir}/{backup_name}")
                print(f"  Temporarily moved: {backup_name}")
        
        try:
            # Run conversion script (only processes our target PDB file)
            cmd = ["python", "import_pdbs_to_dat.py"]
            print(f"Executing command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"❌ PDB to DAT conversion failed:")
                print(f"stdout: {result.stdout}")
                print(f"stderr: {result.stderr}")
                return False
            
            print(f"✅ PDB to DAT conversion completed")
            print(result.stdout)
            
        finally:
            # Restore other PDB files
            for backup_file in glob.glob(f"{temp_dir}/*.pdb"):
                original_name = os.path.basename(backup_file)
                shutil.move(backup_file, f"{self.predicted_structures_dir}/{original_name}")
                print(f"  Restored: {original_name}")
            
            # Remove temporary directory
            if os.path.exists(temp_dir):
                try:
                    os.rmdir(temp_dir)
                except:
                    pass
        
        return True
    
    def get_query_structure_name(self):
        """Get the actual query structure name from DAT files"""
        dat_files_in_query = []
        if os.path.exists("query_structures_DAT"):
            dat_files_in_query = [f for f in os.listdir("query_structures_DAT") if f.endswith('.dat')]
        
        if dat_files_in_query:
            # Use the most recently created DAT file
            actual_query_name = os.path.splitext(dat_files_in_query[-1])[0]
            print(f"Query structure name found: {actual_query_name}")
            return actual_query_name
        else:
            print(f"❌ No DAT files found in query_structures_DAT directory")
            return None
    
    def run_dali_comparison(self):
        """Run DALI structure comparison and save text results and HTML files"""
        print(f"🔍 Starting DALI structure comparison...")
        
        # Get query structure name
        query_name = self.get_query_structure_name()
        if not query_name:
            return False
        
        print(f"Using query structure name: {query_name}")
        
        # Check DAT directory
        dat_dir = "DAT"
        if not os.path.exists(dat_dir):
            print(f"❌ DAT directory does not exist: {dat_dir}")
            return False
        
        # Get all reference DAT files
        ref_dat_files = [f for f in os.listdir(dat_dir) if f.endswith('.dat')]
        if not ref_dat_files:
            print(f"❌ No .dat files found in DAT directory")
            return False
        
        print(f"Found {len(ref_dat_files)} reference structures")
        
        # Create sequence-specific dali results directory with proper path handling
        sequence_dali_dir = os.path.join(self.results_dir, f"{self.pdb_basename}_dali_results")
        sequence_html_dir = os.path.join(self.results_dir, f"{self.pdb_basename}_html_results")
        os.makedirs(sequence_dali_dir, exist_ok=True)
        os.makedirs(sequence_html_dir, exist_ok=True)
        
        print(f"📁 Output directories:")
        print(f"   Text results: {sequence_dali_dir}")
        print(f"   HTML results: {sequence_html_dir}")
        
        # Run DALI comparison for each DAT file
        successful_comparisons = 0
        html_files_generated = 0
        
        # Use all reference DAT files (or limit to max_refs if specified)
        if self.max_refs and self.max_refs < len(ref_dat_files):
            test_refs = ref_dat_files[:self.max_refs]
            print(f"Limiting comparison to first {self.max_refs} references")
        else:
            test_refs = ref_dat_files  # Compare against all references
            print(f"Comparing against all {len(test_refs)} reference structures")
        
        for i, dat_file in enumerate(test_refs):
            ref_id = os.path.splitext(dat_file)[0]
            print(f"\n🔗 [{i+1}/{len(test_refs)}] Comparing {query_name} vs {ref_id}")
            
            # First run: Generate summary (for Z-scores)
            cmd_summary = [
                "/home/wenhao/6tx0/software/dali/DaliLite.v5/bin/dali.pl",
                "--cd1", query_name,
                "--cd2", ref_id,
                "--dat1", "query_structures_DAT",
                "--dat2", "DAT",
                "--outfmt", "summary",
                "--clean"
            ]
            
            print(f"   Running summary command...")
            result_summary = subprocess.run(cmd_summary, capture_output=True, text=True)
            
            # Move summary result file
            default_output = f"{query_name}.txt"
            target_output = os.path.join(sequence_dali_dir, f"{query_name}_vs_{ref_id}.txt")
            
            if os.path.exists(default_output):
                try:
                    shutil.move(default_output, target_output)
                    successful_comparisons += 1
                    print(f"   ✓ Summary saved: {os.path.basename(target_output)}")
                except Exception as e:
                    print(f"   ⚠️ Failed to move summary file: {e}")
            else:
                print(f"   ⚠️ No summary file generated")
            
            # Second run: Generate HTML output (without --clean to preserve intermediate files)
            print(f"   Running HTML generation command...")
            cmd_html = [
                "/home/wenhao/6tx0/software/dali/DaliLite.v5/bin/dali.pl",
                "--cd1", query_name,
                "--cd2", ref_id,
                "--dat1", "query_structures_DAT",
                "--dat2", "DAT"
                # No --clean to preserve HTML and other intermediate files
            ]
            
            result_html = subprocess.run(cmd_html, capture_output=True, text=True)
            
            print(f"   HTML command return code: {result_html.returncode}")
            if result_html.stdout:
                print(f"   HTML stdout: {result_html.stdout[:200]}...")
            if result_html.stderr:
                print(f"   HTML stderr: {result_html.stderr[:200]}...")
            
            # List all files in current directory to see what was generated
            print(f"   Files in current directory after DALI:")
            current_files = [f for f in os.listdir('.') if not os.path.isdir(f)]
            for f in current_files:
                if query_name.lower() in f.lower() or ref_id.lower() in f.lower():
                    print(f"     - {f}")
            
            # Look for generated HTML files and other intermediate files
            possible_output_files = []
            
            # Common DALI output file patterns
            file_patterns = [
                f"{query_name}.html",
                f"{query_name}-{ref_id}.html", 
                f"{query_name}_{ref_id}.html",
                f"{ref_id}.html",
                f"{ref_id}-{query_name}.html",
                f"{ref_id}_{query_name}.html",
                "alignment.html",
                "align.html",
                "output.html",
                "result.html",
                f"{query_name}.ps",  # PostScript files
                f"{query_name}.log", # Log files
                f"{query_name}.out", # Output files
                f"{query_name}.ali", # Alignment files
                f"{query_name}.pdb", # In case PDB files are generated
            ]
            
            # Also check for any files that contain the query or ref names
            all_files = os.listdir('.')
            for file in all_files:
                if (query_name.lower() in file.lower() or 
                    ref_id.lower() in file.lower()):
                    # Skip already processed text files
                    if not (file.endswith('.txt') and file in [f"{query_name}.txt"]):
                        if file not in possible_output_files:
                            possible_output_files.append(file)
            
            # Add standard patterns
            for pattern in file_patterns:
                if pattern not in possible_output_files:
                    possible_output_files.append(pattern)
            
            print(f"   Looking for output files: {possible_output_files}")
            
            # Save all found intermediate files
            files_found = []
            for output_file in possible_output_files:
                if os.path.exists(output_file):
                    print(f"   📋 Found intermediate file: {output_file}")
                    
                    # Determine target directory based on file type
                    if output_file.endswith('.html'):
                        target_dir = sequence_html_dir
                        file_suffix = "_results.html"
                    else:
                        target_dir = sequence_html_dir  # Keep all intermediate files together
                        file_suffix = f"_intermediate.{output_file.split('.')[-1]}"
                    
                    target_file = os.path.join(target_dir, f"{query_name}_vs_{ref_id}{file_suffix}")
                    
                    try:
                        shutil.copy2(output_file, target_file)  # Use copy2 to preserve original
                        files_found.append(output_file)
                        if output_file.endswith('.html'):
                            html_files_generated += 1
                        print(f"   ✓ Intermediate file saved: {os.path.basename(target_file)}")
                    except Exception as e:
                        print(f"   ⚠️ Failed to copy intermediate file: {e}")
            
            if not files_found:
                print(f"   ℹ️ No intermediate files found for {ref_id}")
                
                # Try alternative DALI command
                print(f"   Trying alternative DALI command...")
                cmd_alt = [
                    "/home/wenhao/6tx0/software/dali/DaliLite.v5/bin/dali.pl",
                    "--cd1", query_name,
                    "--cd2", ref_id,
                    "--dat1", "query_structures_DAT", 
                    "--dat2", "DAT",
                    "--outfmt", "html"  # Try explicit HTML format
                ]
                
                result_alt = subprocess.run(cmd_alt, capture_output=True, text=True)
                print(f"   Alternative command return code: {result_alt.returncode}")
                
                # Check again for output files
                for output_file in possible_output_files:
                    if os.path.exists(output_file) and output_file not in files_found:
                        print(f"   📋 Found file (alternative): {output_file}")
                        
                        if output_file.endswith('.html'):
                            target_dir = sequence_html_dir
                            file_suffix = "_results.html"
                        else:
                            target_dir = sequence_html_dir
                            file_suffix = f"_intermediate.{output_file.split('.')[-1]}"
                        
                        target_file = os.path.join(target_dir, f"{query_name}_vs_{ref_id}{file_suffix}")
                        
                        try:
                            shutil.copy2(output_file, target_file)
                            files_found.append(output_file)
                            if output_file.endswith('.html'):
                                html_files_generated += 1
                            print(f"   ✓ File saved (alternative): {os.path.basename(target_file)}")
                        except Exception as e:
                            print(f"   ⚠️ Failed to copy file (alternative): {e}")
            
            print(f"   📊 Files collected for this comparison: {len(files_found)}")
        
        print(f"\n✅ DALI comparison completed:")
        print(f"   - {successful_comparisons} text summaries generated")
        print(f"   - {html_files_generated} HTML files generated")
        print(f"   - Results saved in: {sequence_dali_dir}")
        print(f"   - HTML files saved in: {sequence_html_dir}")
        
        return successful_comparisons > 0
    
    def extract_results(self):
        """Extract Z-score results (same logic as batch script)"""
        print(f"📊 Extracting Z-score results...")
        
        sequence_dali_dir = f"{self.results_dir}/{self.pdb_basename}_dali_results"
        
        if not os.path.exists(sequence_dali_dir):
            print(f"❌ DALI results directory not found: {sequence_dali_dir}")
            return False
        
        # Use the same Z-score extraction logic as batch script
        zscores = {}
        
        def parse_z_score(file_path):
            """Parse Z-score from DALI result file (same as batch script)"""
            try:
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
            except:
                return None
        
        # Parse all result files
        for fname in os.listdir(sequence_dali_dir):
            if fname.endswith(".txt"):
                fpath = os.path.join(sequence_dali_dir, fname)
                z = parse_z_score(fpath)
                if z is not None:
                    zscores[fname] = z
        
        # Sort results by Z-score (descending)
        sorted_results = sorted(zscores.items(), key=lambda x: -x[1])
        
        # Save to CSV file (same format as batch script)
        csv_file = f"{self.results_dir}/{self.pdb_basename}_zscores.csv"
        
        if sorted_results:
            print(f"🎯 Found {len(sorted_results)} Z-score results")
            
            best = sorted_results[0]
            print(f"🏆 Best match: {best[0]} with Z = {best[1]:.2f}")
            
            # Save as CSV
            with open(csv_file, "w") as f:
                f.write("filename,Z-score\n")
                for fname, z in sorted_results:
                    f.write(f"{fname},{z:.2f}\n")
            
            print(f"✅ Z-scores saved to: {csv_file}")
            
            # Display top 10 results
            print(f"\n🏆 Top 10 matches:")
            print(f"{'Rank':<4} {'Reference':<20} {'Z-score':<8}")
            print("-" * 35)
            
            for i, (fname, z) in enumerate(sorted_results[:10]):
                # Extract reference ID from filename
                ref_id = fname.replace(f"{self.get_query_structure_name()}_vs_", "").replace(".txt", "")
                print(f"{i+1:2d}.  {ref_id:<20} {z:<8.2f}")
            
            return True
        else:
            print(f"❌ No valid Z-score results found")
            # Create empty CSV file
            with open(csv_file, "w") as f:
                f.write("filename,Z-score\n")
                f.write("No matches found,0.0\n")
            return False
    
    def generate_summary_report(self):
        """Generate a comprehensive summary report"""
        print(f"📋 Generating summary report...")
        
        # Check what files were generated
        csv_file = f"{self.results_dir}/{self.pdb_basename}_zscores.csv"
        sequence_dali_dir = f"{self.results_dir}/{self.pdb_basename}_dali_results"
        sequence_html_dir = f"{self.results_dir}/{self.pdb_basename}_html_results"
        
        # Count results
        total_comparisons = 0
        html_files = 0
        intermediate_files = 0
        
        if os.path.exists(sequence_dali_dir):
            total_comparisons = len([f for f in os.listdir(sequence_dali_dir) if f.endswith('.txt')])
        if os.path.exists(sequence_html_dir):
            all_files = os.listdir(sequence_html_dir)
            html_files = len([f for f in all_files if f.endswith('.html')])
            intermediate_files = len([f for f in all_files if not f.endswith('.html')])
        
        # Read CSV results
        significant_matches = 0
        high_scoring_matches = 0
        best_zscore = 0.0
        best_match = "None"
        
        if os.path.exists(csv_file):
            try:
                with open(csv_file, 'r') as f:
                    lines = f.readlines()[1:]  # Skip header
                    for line in lines:
                        if line.strip() and "No matches found" not in line:
                            parts = line.strip().split(',')
                            if len(parts) >= 2:
                                zscore = float(parts[1])
                                if zscore > 2.0:
                                    significant_matches += 1
                                if zscore > 9.0:
                                    high_scoring_matches += 1
                                if zscore > best_zscore:
                                    best_zscore = zscore
                                    best_match = parts[0]
            except:
                pass
        
        # Generate summary
        summary_file = f"{self.results_dir}/{self.pdb_basename}_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write(f"DALI Structure Comparison Summary\n")
            f.write(f"=" * 40 + "\n\n")
            f.write(f"Query structure: {self.pdb_filename}\n")
            f.write(f"Query DAT name: {self.get_query_structure_name() or 'Unknown'}\n")
            f.write(f"Analysis date: {os.popen('date').read().strip()}\n\n")
            
            f.write(f"Results Statistics:\n")
            f.write(f"- Total comparisons: {total_comparisons}\n")
            f.write(f"- HTML files generated: {html_files}\n")
            f.write(f"- Other intermediate files: {intermediate_files}\n")
            f.write(f"- Significant matches (Z > 2.0): {significant_matches}\n")
            f.write(f"- High-scoring matches (Z > 9.0): {high_scoring_matches}\n")
            f.write(f"- Best match: {best_match} (Z = {best_zscore:.2f})\n\n")
            
            f.write(f"Output Files:\n")
            f.write(f"- Detailed results: {sequence_dali_dir}/\n")
            f.write(f"- HTML & intermediate files: {sequence_html_dir}/\n")
            f.write(f"- Z-score CSV: {csv_file}\n")
            f.write(f"- This summary: {summary_file}\n\n")
            
            f.write(f"File Locations:\n")
            f.write(f"- Text summaries: {sequence_dali_dir}/*_vs_*.txt\n")
            f.write(f"- HTML results: {sequence_html_dir}/*_vs_*_results.html\n")
            f.write(f"- Intermediate files: {sequence_html_dir}/*_vs_*_intermediate.*\n")
        
        print(f"📄 Summary report saved to: {summary_file}")
        print(f"📊 Generated {html_files} HTML files and {intermediate_files} intermediate files")
        return summary_file
    
    def cleanup_temp_files(self):
        """Clean up temporary files (but preserve HTML and intermediate files)"""
        print(f"🗑️ Cleaning up temporary files...")
        
        # Remove query DAT files to clean up for next run
        if os.path.exists("query_structures_DAT"):
            for dat_file in glob.glob("query_structures_DAT/*.dat"):
                try:
                    os.remove(dat_file)
                    print(f"  Removed: {dat_file}")
                except:
                    pass
        
        # Remove leftover DALI temporary files in working directory
        # But be more careful - only remove files we're sure are temporary
        query_name = self.get_query_structure_name()
        if query_name:
            temp_patterns = [
                f"{query_name}.tmp",
                f"{query_name}.log",
                "*.tmp",
                "dali*.log"
            ]
            
            for pattern in temp_patterns:
                for temp_file in glob.glob(pattern):
                    try:
                        os.remove(temp_file)
                        print(f"  Removed temp file: {temp_file}")
                    except:
                        pass
    
    def check_prerequisites(self):
        """Check if required files and directories exist"""
        print(f"🔍 Checking prerequisites...")
        
        # Check if we're in the right directory
        if not os.path.exists("import_pdbs_to_dat.py"):
            print(f"❌ import_pdbs_to_dat.py not found. Are you in the correct directory?")
            print(f"Current directory: {os.getcwd()}")
            return False
        
        # Check if PDB file exists
        pdb_path = os.path.join(self.predicted_structures_dir, self.pdb_filename)
        if not os.path.exists(pdb_path):
            print(f"❌ PDB file not found: {pdb_path}")
            return False
        
        # Check for required directories
        if not os.path.exists("DAT"):
            print(f"❌ DAT directory not found")
            return False
        
        # Check DALI executable
        dali_executable = "/home/wenhao/6tx0/software/dali/DaliLite.v5/bin/dali.pl"
        if not os.path.exists(dali_executable):
            print(f"❌ DALI executable not found: {dali_executable}")
            return False
        
        # Create required directories
        os.makedirs("query_structures_DAT", exist_ok=True)
        
        print(f"✅ Prerequisites check passed")
        return True
    
    def run_full_pipeline(self):
        """Run the complete pipeline from PDB to DALI results"""
        print(f"🚀 Starting PDB to DALI pipeline...")
        print(f"📋 PDB file: {self.pdb_filename}")
        print(f"📂 Working directory: {self.base_dir}")
        print(f"📁 Results directory: {self.results_dir}")
        print("="*60)
        
        try:
            # Step 1: Check prerequisites
            if not self.check_prerequisites():
                return False
            
            # Step 2: Convert PDB to DAT
            if not self.convert_pdb_to_dat():
                print(f"❌ PDB to DAT conversion failed")
                return False
            
            # Step 3: Run DALI comparison (preserving HTML and intermediate files)
            if not self.run_dali_comparison():
                print(f"❌ DALI comparison failed")
                return False
            
            # Step 4: Extract Z-score results (same logic as batch script)
            self.extract_results()
            
            # Step 5: Generate summary report
            summary_file = self.generate_summary_report()
            
            # Step 6: Clean up temporary files (but keep HTML/intermediate files)
            self.cleanup_temp_files()
            
            print("="*60)
            print("🎉 Pipeline completed successfully!")
            print(f"📁 Results saved in:")
            print(f"   - {self.results_dir}/{self.pdb_basename}_zscores.csv (Z-scores)")
            print(f"   - {self.results_dir}/{self.pdb_basename}_dali_results/ (text summaries)")
            print(f"   - {self.results_dir}/{self.pdb_basename}_html_results/ (HTML & intermediate files)")
            print(f"   - {summary_file} (summary report)")
            
            # Show file counts
            html_dir = f"{self.results_dir}/{self.pdb_basename}_html_results"
            if os.path.exists(html_dir):
                all_files = os.listdir(html_dir)
                html_count = len([f for f in all_files if f.endswith('.html')])
                other_count = len([f for f in all_files if not f.endswith('.html')])
                print(f"📊 Generated {html_count} HTML files and {other_count} other intermediate files")
                
                if len(all_files) > 0:
                    print(f"\n🧬 HTML & intermediate files location:")
                    print(f"   {html_dir}/")
                    print(f"   HTML files: {self.get_query_structure_name() or 'query'}_vs_[reference]_results.html")
                    print(f"   Intermediate files: {self.get_query_structure_name() or 'query'}_vs_[reference]_intermediate.*")
            
            return True
            
        except Exception as e:
            print(f"❌ Pipeline execution failed: {str(e)}")
            return False


def main():
    parser = argparse.ArgumentParser(
        description='Convert PDB to DAT and run DALI comparison (preserving HTML intermediate files)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pdb_to_dali.py sequence_1_sample_0.pdb     # Run full pipeline
  python pdb_to_dali.py --check sequence_1_sample_0.pdb  # Check prerequisites only
  python pdb_to_dali.py --convert-only sequence_1_sample_0.pdb  # Convert PDB to DAT only
  python pdb_to_dali.py --dali-only sequence_1_sample_0.pdb     # Run DALI only
        """
    )
    
    parser.add_argument('pdb_filename', help='PDB filename (should be in predicted_structures directory)')
    parser.add_argument('--check', '-c', action='store_true', help='Only check prerequisites')
    parser.add_argument('--convert-only', action='store_true', help='Only convert PDB to DAT')
    parser.add_argument('--dali-only', action='store_true', help='Only run DALI (assumes DAT file exists)')
    parser.add_argument('--no-cleanup', action='store_true', help='Skip cleanup of temporary files')
    parser.add_argument('--max-refs', type=int, default=None, help='Maximum number of reference structures to compare (default: all)')
    
    args = parser.parse_args()
    
    processor = PDBToDaliProcessor(args.pdb_filename, args.max_refs)
    
    if args.check:
        print("🔍 Check mode")
        success = processor.check_prerequisites()
        if success:
            print("✅ All prerequisites satisfied")
        else:
            print("❌ Prerequisites check failed")
        return
    
    if args.convert_only:
        print("🔄 Convert only mode")
        if processor.check_prerequisites() and processor.convert_pdb_to_dat():
            print("✅ PDB to DAT conversion completed")
        else:
            print("❌ PDB to DAT conversion failed")
        return
    
    if args.dali_only:
        print("🔍 DALI only mode")
        try:
            if processor.run_dali_comparison():
                processor.extract_results()
                processor.generate_summary_report()
                if not args.no_cleanup:
                    processor.cleanup_temp_files()
                print("✅ DALI comparison completed")
            else:
                print("❌ DALI comparison failed")
        except Exception as e:
            print(f"❌ DALI comparison failed: {str(e)}")
        return
    
    # Run full pipeline
    success = processor.run_full_pipeline()
    
    if success:
        print("\n🎊 All operations completed successfully!")
        sys.exit(0)
    else:
        print("\n💥 Pipeline execution failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()