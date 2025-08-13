#!/usr/bin/env python3
"""
Remote processing script - Run all AlphaFold predictions and subsequent processing steps on remote server
Please place this script in the ~/student/students_webserver/zhijing/ directory on the remote server

Modified to generate both Z-scores and aligned PDB files from DALI comparisons
"""

import os
import sys
import subprocess
import argparse
import time
import glob
import shutil
import re

class RemotePipeline:
    def __init__(self, json_filename):
        self.json_filename = json_filename
        self.output_name = os.path.splitext(json_filename)[0]
        self.base_dir = os.path.expanduser("~/student/students_webserver/zhijing")
        self.input_dir = os.path.join(self.base_dir, "input_jsons")
        self.output_dir = os.path.join(self.base_dir, "predicted_structures")
        
    def run_alphafold_prediction(self):
        """Run AlphaFold prediction"""
        print(f"🔬 Starting AlphaFold prediction...")
        
        # Ensure we're in the correct directory
        os.chdir(self.base_dir)
        
        # Build command
        cmd = [
            "protenix", "predict",
            "--input", f"{self.input_dir}/{self.json_filename}",
            "--out_dir", f"./{self.output_name}_output",
            "--seeds", "101",
            "--use_msa_server"
        ]
        
        print(f"Executing command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"❌ AlphaFold prediction failed:")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")
            return False
        
        print(f"✅ AlphaFold prediction completed")
        return True
    
    def convert_cif_to_pdb(self):
        """Convert CIF file to PDB file"""
        print(f"🔄 Converting CIF file to PDB file...")
        
        # The actual output path is sequence_1_output instead of predicted_structures
        actual_output_dir = f"{self.base_dir}/{self.output_name}_output"
        cif_pattern = f"{actual_output_dir}/{self.output_name}/seed_101/predictions/{self.output_name}_seed_101_sample_0.cif"
        
        print(f"Looking for CIF file: {cif_pattern}")
        
        if not os.path.exists(cif_pattern):
            print(f"❌ CIF file not found: {cif_pattern}")
            # Try to find with wildcard
            cif_files = glob.glob(f"{actual_output_dir}/{self.output_name}/seed_101/predictions/*sample_0.cif")
            if cif_files:
                cif_pattern = cif_files[0]
                print(f"Found CIF file: {cif_pattern}")
            else:
                return False
        
        # Generate output PDB filename to predicted_structures directory
        os.makedirs(self.output_dir, exist_ok=True)
        pdb_file = f"{self.output_dir}/{self.output_name}_sample_0.pdb"
        
        # Run conversion script
        cmd = ["python", "convert_cif_to_pdb.py", cif_pattern, pdb_file]
        print(f"Executing command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"❌ CIF to PDB conversion failed:")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")
            return False
        
        print(f"✅ CIF to PDB conversion completed")
        print(result.stdout)
        return True
    
    def convert_pdb_to_dat(self):
        """Convert PDB file to DAT file"""
        print(f"🔄 Converting PDB file to DAT file...")
        
        # Look for the converted PDB file
        pdb_file = f"{self.output_dir}/{self.output_name}_sample_0.pdb"
        
        if not os.path.exists(pdb_file):
            print(f"❌ PDB file does not exist: {pdb_file}")
            return False
        
        print(f"Found PDB file: {pdb_file}")
        
        # Run conversion script (your script will automatically process all PDB files in predicted_structures directory)
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
        return True
    
    def run_dali_comparison(self):
        """Run DALI structure comparison and generate both Z-scores and aligned PDB files"""
        print(f"🔍 Starting DALI structure comparison...")
        
        # Ensure we're in the correct directory
        os.chdir(self.base_dir)
        
        # Create output directories
        os.makedirs("dali_results", exist_ok=True)
        os.makedirs("dali_alignments", exist_ok=True)  # 新增：存储对齐后的PDB文件
        
        # Get query structure name
        query_name = (self.output_name[:4].lower() + "xxx")[:4]
        
        # Check for actual files in query_structures_DAT directory
        dat_files_in_query = []
        if os.path.exists("query_structures_DAT"):
            dat_files_in_query = [f for f in os.listdir("query_structures_DAT") if f.endswith('.dat')]
        
        if dat_files_in_query:
            actual_query_name = os.path.splitext(dat_files_in_query[0])[0]
            print(f"Query structure name found in query_structures_DAT: {actual_query_name}")
            query_name = actual_query_name
        else:
            print(f"⚠️ No DAT files found in query_structures_DAT directory, using default name: {query_name}")
        
        print(f"Using query structure name: {query_name}")
        
        # Check DAT directory
        dat_dir = "DAT"
        if not os.path.exists(dat_dir):
            print(f"❌ DAT directory does not exist: {dat_dir}")
            return False
        
        # Get all DAT files
        dat_files = [f for f in os.listdir(dat_dir) if f.endswith('.dat')]
        if not dat_files:
            print(f"❌ No .dat files found in DAT directory")
            return False
        
        print(f"Found {len(dat_files)} reference structures")
        
        # Run DALI comparison for each DAT file
        successful_comparisons = 0
        
        for dat_file in dat_files:
            ref_id = os.path.splitext(dat_file)[0]
            print(f"🔗 Comparing {query_name} vs {ref_id}")
            
            # 修改DALI命令以生成更详细的输出
            cmd = [
                "/home/wenhao/6tx0/software/dali/DaliLite.v5/bin/dali.pl",
                "--cd1", query_name,
                "--cd2", ref_id,
                "--dat1", "query_structures_DAT",
                "--dat2", "DAT",
                "--outfmt", "full",  # 改为full格式以获得更多信息
                "--pairwise"  # 生成成对对齐
            ]
            
            # 为每个比较创建单独的工作目录
            comparison_dir = f"dali_temp_{query_name}_vs_{ref_id}"
            os.makedirs(comparison_dir, exist_ok=True)
            
            # 在临时目录中运行DALI
            original_dir = os.getcwd()
            
            try:
                os.chdir(comparison_dir)
                
                # 复制必要的文件到临时目录
                shutil.copy(f"../query_structures_DAT/{query_name}.dat", f"{query_name}.dat")
                shutil.copy(f"../DAT/{ref_id}.dat", f"{ref_id}.dat")
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                os.chdir(original_dir)
                
                if result.returncode != 0:
                    print(f"⚠️ DALI comparison failed for {ref_id}: {result.stderr}")
                    continue
                
                # 处理DALI输出文件
                self.process_dali_output(comparison_dir, query_name, ref_id)
                
                successful_comparisons += 1
                print(f"✅ Successfully compared {ref_id}")
                
            except Exception as e:
                os.chdir(original_dir)
                print(f"⚠️ Error processing {ref_id}: {str(e)}")
                continue
            finally:
                # 清理临时目录
                if os.path.exists(comparison_dir):
                    try:
                        shutil.rmtree(comparison_dir)
                    except:
                        pass  # 忽略清理错误
        
        print(f"✅ DALI comparison completed, successfully compared {successful_comparisons} structures")
        
        # 生成汇总报告
        if successful_comparisons > 0:
            self.generate_dali_summary()
        
        return successful_comparisons > 0
    
    def process_dali_output(self, temp_dir, query_name, ref_id):
        """处理单个DALI比较的输出文件"""
        
        # DALI生成的主要输出文件
        summary_file = f"{temp_dir}/{query_name}.txt"
        
        # 目标路径
        target_summary = f"dali_results/{query_name}_vs_{ref_id}.txt"
        target_alignment = f"dali_alignments/{query_name}_vs_{ref_id}_aligned.pdb"
        
        # 移动摘要文件
        if os.path.exists(summary_file):
            shutil.move(summary_file, target_summary)
            print(f"📄 Summary saved to: {target_summary}")
        
        # 查找并移动对齐文件 - DALI可能生成多种命名格式的对齐文件
        possible_alignment_files = [
            f"{temp_dir}/{query_name}-{ref_id}.pdb",
            f"{temp_dir}/{query_name}_{ref_id}.pdb",
            f"{temp_dir}/alignment.pdb",
            f"{temp_dir}/{query_name}.pdb",
            f"{temp_dir}/pairwise.pdb",
            f"{temp_dir}/align.pdb"
        ]
        
        alignment_found = False
        
        # 首先尝试找到明确的对齐文件
        for align_file in possible_alignment_files:
            if os.path.exists(align_file):
                shutil.copy(align_file, target_alignment)  # 使用copy而不是move，避免丢失原文件
                print(f"🧬 Aligned structure saved to: {target_alignment}")
                alignment_found = True
                break
        
        # 如果没有找到，列出所有PDB文件并尝试识别
        if not alignment_found:
            if os.path.exists(temp_dir):
                all_files = os.listdir(temp_dir)
                pdb_files = [f for f in all_files if f.endswith('.pdb')]
                
                if pdb_files:
                    # 取第一个PDB文件作为对齐结果
                    source_pdb = f"{temp_dir}/{pdb_files[0]}"
                    shutil.copy(source_pdb, target_alignment)
                    print(f"🧬 Aligned structure (auto-detected) saved to: {target_alignment}")
                    alignment_found = True
                else:
                    print(f"⚠️ No PDB alignment file found for {ref_id}")
                    print(f"Files in temp directory: {all_files}")
        
        return alignment_found
    
    def generate_dali_summary(self):
        """生成DALI结果汇总"""
        print(f"📊 Generating DALI summary...")
        
        summary_data = []
        query_name_short = self.output_name[:4].lower()
        
        # 读取所有结果文件
        if os.path.exists("dali_results"):
            for result_file in os.listdir("dali_results"):
                if result_file.endswith(".txt") and not result_file.endswith("_summary.tsv"):
                    file_path = f"dali_results/{result_file}"
                    try:
                        zscore, rmsd, align_len, lali, identity = self.parse_dali_result(file_path)
                        
                        # 从文件名提取参考结构ID
                        ref_id = result_file.replace(f"{query_name_short}xxx_vs_", "").replace(".txt", "")
                        if ref_id == result_file:  # 如果替换失败，尝试其他模式
                            ref_id = result_file.replace(".txt", "").split("_vs_")[-1]
                        
                        # 检查对应的对齐文件是否存在
                        alignment_file = f"dali_alignments/{query_name_short}xxx_vs_{ref_id}_aligned.pdb"
                        has_alignment = os.path.exists(alignment_file)
                        
                        summary_data.append({
                            'reference_id': ref_id,
                            'zscore': zscore,
                            'rmsd': rmsd,
                            'alignment_length': align_len,
                            'lali': lali,
                            'identity': identity,
                            'result_file': result_file,
                            'alignment_pdb': f"{query_name_short}xxx_vs_{ref_id}_aligned.pdb" if has_alignment else "N/A",
                            'has_alignment': has_alignment
                        })
                    except Exception as e:
                        print(f"⚠️ Error parsing {result_file}: {str(e)}")
        
        if not summary_data:
            print(f"⚠️ No valid DALI results found to summarize")
            return
        
        # 按Z-score排序
        summary_data.sort(key=lambda x: float(x['zscore']) if x['zscore'] != 'N/A' and x['zscore'].replace('.','').replace('-','').isdigit() else -999, reverse=True)
        
        # 写入汇总文件
        summary_file = f"dali_results/{self.output_name}_summary.tsv"
        with open(summary_file, 'w') as f:
            f.write("Reference_ID\tZ-score\tRMSD\tAlignment_Length\tLali\tIdentity\tResult_File\tAlignment_PDB\tHas_Alignment\n")
            for item in summary_data:
                f.write(f"{item['reference_id']}\t{item['zscore']}\t{item['rmsd']}\t{item['alignment_length']}\t{item['lali']}\t{item['identity']}\t{item['result_file']}\t{item['alignment_pdb']}\t{item['has_alignment']}\n")
        
        print(f"📋 Summary saved to: {summary_file}")
        
        # 显示前10个最佳匹配
        print(f"\n🏆 Top 10 matches:")
        print(f"{'Rank':<4} {'Reference':<12} {'Z-score':<8} {'RMSD':<8} {'Lali':<6} {'Identity':<10} {'Alignment'}")
        print("-" * 70)
        
        for i, item in enumerate(summary_data[:10]):
            alignment_status = "✓" if item['has_alignment'] else "✗"
            print(f"{i+1:2d}.  {item['reference_id']:<12} {item['zscore']:<8} {item['rmsd']:<8} {item['lali']:<6} {item['identity']:<10} {alignment_status}")
        
        # 统计信息
        total_comparisons = len(summary_data)
        with_alignments = len([x for x in summary_data if x['has_alignment']])
        significant_matches = len([x for x in summary_data if x['zscore'] != 'N/A' and float(x['zscore']) > 2.0])
        
        print(f"\n📈 Statistics:")
        print(f"Total comparisons: {total_comparisons}")
        print(f"With alignment files: {with_alignments}")
        print(f"Significant matches (Z-score > 2.0): {significant_matches}")
        
        return summary_file
    
    def parse_dali_result(self, result_file):
        """解析DALI结果文件，提取Z-score、RMSD等信息"""
        zscore = "N/A"
        rmsd = "N/A"
        align_len = "N/A"
        lali = "N/A"
        identity = "N/A"
        
        try:
            with open(result_file, 'r') as f:
                content = f.read()
                
                # 查找Z-score (通常在"Z-score"或"Z="后面)
                zscore_patterns = [
                    r'Z[=-]?\s*([0-9]+\.?[0-9]*)',
                    r'Z-score[:\s]*([0-9]+\.?[0-9]*)',
                    r'Z\s*=\s*([0-9]+\.?[0-9]*)'
                ]
                
                for pattern in zscore_patterns:
                    zscore_match = re.search(pattern, content, re.IGNORECASE)
                    if zscore_match:
                        zscore = zscore_match.group(1)
                        break
                
                # 查找RMSD
                rmsd_patterns = [
                    r'RMSD[=:\s]*([0-9]+\.?[0-9]*)',
                    r'rmsd[=:\s]*([0-9]+\.?[0-9]*)',
                    r'RMS[=:\s]*([0-9]+\.?[0-9]*)'
                ]
                
                for pattern in rmsd_patterns:
                    rmsd_match = re.search(pattern, content, re.IGNORECASE)
                    if rmsd_match:
                        rmsd = rmsd_match.group(1)
                        break
                
                # 查找对齐长度 (Lali)
                lali_patterns = [
                    r'Lali[=:\s]*([0-9]+)',
                    r'lali[=:\s]*([0-9]+)',
                    r'aligned\s+residues[:\s]*([0-9]+)',
                    r'alignment\s+length[:\s]*([0-9]+)'
                ]
                
                for pattern in lali_patterns:
                    lali_match = re.search(pattern, content, re.IGNORECASE)
                    if lali_match:
                        lali = lali_match.group(1)
                        align_len = lali  # 使用lali作为对齐长度
                        break
                
                # 查找序列一致性
                identity_patterns = [
                    r'Identity[=:\s]*([0-9]+\.?[0-9]*%?)',
                    r'identity[=:\s]*([0-9]+\.?[0-9]*%?)',
                    r'Id[=:\s]*([0-9]+\.?[0-9]*%?)'
                ]
                
                for pattern in identity_patterns:
                    identity_match = re.search(pattern, content, re.IGNORECASE)
                    if identity_match:
                        identity = identity_match.group(1)
                        break
                
        except Exception as e:
            print(f"Error parsing {result_file}: {str(e)}")
        
        return zscore, rmsd, align_len, lali, identity
    
    def extract_results(self):
        """Extract Z-score results using existing script"""
        print(f"📊 Extracting results using existing script...")
        
        # Check if dali_results directory exists
        if not os.path.exists("dali_results"):
            print(f"❌ dali_results directory does not exist, please run DALI comparison first")
            return False
        
        cmd = ["python", "extract_zscores.py"]
        print(f"Executing command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"⚠️ Original result extraction script failed (this is expected with new format):")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")
            print(f"✅ Using built-in summary generation instead")
            return True  # 不视为失败，因为我们已经有了自己的汇总功能
        
        print(f"✅ Original result extraction completed")
        print(result.stdout)
        return True
    
    def check_prerequisites(self):
        """Check required files and directories"""
        print(f"🔍 Checking environment...")
        
        # Check input file
        input_file = f"{self.input_dir}/{self.json_filename}"
        if not os.path.exists(input_file):
            print(f"❌ Input file does not exist: {input_file}")
            return False
        
        # Check script files
        required_scripts = [
            "convert_cif_to_pdb.py",
            "import_pdbs_to_dat.py"
        ]
        
        for script in required_scripts:
            if not os.path.exists(script):
                print(f"❌ Missing required script: {script}")
                return False
        
        # Check directories
        required_dirs = [
            "DAT",
            "query_structures_DAT"
        ]
        
        for dir_name in required_dirs:
            if not os.path.exists(dir_name):
                print(f"❌ Missing required directory: {dir_name}")
                return False
        
        # Create output directories if they don't exist
        os.makedirs("predicted_structures", exist_ok=True)
        os.makedirs("dali_results", exist_ok=True)
        os.makedirs("dali_alignments", exist_ok=True)
        
        print(f"✅ Environment check passed")
        return True
    
    def debug_find_output_files(self):
        """Debug: Find prediction output files"""
        print(f"🔍 Finding prediction output files...")
        
        base_dir = os.path.expanduser("~/student/students_webserver/zhijing")
        
        debug_commands = [
            f"find {base_dir} -name '*{self.output_name}*' -type d",
            f"find {base_dir} -name '*{self.output_name}*' -type f",
            f"find {base_dir} -name '*.cif' | grep {self.output_name}",
            f"ls -la {base_dir}/predicted_structures/",
            f"find {base_dir} -name '*sequence_1*' | head -20",
            f"find {base_dir} -name '*.cif' | head -20"
        ]
        
        for cmd in debug_commands:
            print(f"Executing: {cmd}")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            print(f"Output: {result.stdout.strip()}")
            if result.stderr:
                print(f"Error: {result.stderr.strip()}")
            print("-" * 40)
    
    def run_pipeline(self):
        """Run complete remote processing pipeline"""
        print(f"🚀 Starting remote processing pipeline...")
        print(f"📋 JSON file: {self.json_filename}")
        print(f"📋 Output name: {self.output_name}")
        print(f"📋 Working directory: {self.base_dir}")
        print("="*50)
        
        try:
            # Step 1: Check environment
            if not self.check_prerequisites():
                return False
            
            # Step 2: Run AlphaFold prediction
            if not self.run_alphafold_prediction():
                return False
            
            # Step 3: Convert CIF to PDB
            if not self.convert_cif_to_pdb():
                return False
            
            # Step 4: Convert PDB to DAT
            if not self.convert_pdb_to_dat():
                return False
            
            # Step 5: Run DALI comparison (now generates both Z-scores and alignments)
            if not self.run_dali_comparison():
                return False
            
            # Step 6: Extract results (optional, as we now have built-in summary)
            self.extract_results()  # Don't fail if this doesn't work
            
            print("="*50)
            print("🎉 All steps completed!")
            print(f"📁 Results saved in:")
            print(f"   - dali_results/{self.output_name}_summary.tsv (summary)")
            print(f"   - dali_results/ (detailed comparison results)")
            print(f"   - dali_alignments/ (aligned PDB structures)")
            
        except Exception as e:
            print(f"❌ Pipeline execution failed: {str(e)}")
            return False
        
        return True


def main():
    parser = argparse.ArgumentParser(description='Remote AlphaFold prediction and structure comparison pipeline')
    parser.add_argument('json_filename', help='JSON filename (not path, just filename)')
    parser.add_argument('--check', '-c', action='store_true', help='Only check environment, do not run pipeline')
    parser.add_argument('--find-files', '-f', action='store_true', help='Find output files')
    parser.add_argument('--skip-prediction', '-s', action='store_true', help='Skip prediction step, directly process existing files')
    parser.add_argument('--dali-only', '-d', action='store_true', help='Only run DALI comparison (assumes PDB and DAT files exist)')
    
    args = parser.parse_args()
    
    pipeline = RemotePipeline(args.json_filename)
    
    if args.check:
        print("🔍 Check mode")
        success = pipeline.check_prerequisites()
        if success:
            print("✅ Environment check passed, pipeline can be run")
        return
    
    if args.find_files:
        print("🔍 Find files mode")
        pipeline.debug_find_output_files()
        return
    
    if args.dali_only:
        print("🔍 DALI comparison only mode")
        try:
            if not pipeline.run_dali_comparison():
                print("💥 DALI comparison failed!")
                sys.exit(1)
            
            pipeline.extract_results()  # Optional
            print("🎉 DALI processing completed!")
            
        except Exception as e:
            print(f"❌ DALI processing failed: {str(e)}")
            sys.exit(1)
        return
    
    if args.skip_prediction:
        print("⏭️ Skipping prediction step")
        try:
            # Skip prediction, start directly from conversion
            if not pipeline.convert_cif_to_pdb():
                print("💥 CIF to PDB conversion failed!")
                sys.exit(1)
            
            if not pipeline.convert_pdb_to_dat():
                print("💥 PDB to DAT conversion failed!")
                sys.exit(1)
            
            if not pipeline.run_dali_comparison():
                print("💥 DALI comparison failed!")
                sys.exit(1)
            
            pipeline.extract_results()  # Optional
            print("🎉 Processing completed!")
            
        except Exception as e:
            print(f"❌ Processing failed: {str(e)}")
            sys.exit(1)
        return
    
    # Run complete pipeline
    success = pipeline.run_pipeline()
    
    if success:
        print("\n🎊 Pipeline executed successfully!")
        sys.exit(0)
    else:
        print("\n💥 Pipeline execution failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()