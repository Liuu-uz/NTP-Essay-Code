#!/usr/bin/env python3
"""
批量远程处理脚本 - 在远程服务器上批量运行所有AlphaFold预测和后续处理步骤
请将此脚本放在远程服务器的 ~/student/students_webserver/zhijing/ 目录下
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
        
        # 创建结果目录
        os.makedirs(self.results_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)
        
    def get_json_files(self):
        """获取所有JSON文件"""
        if not os.path.exists(self.input_dir):
            print(f"❌ 输入目录不存在: {self.input_dir}")
            return []
        
        json_files = glob.glob(os.path.join(self.input_dir, "*.json"))
        
        if not json_files:
            print(f"❌ 在输入目录中没有找到JSON文件: {self.input_dir}")
            return []
        
        # 只返回文件名，不是完整路径
        json_filenames = [os.path.basename(f) for f in json_files]
        print(f"📁 找到 {len(json_filenames)} 个JSON文件")
        return json_filenames
    
    def run_single_alphafold_prediction(self, json_filename):
        """运行单个AlphaFold预测"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"🔬 开始AlphaFold预测: {json_filename}")
        
        # 确保在正确的目录
        os.chdir(self.base_dir)
        
        # 构建命令
        cmd = [
            "protenix", "predict",
            "--input", f"{self.input_dir}/{json_filename}",
            "--out_dir", f"./{output_name}_output",
            "--seeds", "101",
            "--use_msa_server"
        ]
        
        print(f"执行命令: {' '.join(cmd)}")
        print("=" * 50)
        
        # 使用实时输出模式
        try:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                bufsize=1
            )
            
            # 实时显示输出
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    print(f"  {output.strip()}")
            
            # 等待进程完成
            return_code = process.poll()
            
            if return_code != 0:
                print(f"❌ AlphaFold预测失败 {json_filename} (返回码: {return_code})")
                return False
            
            print("=" * 50)
            print(f"✅ AlphaFold预测完成: {json_filename}")
            return True
            
        except Exception as e:
            print(f"❌ AlphaFold预测异常 {json_filename}: {str(e)}")
            return False
    
    def convert_single_cif_to_pdb(self, json_filename):
        """转换单个CIF文件为PDB文件"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"🔄 转换CIF为PDB: {json_filename}")
        
        # 查找CIF文件
        actual_output_dir = f"{self.base_dir}/{output_name}_output"
        cif_pattern = f"{actual_output_dir}/{output_name}/seed_101/predictions/{output_name}_seed_101_sample_0.cif"
        
        if not os.path.exists(cif_pattern):
            # 尝试用通配符查找
            cif_files = glob.glob(f"{actual_output_dir}/{output_name}/seed_101/predictions/*sample_0.cif")
            if cif_files:
                cif_pattern = cif_files[0]
                print(f"  找到CIF文件: {cif_pattern}")
            else:
                print(f"❌ 找不到CIF文件: {json_filename}")
                return False
        else:
            print(f"  找到CIF文件: {cif_pattern}")
        
        # 生成输出PDB文件名
        pdb_file = f"{self.output_dir}/{output_name}_sample_0.pdb"
        
        # 运行转换脚本
        cmd = ["python", "convert_cif_to_pdb.py", cif_pattern, pdb_file]
        print(f"  执行: {' '.join(cmd)}")
        
        try:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True
            )
            
            stdout, _ = process.communicate()
            
            if process.returncode != 0:
                print(f"❌ CIF转PDB失败 {json_filename}:")
                print(f"  输出: {stdout}")
                return False
            
            print(f"  转换输出: {stdout.strip()}")
            print(f"✅ CIF转PDB完成: {json_filename}")
            return True
            
        except Exception as e:
            print(f"❌ CIF转PDB异常 {json_filename}: {str(e)}")
            return False
    
    def convert_single_pdb_to_dat(self, json_filename):
        """转换单个PDB文件为DAT文件"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"🔄 转换PDB为DAT: {json_filename}")
        
        # 检查PDB文件是否存在
        pdb_file = f"{self.output_dir}/{output_name}_sample_0.pdb"
        if not os.path.exists(pdb_file):
            print(f"❌ PDB文件不存在: {json_filename}")
            return False
        
        # 临时移动其他PDB文件，只处理当前文件
        temp_dir = f"{self.base_dir}/temp_pdb_backup"
        os.makedirs(temp_dir, exist_ok=True)
        
        # 备份其他PDB文件
        other_pdb_files = []
        for f in glob.glob(f"{self.output_dir}/*.pdb"):
            if f != pdb_file:
                other_pdb_files.append(f)
                backup_name = os.path.basename(f)
                os.rename(f, f"{temp_dir}/{backup_name}")
        
        # 运行转换脚本（只处理当前PDB文件）
        cmd = ["python", "import_pdbs_to_dat.py"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # 恢复其他PDB文件
        for backup_file in glob.glob(f"{temp_dir}/*.pdb"):
            original_name = os.path.basename(backup_file)
            os.rename(backup_file, f"{self.output_dir}/{original_name}")
        
        if result.returncode != 0:
            print(f"❌ PDB转DAT失败 {json_filename}: {result.stderr}")
            return False
        
        print(f"✅ PDB转DAT完成: {json_filename}")
        return True
    
    def run_single_dali_comparison(self, json_filename):
        """运行单个序列的DALI结构对比"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"🔍 开始DALI结构对比: {json_filename}")
        
        # 确保在正确的目录
        os.chdir(self.base_dir)
        
        # 创建该序列专用的dali_results目录
        sequence_dali_dir = f"{self.results_dir}/{output_name}_dali_results"
        os.makedirs(sequence_dali_dir, exist_ok=True)
        
        # 查找生成的DAT文件
        dat_files_in_query = []
        if os.path.exists("query_structures_DAT"):
            dat_files_in_query = [f for f in os.listdir("query_structures_DAT") if f.endswith('.dat')]
        
        if not dat_files_in_query:
            print(f"❌ 没有找到查询结构DAT文件: {json_filename}")
            return False
        
        # 使用最新生成的DAT文件作为查询结构
        query_name = os.path.splitext(dat_files_in_query[-1])[0]  # 使用最后一个文件
        print(f"使用查询结构: {query_name}")
        
        # 检查参考DAT目录
        dat_dir = "DAT"
        if not os.path.exists(dat_dir):
            print(f"❌ DAT目录不存在: {dat_dir}")
            return False
        
        # 获取所有参考DAT文件
        ref_dat_files = [f for f in os.listdir(dat_dir) if f.endswith('.dat')]
        if not ref_dat_files:
            print(f"❌ 参考DAT目录中没有找到.dat文件")
            return False
        
        print(f"找到 {len(ref_dat_files)} 个参考结构")
        
        # 对每个参考DAT文件运行DALI对比
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
            
            # 移动结果文件到序列专用目录
            default_output = f"{query_name}.txt"
            target_output = f"{sequence_dali_dir}/{query_name}_vs_{ref_id}.txt"
            
            if os.path.exists(default_output):
                try:
                    os.rename(default_output, target_output)
                    successful_comparisons += 1
                except Exception as e:
                    print(f"⚠️ 移动输出文件失败: {e}")
        
        print(f"✅ DALI对比完成 {json_filename}: {successful_comparisons} 个成功对比")
        return successful_comparisons > 0
    
    def extract_single_results(self, json_filename):
        """提取单个序列的Z-score结果"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"📊 提取结果: {json_filename}")
        
        sequence_dali_dir = f"{self.results_dir}/{output_name}_dali_results"
        
        if not os.path.exists(sequence_dali_dir):
            print(f"❌ DALI结果目录不存在: {json_filename}")
            return False
        
        # 修改的extract_zscores.py逻辑，专门处理这个序列
        zscores = {}
        
        def parse_z_score(file_path):
            with open(file_path) as f:
                for line in f:
                    if line.strip().startswith("1:"):
                        parts = line.strip().split()
                        try:
                            z = float(parts[2])  # 第三个字段是 Z-score
                            return z
                        except:
                            return None
            return None
        
        # 遍历该序列的结果目录
        for fname in os.listdir(sequence_dali_dir):
            if fname.endswith(".txt"):
                fpath = os.path.join(sequence_dali_dir, fname)
                z = parse_z_score(fpath)
                if z is not None:
                    zscores[fname] = z
        
        # 排序输出
        sorted_results = sorted(zscores.items(), key=lambda x: -x[1])
        
        # 保存到序列专用的CSV文件
        csv_file = f"{self.results_dir}/{output_name}_zscores.csv"
        
        if sorted_results:
            print(f"🎯 找到 {len(sorted_results)} 个Z-score结果")
            
            best = sorted_results[0]
            print(f"🏆 最佳匹配: {best[0]} with Z = {best[1]:.2f}")
            
            # 保存为 CSV
            with open(csv_file, "w") as f:
                f.write("filename,Z-score\n")
                for fname, z in sorted_results:
                    f.write(f"{fname},{z:.2f}\n")
            
            print(f"✅ Z-scores保存到: {csv_file}")
            return True
        else:
            print(f"❌ 没有找到有效的Z-score结果: {json_filename}")
            # 创建空的CSV文件表示处理完成但无结果
            with open(csv_file, "w") as f:
                f.write("filename,Z-score\n")
                f.write("No matches found,0.0\n")
            return False
    
    def cleanup_intermediate_files(self, json_filename):
        """清理中间文件以节省空间"""
        output_name = os.path.splitext(json_filename)[0]
        
        # 删除AlphaFold输出目录（保留PDB文件）
        alphafold_output_dir = f"{self.base_dir}/{output_name}_output"
        if os.path.exists(alphafold_output_dir):
            try:
                subprocess.run(["rm", "-rf", alphafold_output_dir], check=True)
                print(f"🗑️ 清理AlphaFold输出目录: {output_name}")
            except:
                print(f"⚠️ 清理失败: {alphafold_output_dir}")
        
        # 清理query_structures_DAT中的临时文件
        if os.path.exists("query_structures_DAT"):
            for dat_file in glob.glob("query_structures_DAT/*.dat"):
                try:
                    os.remove(dat_file)
                except:
                    pass
    
    def process_single_sequence(self, json_filename):
        """处理单个序列的完整流程"""
        output_name = os.path.splitext(json_filename)[0]
        print(f"\n{'='*60}")
        print(f"🎯 处理序列: {json_filename}")
        print(f"📋 输出名称: {output_name}")
        print(f"{'='*60}")
        
        try:
            # 步骤1: AlphaFold预测
            if not self.run_single_alphafold_prediction(json_filename):
                return False
            
            # 步骤2: 转换CIF到PDB
            if not self.convert_single_cif_to_pdb(json_filename):
                return False
            
            # 步骤3: 转换PDB到DAT
            if not self.convert_single_pdb_to_dat(json_filename):
                return False
            
            # 步骤4: DALI对比
            if not self.run_single_dali_comparison(json_filename):
                return False
            
            # 步骤5: 提取结果
            if not self.extract_single_results(json_filename):
                print(f"⚠️ Z-score提取失败，但继续处理")
            
            # 步骤6: 清理中间文件
            self.cleanup_intermediate_files(json_filename)
            
            print(f"✅ 序列处理完成: {json_filename}")
            return True
            
        except Exception as e:
            print(f"❌ 序列处理失败 {json_filename}: {str(e)}")
            return False
    
    def run_batch_pipeline(self):
        """运行批量处理流水线"""
        print(f"🚀 开始批量远程处理流水线...")
        print(f"📂 工作目录: {self.base_dir}")
        print(f"📁 输入目录: {self.input_dir}")
        print(f"📁 结果目录: {self.results_dir}")
        print("="*60)
        
        # 获取所有JSON文件
        json_files = self.get_json_files()
        if not json_files:
            print("❌ 没有找到JSON文件")
            return False
        
        # 统计变量
        total_files = len(json_files)
        successful_processes = 0
        failed_processes = 0
        
        print(f"📊 开始处理 {total_files} 个序列...")
        
        # 逐个处理序列
        for i, json_filename in enumerate(json_files, 1):
            print(f"\n[{i}/{total_files}] 当前处理: {json_filename}")
            
            start_time = time.time()
            
            if self.process_single_sequence(json_filename):
                successful_processes += 1
                elapsed = time.time() - start_time
                print(f"✅ 序列完成 ({elapsed:.1f}s): {json_filename}")
            else:
                failed_processes += 1
                elapsed = time.time() - start_time
                print(f"❌ 序列失败 ({elapsed:.1f}s): {json_filename}")
            
            # 显示进度
            remaining = total_files - i
            if remaining > 0:
                print(f"📈 进度: {i}/{total_files} 完成, 剩余 {remaining} 个")
        
        # 最终统计
        print("\n" + "="*60)
        print("📊 批量处理统计:")
        print(f"   总序列数: {total_files}")
        print(f"   成功处理: {successful_processes}")
        print(f"   失败数量: {failed_processes}")
        print(f"   成功率: {successful_processes/total_files*100:.1f}%")
        
        # 显示结果文件
        csv_files = glob.glob(f"{self.results_dir}/*_zscores.csv")
        print(f"\n📋 生成的结果文件: {len(csv_files)} 个")
        print(f"📁 结果目录: {self.results_dir}")
        
        if csv_files:
            print("🎉 批量处理完成!")
            print(f"📝 后续步骤:")
            print(f"   1. 检查结果: ls -la {self.results_dir}/")
            print(f"   2. 下载结果: scp -r webserver@coulomb.phys.ucl.ac.uk:{self.results_dir}/ ./")
            return True
        else:
            print("❌ 没有生成任何结果文件")
            return False

def main():
    parser = argparse.ArgumentParser(
        description='批量远程AlphaFold预测和结构对比流水线',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python batch_remote_pipeline.py                 # 批量处理所有JSON文件
  python batch_remote_pipeline.py --check         # 检查环境
  python batch_remote_pipeline.py --results-only  # 只提取已有结果
        """
    )
    
    parser.add_argument('--check', '-c', action='store_true', 
                       help='只检查环境，不运行流水线')
    parser.add_argument('--results-only', '-r', action='store_true', 
                       help='只提取结果，跳过预测步骤')
    
    args = parser.parse_args()
    
    pipeline = BatchRemotePipeline()
    
    if args.check:
        print("🔍 环境检查模式")
        # 基本环境检查
        required_scripts = [
            "convert_cif_to_pdb.py",
            "import_pdbs_to_dat.py"
        ]
        
        all_good = True
        for script in required_scripts:
            if os.path.exists(script):
                print(f"✅ 找到脚本: {script}")
            else:
                print(f"❌ 缺少脚本: {script}")
                all_good = False
        
        if os.path.exists(pipeline.input_dir):
            json_count = len(glob.glob(os.path.join(pipeline.input_dir, "*.json")))
            print(f"✅ 输入目录存在，包含 {json_count} 个JSON文件")
        else:
            print(f"❌ 输入目录不存在: {pipeline.input_dir}")
            all_good = False
        
        if all_good:
            print("✅ 环境检查通过，可以运行批量流水线")
        else:
            print("❌ 环境检查失败")
        return
    
    if args.results_only:
        print("📊 仅提取结果模式")
        # 这里可以添加只提取结果的逻辑
        print("功能开发中...")
        return
    
    # 运行批量流水线
    success = pipeline.run_batch_pipeline()
    
    if success:
        print("\n🎊 批量流水线执行成功!")
        sys.exit(0)
    else:
        print("\n💥 批量流水线执行失败!")
        sys.exit(1)

if __name__ == "__main__":
    main()
    