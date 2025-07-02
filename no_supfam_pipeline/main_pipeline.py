#!/usr/bin/env python3
"""
蛋白质结构预测和分析主流程脚本
Date: 2025-07-02

功能：
1. 上传JSON文件到远程服务器
2. 远程执行AlphaFold预测
3. 转换CIF到PDB格式
4. 执行Dali结构比对

用法：
python main_pipeline.py <json_file_path> [--protein_name <name>]
"""

import os
import sys
import argparse
import subprocess
import time
from pathlib import Path

class ProteinPredictionPipeline:
    def __init__(self, json_file_path, protein_name=None):
        self.json_file_path = Path(json_file_path)
        self.protein_name = protein_name or self.json_file_path.stem
        self.remote_host = "webserver@coulomb.phys.ucl.ac.uk"
        self.remote_base_dir = "/home/webserver/student/zhijing"
        self.input_dir = f"{self.remote_base_dir}/input_jsons"
        self.output_dir = f"{self.remote_base_dir}/predicted_structures"
        
    def upload_json_file(self):
        """上传JSON文件到远程服务器"""
        print(f"正在上传 {self.json_file_path} 到远程服务器...")
        
        # 确保远程目录存在
        mkdir_cmd = f'ssh {self.remote_host} "mkdir -p {self.input_dir}"'
        result = subprocess.run(mkdir_cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"创建远程目录失败: {result.stderr}")
            return False
            
        # 上传文件
        scp_cmd = f"scp {self.json_file_path} {self.remote_host}:{self.input_dir}/{self.protein_name}.json"
        result = subprocess.run(scp_cmd, shell=True)
        
        if result.returncode == 0:
            print("文件上传成功!")
            return True
        else:
            print("文件上传失败!")
            return False
    
    def run_alphafold_prediction(self):
        """在远程服务器上运行AlphaFold预测"""
        print(f"正在为 {self.protein_name} 运行AlphaFold预测...")
        
        # 构建远程命令
        remote_cmd = f"""
        cd {self.remote_base_dir} && \
        conda activate protoneix && \
        protenix predict \
        --input {self.input_dir}/{self.protein_name}.json \
        --out_dir {self.output_dir}/{self.protein_name}_output \
        --seeds 101 \
        --use_msa_server
        """
        
        # 执行远程命令
        ssh_cmd = f'ssh {self.remote_host} "{remote_cmd}"'
        print(f"执行命令: {ssh_cmd}")
        
        result = subprocess.run(ssh_cmd, shell=True)
        
        if result.returncode == 0:
            print("AlphaFold预测完成!")
            return True
        else:
            print("AlphaFold预测失败!")
            return False
    
    def convert_cif_to_pdb(self):
        """转换CIF文件为PDB格式"""
        print(f"正在转换 {self.protein_name} 的CIF文件为PDB格式...")
        
        cif_path = f"{self.output_dir}/{self.protein_name}_output/{self.protein_name}/seed_101/predictions/{self.protein_name}_seed_101_sample_0.cif"
        pdb_path = f"{self.output_dir}/{self.protein_name}.pdb"
        
        remote_cmd = f"""
        cd {self.remote_base_dir} && \
        python3 convert_cif_to_pdb.py {cif_path} {pdb_path}
        """
        
        ssh_cmd = f'ssh {self.remote_host} "{remote_cmd}"'
        result = subprocess.run(ssh_cmd, shell=True)
        
        if result.returncode == 0:
            print("CIF到PDB转换完成!")
            return True
        else:
            print("CIF到PDB转换失败!")
            return False
    
    def run_dali_comparison(self):
        """运行Dali结构比对"""
        print(f"正在为 {self.protein_name} 运行Dali结构比对...")
        
        pdb_path = f"{self.output_dir}/{self.protein_name}.pdb"
        dali_db = f"{self.remote_base_dir}/dali_database"
        
        # 注意：Dali在base环境中运行，需要切换环境
        remote_cmd = f"""
        cd {self.remote_base_dir} && \
        conda deactivate && \
        conda activate base && \
        /home/wenhao/6tx0/software/dali/DaliLite.v5/bin/dali.pl \
        -query {pdb_path} \
        -db {dali_db} \
        -title {self.protein_name} \
        -outfmt summary
        """
        
        ssh_cmd = f'ssh {self.remote_host} "{remote_cmd}"'
        result = subprocess.run(ssh_cmd, shell=True)
        
        if result.returncode == 0:
            print("Dali结构比对完成!")
            return True
        else:
            print("Dali结构比对失败!")
            return False
    
    def download_results(self, local_output_dir="./results"):
        """下载结果文件到本地"""
        print("正在下载结果文件...")
        
        os.makedirs(local_output_dir, exist_ok=True)
        
        # 下载PDB文件
        pdb_remote = f"{self.remote_host}:{self.output_dir}/{self.protein_name}.pdb"
        subprocess.run(f"scp {pdb_remote} {local_output_dir}/", shell=True)
        
        # 下载Dali结果文件
        dali_remote = f"{self.remote_host}:{self.remote_base_dir}/{self.protein_name}.*"
        subprocess.run(f"scp {dali_remote} {local_output_dir}/ 2>/dev/null", shell=True)
        
        print(f"结果文件已下载到: {local_output_dir}")
    
    def run_full_pipeline(self):
        """运行完整流程"""
        print(f"开始处理蛋白质: {self.protein_name}")
        print("=" * 50)
        
        steps = [
            ("上传JSON文件", self.upload_json_file),
            ("运行AlphaFold预测", self.run_alphafold_prediction),
            ("转换CIF为PDB", self.convert_cif_to_pdb),
            ("运行Dali比对", self.run_dali_comparison),
            ("下载结果", self.download_results)
        ]
        
        for step_name, step_func in steps:
            print(f"\n步骤: {step_name}")
            print("-" * 30)
            
            if not step_func():
                print(f"错误: {step_name} 失败，停止流程")
                return False
            
            print(f"{step_name} 完成")
        
        print("\n" + "=" * 50)
        print("完整流程执行完成!")
        return True

def main():
    parser = argparse.ArgumentParser(description="蛋白质结构预测和分析流程")
    parser.add_argument("json_file", help="输入的JSON文件路径")
    parser.add_argument("--protein_name", help="蛋白质名称（默认使用文件名）")
    parser.add_argument("--download_only", action="store_true", help="仅下载结果文件")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.json_file):
        print(f"错误: 文件 {args.json_file} 不存在")
        sys.exit(1)
    
    pipeline = ProteinPredictionPipeline(args.json_file, args.protein_name)
    
    if args.download_only:
        pipeline.download_results()
    else:
        success = pipeline.run_full_pipeline()
        sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()