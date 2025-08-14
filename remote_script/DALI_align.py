#!/usr/bin/env python3
"""
Use DALI's applymatrix.pl to transform PDB coordinates using DALI alignment results
Usage: python dali_txt_to_pdb.py <original_pdb> <dali_txt_file> [--output output.pdb]
"""

import os
import sys
import subprocess
import argparse
import glob
from pathlib import Path

class DaliMatrixApplier:
    def __init__(self, original_pdb, dali_txt_file, output_pdb=None):
        self.original_pdb = original_pdb
        self.dali_txt_file = dali_txt_file
        
        if output_pdb:
            self.output_pdb = output_pdb
        else:
            # 自动生成输出文件名
            base_name = os.path.splitext(os.path.basename(original_pdb))[0]
            txt_base = os.path.splitext(os.path.basename(dali_txt_file))[0]
            self.output_pdb = f"{base_name}_transformed_by_{txt_base}.pdb"
        
        # DALI applymatrix.pl 工具路径
        self.applymatrix_path = "/home/wenhao/6tx0/software/dali/DaliLite.v5/bin/applymatrix.pl"
        
        print(f"🔧 Initialized DaliMatrixApplier:")
        print(f"   Original PDB: {self.original_pdb}")
        print(f"   DALI TXT file: {self.dali_txt_file}")
        print(f"   Output PDB: {self.output_pdb}")
        print(f"   applymatrix.pl: {self.applymatrix_path}")
    
    def check_prerequisites(self):
        """检查必需文件是否存在"""
        print(f"🔍 Checking prerequisites...")
        
        # 检查原始PDB文件
        if not os.path.exists(self.original_pdb):
            print(f"❌ Original PDB file not found: {self.original_pdb}")
            return False
        
        # 检查DALI TXT文件
        if not os.path.exists(self.dali_txt_file):
            print(f"❌ DALI TXT file not found: {self.dali_txt_file}")
            return False
        
        # 检查applymatrix.pl工具
        if not os.path.exists(self.applymatrix_path):
            print(f"❌ applymatrix.pl not found: {self.applymatrix_path}")
            return False
        
        print(f"✅ All prerequisites satisfied")
        return True
    
    def analyze_dali_txt(self):
        """分析DALI TXT文件内容"""
        print(f"📊 Analyzing DALI TXT file...")
        
        try:
            with open(self.dali_txt_file, 'r') as f:
                content = f.read()
            
            print(f"📄 DALI TXT file content preview:")
            lines = content.strip().split('\n')
            for i, line in enumerate(lines[:20]):  # 显示前20行
                print(f"   {i+1:2d}: {line}")
            
            if len(lines) > 20:
                print(f"   ... (total {len(lines)} lines)")
            
            # 查找特定的模式
            if "No:  Chain   Z    rmsd lali nres  %id PDB" in content:
                print(f"✅ Found alignment summary table")
            
            if "1:" in content:
                print(f"✅ Found alignment entry")
            
            return True
            
        except Exception as e:
            print(f"❌ Failed to read DALI TXT file: {e}")
            return False
    
    def apply_transformation(self):
        """使用DALI的applymatrix.pl应用变换"""
        print(f"🔄 Applying DALI transformation...")
        
        # 构建命令
        cmd = [
            self.applymatrix_path,
            self.original_pdb
        ]
        
        print(f"🖥️  Executing command:")
        print(f"   {' '.join(cmd)} < {self.dali_txt_file} > {self.output_pdb}")
        
        try:
            # 读取DALI TXT文件内容
            with open(self.dali_txt_file, 'r') as txt_file:
                txt_content = txt_file.read()
            
            # 执行applymatrix.pl命令
            process = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            stdout, stderr = process.communicate(input=txt_content)
            
            if process.returncode == 0:
                # 写入输出文件
                with open(self.output_pdb, 'w') as output_file:
                    output_file.write(stdout)
                
                print(f"✅ Transformation completed successfully")
                print(f"📁 Output saved to: {self.output_pdb}")
                return True
            else:
                print(f"❌ applymatrix.pl failed with return code {process.returncode}")
                print(f"STDERR: {stderr}")
                return False
                
        except Exception as e:
            print(f"❌ Failed to execute applymatrix.pl: {e}")
            return False
    
    def analyze_output(self):
        """分析输出的PDB文件"""
        print(f"📊 Analyzing output PDB file...")
        
        if not os.path.exists(self.output_pdb):
            print(f"❌ Output PDB file not found: {self.output_pdb}")
            return False
        
        try:
            atom_count = 0
            chains = set()
            
            with open(self.output_pdb, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        atom_count += 1
                        chain = line[21].strip()
                        if chain:
                            chains.add(chain)
            
            print(f"📈 Output statistics:")
            print(f"   Total atoms: {atom_count}")
            print(f"   Chains: {sorted(chains) if chains else 'None'}")
            
            # 比较原始文件
            original_atom_count = 0
            with open(self.original_pdb, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        original_atom_count += 1
            
            print(f"   Original atoms: {original_atom_count}")
            
            if atom_count == original_atom_count:
                print(f"✅ Atom count matches original file")
            else:
                print(f"⚠️  Atom count differs from original")
            
            return True
            
        except Exception as e:
            print(f"❌ Failed to analyze output: {e}")
            return False
    
    def generate_summary(self):
        """生成变换摘要报告"""
        print(f"📋 Generating transformation summary...")
        
        summary_file = f"{os.path.splitext(self.output_pdb)[0]}_summary.txt"
        
        try:
            with open(summary_file, 'w') as f:
                f.write(f"DALI Matrix Transformation Summary\n")
                f.write(f"=" * 40 + "\n\n")
                f.write(f"Original PDB file: {self.original_pdb}\n")
                f.write(f"DALI TXT file: {self.dali_txt_file}\n")
                f.write(f"Transformed PDB file: {self.output_pdb}\n")
                f.write(f"Transformation date: {os.popen('date').read().strip()}\n\n")
                
                # 添加DALI TXT文件的前几行
                f.write(f"DALI TXT file content:\n")
                f.write(f"-" * 20 + "\n")
                with open(self.dali_txt_file, 'r') as txt_f:
                    lines = txt_f.readlines()[:20]
                    for line in lines:
                        f.write(line)
                if len(lines) >= 20:
                    f.write("...\n")
                
                f.write(f"\n" + "-" * 20 + "\n")
                f.write(f"Command used:\n")
                f.write(f"{self.applymatrix_path} {self.original_pdb} < {self.dali_txt_file} > {self.output_pdb}\n")
            
            print(f"📄 Summary saved to: {summary_file}")
            return summary_file
            
        except Exception as e:
            print(f"❌ Failed to generate summary: {e}")
            return None
    
    def run_transformation(self):
        """运行完整的变换流程"""
        print(f"🚀 Starting DALI matrix transformation...")
        print("=" * 60)
        
        try:
            # 检查先决条件
            if not self.check_prerequisites():
                return False
            
            # 分析DALI TXT文件
            if not self.analyze_dali_txt():
                return False
            
            # 应用变换
            if not self.apply_transformation():
                return False
            
            # 分析输出
            self.analyze_output()
            
            # 生成摘要
            summary_file = self.generate_summary()
            
            print("=" * 60)
            print("🎉 Transformation completed successfully!")
            print(f"📁 Files generated:")
            print(f"   - Transformed PDB: {self.output_pdb}")
            if summary_file:
                print(f"   - Summary report: {summary_file}")
            
            return True
            
        except Exception as e:
            print(f"❌ Transformation failed: {str(e)}")
            return False

def find_dali_txt_files(directory, pattern="*_vs_*.txt"):
    """查找DALI结果TXT文件"""
    search_path = os.path.join(directory, pattern)
    txt_files = glob.glob(search_path)
    return sorted(txt_files)

def main():
    parser = argparse.ArgumentParser(
        description='Apply DALI transformation matrix to PDB coordinates using applymatrix.pl',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 基本用法
  python dali_txt_to_pdb.py original.pdb a0l1A_vs_7tgkD.txt
  
  # 指定输出文件
  python dali_txt_to_pdb.py original.pdb alignment_result.txt --output transformed.pdb
  
  # 查找和处理多个TXT文件
  python dali_txt_to_pdb.py original.pdb --search-dir ./dali_results --pattern "*_vs_7tgkD.txt"
  
  # 批量处理目录中的所有TXT文件
  python dali_txt_to_pdb.py original.pdb --batch-dir ./dali_results
        """
    )
    
    parser.add_argument('original_pdb', help='Original PDB file to transform')
    parser.add_argument('dali_txt_file', nargs='?', help='DALI alignment result TXT file')
    parser.add_argument('--output', '-o', help='Output PDB file name')
    parser.add_argument('--search-dir', help='Directory to search for DALI TXT files')
    parser.add_argument('--pattern', default='*_vs_*.txt', help='Pattern to match TXT files (default: *_vs_*.txt)')
    parser.add_argument('--batch-dir', help='Process all TXT files in directory')
    parser.add_argument('--target', help='Specific target to look for (e.g., 7tgkD)')
    
    args = parser.parse_args()
    
    # 检查原始PDB文件
    if not os.path.exists(args.original_pdb):
        print(f"❌ Original PDB file not found: {args.original_pdb}")
        sys.exit(1)
    
    # 处理不同的使用模式
    txt_files_to_process = []
    
    if args.batch_dir:
        # 批量处理模式
        print(f"🔍 Searching for TXT files in: {args.batch_dir}")
        txt_files = find_dali_txt_files(args.batch_dir, args.pattern)
        
        if args.target:
            # 过滤特定目标
            txt_files = [f for f in txt_files if args.target in os.path.basename(f)]
            print(f"🎯 Filtering for target: {args.target}")
        
        txt_files_to_process = txt_files
        
    elif args.search_dir:
        # 搜索模式
        print(f"🔍 Searching for TXT files in: {args.search_dir}")
        txt_files = find_dali_txt_files(args.search_dir, args.pattern)
        
        if not txt_files:
            print(f"❌ No TXT files found matching pattern: {args.pattern}")
            sys.exit(1)
        
        print(f"📄 Found {len(txt_files)} TXT files:")
        for i, txt_file in enumerate(txt_files):
            print(f"   {i+1}. {os.path.basename(txt_file)}")
        
        if args.target:
            # 查找特定目标
            target_files = [f for f in txt_files if args.target in os.path.basename(f)]
            if target_files:
                txt_files_to_process = target_files
                print(f"🎯 Found target files: {[os.path.basename(f) for f in target_files]}")
            else:
                print(f"❌ No files found containing: {args.target}")
                sys.exit(1)
        else:
            # 处理第一个文件
            txt_files_to_process = [txt_files[0]]
            print(f"📝 Processing first file: {os.path.basename(txt_files[0])}")
            
    elif args.dali_txt_file:
        # 单文件模式
        if not os.path.exists(args.dali_txt_file):
            print(f"❌ DALI TXT file not found: {args.dali_txt_file}")
            sys.exit(1)
        txt_files_to_process = [args.dali_txt_file]
        
    else:
        print(f"❌ Must provide either dali_txt_file, --search-dir, or --batch-dir")
        sys.exit(1)
    
    # 处理文件
    success_count = 0
    total_count = len(txt_files_to_process)
    
    for i, txt_file in enumerate(txt_files_to_process):
        print(f"\n{'='*60}")
        print(f"📁 Processing {i+1}/{total_count}: {os.path.basename(txt_file)}")
        print(f"{'='*60}")
        
        # 为批量处理生成输出文件名
        if total_count > 1:
            base_name = os.path.splitext(os.path.basename(args.original_pdb))[0]
            txt_base = os.path.splitext(os.path.basename(txt_file))[0]
            output_pdb = f"{base_name}_transformed_by_{txt_base}.pdb"
        else:
            output_pdb = args.output
        
        # 创建变换器并运行
        transformer = DaliMatrixApplier(args.original_pdb, txt_file, output_pdb)
        
        if transformer.run_transformation():
            success_count += 1
        else:
            print(f"❌ Failed to process: {txt_file}")
    
    # 最终报告
    print(f"\n{'='*60}")
    print(f"🎊 Batch processing completed!")
    print(f"📊 Successfully processed: {success_count}/{total_count} files")
    
    if success_count == total_count:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()