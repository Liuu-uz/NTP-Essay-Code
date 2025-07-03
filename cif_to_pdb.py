#!/usr/bin/env python3
"""
CIF to PDB Converter Script
将CIF格式的结构文件转换为PDB格式

依赖库: biopython
安装: pip install biopython
"""

import os
import sys
import argparse
from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio.PDB.PDBIO import Select


class CifToPdbConverter:
    def __init__(self):
        self.cif_parser = MMCIFParser()
        self.pdb_io = PDBIO()
    
    def convert_file(self, cif_file, pdb_file=None, structure_id=None):
        """
        转换单个CIF文件为PDB格式
        
        Args:
            cif_file (str): 输入的CIF文件路径
            pdb_file (str): 输出的PDB文件路径（可选）
            structure_id (str): 结构ID（可选）
        """
        try:
            # 检查输入文件是否存在
            if not os.path.exists(cif_file):
                raise FileNotFoundError(f"CIF文件不存在: {cif_file}")
            
            # 生成输出文件名
            if pdb_file is None:
                base_name = os.path.splitext(os.path.basename(cif_file))[0]
                pdb_file = f"{base_name}.pdb"
            
            # 生成结构ID
            if structure_id is None:
                structure_id = os.path.splitext(os.path.basename(cif_file))[0]
            
            print(f"正在转换: {cif_file} -> {pdb_file}")
            
            # 解析CIF文件
            structure = self.cif_parser.get_structure(structure_id, cif_file)
            
            # 设置结构并保存为PDB
            self.pdb_io.set_structure(structure)
            self.pdb_io.save(pdb_file)
            
            print(f"转换成功: {pdb_file}")
            return True
            
        except Exception as e:
            print(f"转换失败: {e}")
            return False
    
    def convert_batch(self, input_dir, output_dir=None):
        """
        批量转换目录中的所有CIF文件
        
        Args:
            input_dir (str): 包含CIF文件的输入目录
            output_dir (str): 输出目录（可选）
        """
        if not os.path.exists(input_dir):
            raise FileNotFoundError(f"输入目录不存在: {input_dir}")
        
        if output_dir is None:
            output_dir = input_dir
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 查找所有CIF文件
        cif_files = [f for f in os.listdir(input_dir) 
                    if f.lower().endswith(('.cif', '.mmcif'))]
        
        if not cif_files:
            print(f"在目录 {input_dir} 中没有找到CIF文件")
            return
        
        print(f"找到 {len(cif_files)} 个CIF文件")
        
        success_count = 0
        for cif_file in cif_files:
            input_path = os.path.join(input_dir, cif_file)
            base_name = os.path.splitext(cif_file)[0]
            output_path = os.path.join(output_dir, f"{base_name}.pdb")
            
            if self.convert_file(input_path, output_path):
                success_count += 1
        
        print(f"\n批量转换完成: {success_count}/{len(cif_files)} 个文件转换成功")


def main():
    parser = argparse.ArgumentParser(description='将CIF格式文件转换为PDB格式')
    parser.add_argument('input', help='输入CIF文件或包含CIF文件的目录')
    parser.add_argument('-o', '--output', help='输出PDB文件路径或输出目录')
    parser.add_argument('-b', '--batch', action='store_true', 
                       help='批量转换模式（当输入为目录时）')
    parser.add_argument('--structure-id', help='结构ID（仅用于单文件转换）')
    
    args = parser.parse_args()
    
    converter = CifToPdbConverter()
    
    try:
        if os.path.isfile(args.input):
            # 单文件转换
            converter.convert_file(args.input, args.output, args.structure_id)
        elif os.path.isdir(args.input):
            # 目录批量转换
            converter.convert_batch(args.input, args.output)
        else:
            print(f"错误: 输入路径不存在: {args.input}")
            sys.exit(1)
            
    except Exception as e:
        print(f"程序执行错误: {e}")
        sys.exit(1)


if __name__ == "__main__":
    # 如果直接运行脚本，可以在这里添加示例用法
    if len(sys.argv) == 1:
        print("CIF to PDB Converter")
        print("\n使用方法:")
        print("  python cif_to_pdb.py input.cif                    # 转换单个文件")
        print("  python cif_to_pdb.py input.cif -o output.pdb      # 指定输出文件")
        print("  python cif_to_pdb.py /path/to/cif/files/          # 批量转换目录")
        print("  python cif_to_pdb.py /path/to/cif/ -o /output/    # 指定输出目录")
        print("\n简单示例:")
        
        # 创建转换器实例供交互使用
        converter = CifToPdbConverter()
        
        # 示例：如果用户想要交互式使用
        print("\n如果您想要交互式使用，请取消注释下面的代码:")
        print("# input_file = input('请输入CIF文件路径: ')")
        print("# converter.convert_file(input_file)")
        
    else:
        main()
        