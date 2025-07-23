#!/usr/bin/env python3
"""
DALI序列提取器 - 从Excel文件提取DALI选中的序列并转换为JSON
固定路径版本 - 直接运行即可
"""

import json
import sys
import os
import pandas as pd
from typing import Dict, List, Optional
from pathlib import Path

class DALISequenceExtractor:
    def __init__(self):
        # 固定的绝对路径
        self.excel_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/NTP_Analysis_Report_20250723_133848.xlsx"
        self.fasta_base_path = Path("/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files")
        
    def validate_sequence(self, sequence: str) -> bool:
        """验证氨基酸序列是否有效 - 包含扩展IUPAC代码"""
        # 标准20种氨基酸 + 扩展IUPAC代码
        valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWYXBZJUO*')
        # A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y - 标准20种
        # X - 未知氨基酸
        # B - Asn或Asp
        # Z - Gln或Glu  
        # J - Leu或Ile
        # U - 硒代半胱氨酸
        # O - 吡咯赖氨酸
        # * - 终止密码子
        cleaned_sequence = sequence.replace(' ', '').replace('\n', '').upper()
        return all(aa in valid_amino_acids for aa in cleaned_sequence)
    
    def clean_sequence(self, sequence: str) -> str:
        """清理序列字符串"""
        return sequence.replace(' ', '').replace('\n', '').replace('\r', '').upper().strip()
    
    def read_dali_records(self) -> List[Dict]:
        """从Excel文件读取DALI记录"""
        try:
            print(f"读取Excel文件: {self.excel_path}")
            df = pd.read_excel(self.excel_path, sheet_name='Detailed_Results')
            
            # 筛选Result_Category为DALI的记录
            dali_records = df[df['Result_Category'] == 'DALI']
            print(f"找到 {len(dali_records)} 条DALI记录")
            
            records = []
            for _, row in dali_records.iterrows():
                records.append({
                    'folder': str(row['Folder']),
                    'sequence_id': str(row['Sequence_ID']),
                    'html_file': str(row['HTML_File']),
                    'found_superfamily': str(row['Found_Superfamily']) if pd.notna(row['Found_Superfamily']) else 'N/A'
                })
            
            return records
            
        except FileNotFoundError:
            print(f"错误: Excel文件不存在 - {self.excel_path}")
            sys.exit(1)
        except Exception as e:
            print(f"读取Excel文件时出错: {e}")
            sys.exit(1)
    
    def find_fasta_file(self, folder: str, sequence_id: str) -> Optional[Path]:
        """查找对应的FASTA文件"""
        folder_path = self.fasta_base_path / folder
        
        if not folder_path.exists():
            print(f"警告: 文件夹不存在 - {folder_path}")
            return None
        
        # 可能的FASTA文件扩展名
        extensions = ['.fasta', '.fa', '.fas', '.seq', '.txt']
        
        # 精确匹配
        for ext in extensions:
            fasta_file = folder_path / f"{sequence_id}{ext}"
            if fasta_file.exists():
                return fasta_file
        
        # 模糊匹配
        for file in folder_path.iterdir():
            if file.is_file() and sequence_id in file.stem:
                return file
        
        print(f"警告: 未找到FASTA文件 - {folder}/{sequence_id}")
        return None
    
    def read_fasta_sequence(self, fasta_file: Path) -> Optional[str]:
        """从FASTA文件读取序列"""
        try:
            with open(fasta_file, 'r', encoding='utf-8') as f:
                content = f.read().strip()
            
            # 处理FASTA格式 - 跳过以>开头的行
            lines = content.split('\n')
            sequence = ""
            
            for line in lines:
                line = line.strip()
                if not line.startswith('>') and line:
                    sequence += line
            
            return self.clean_sequence(sequence) if sequence else None
                
        except Exception as e:
            print(f"读取FASTA文件出错 {fasta_file}: {e}")
            return None
    
    def extract_dali_sequences(self) -> List[Dict]:
        """提取所有DALI选中的序列"""
        dali_records = self.read_dali_records()
        
        if not dali_records:
            print("没有找到DALI记录")
            return []
        
        results = []
        success_count = 0
        
        print(f"\n开始处理 {len(dali_records)} 个DALI记录...")
        
        for i, record in enumerate(dali_records, 1):
            folder = record['folder']
            sequence_id = record['sequence_id']
            found_superfamily = record['found_superfamily']
            
            print(f"[{i}/{len(dali_records)}] 处理: {sequence_id} (文件夹: {folder})")
            
            # 查找FASTA文件
            fasta_file = self.find_fasta_file(folder, sequence_id)
            if not fasta_file:
                continue
            
            # 读取序列
            sequence = self.read_fasta_sequence(fasta_file)
            if not sequence:
                print(f"  -> 序列为空，跳过")
                continue
            
            # 验证序列
            if not self.validate_sequence(sequence):
                print(f"  -> 序列包含无效字符，跳过")
                continue
            
            # 创建JSON条目 - 简单格式
            json_entry = {
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": sequence,
                            "count": 1
                        }
                    }
                ],
                "name": sequence_id  # 直接使用序列ID作为名称
            }
            
            results.append(json_entry)
            success_count += 1
            print(f"  -> 成功 (长度: {len(sequence)})")
        
        print(f"\n处理完成: {success_count}/{len(dali_records)} 个序列成功")
        return results
    
    def save_sequences(self, sequences: List[Dict], output_path: str, single_file: bool = False):
        """保存序列到JSON文件"""
        if not sequences:
            print("没有序列需要保存")
            return
        
        if single_file:
            # 保存为单个文件
            output_file = Path(output_path)
            if not output_file.suffix:
                output_file = output_file.with_suffix('.json')
            
            try:
                with open(output_file, 'w', encoding='utf-8') as f:
                    json.dump(sequences, f, indent=2, ensure_ascii=False)
                print(f"保存到单个文件: {output_file} ({len(sequences)} 个序列)")
            except Exception as e:
                print(f"保存文件出错: {e}")
        else:
            # 保存为多个文件到文件夹
            output_dir = Path(output_path)
            output_dir.mkdir(exist_ok=True)
            
            saved_count = 0
            for sequence in sequences:
                filename = f"{sequence['name']}.json"  # 使用序列名作为文件名
                filepath = output_dir / filename
                
                try:
                    with open(filepath, 'w', encoding='utf-8') as f:
                        json.dump([sequence], f, indent=2, ensure_ascii=False)  # 保持数组格式
                    saved_count += 1
                except Exception as e:
                    print(f"保存文件出错 {filepath}: {e}")
            
            print(f"保存到文件夹: {output_dir} ({saved_count} 个文件)")

def main():
    # 固定的文件路径
    excel_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/NTP_Analysis_Report_20250723_133848.xlsx"
    fasta_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files"
    output_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/dali_sequences"
    
    # 检查文件是否存在
    if not os.path.exists(excel_path):
        print(f"错误: Excel文件不存在 - {excel_path}")
        sys.exit(1)
    
    if not os.path.exists(fasta_path):
        print(f"错误: FASTA基础路径不存在 - {fasta_path}")
        sys.exit(1)
    
    print("=" * 60)
    print("DALI序列提取器")
    print("=" * 60)
    print(f"Excel文件: {excel_path}")
    print(f"FASTA路径: {fasta_path}")
    print(f"输出路径: {output_path}")
    print("=" * 60)
    
    # 创建提取器并处理
    extractor = DALISequenceExtractor()
    sequences = extractor.extract_dali_sequences()
    
    if sequences:
        extractor.save_sequences(sequences, output_path, single_file=False)
        print(f"\n✅ 处理完成! 成功提取 {len(sequences)} 个DALI序列")
        print(f"JSON文件已保存到: {output_path}")
    else:
        print("\n❌ 没有成功提取任何序列")

if __name__ == "__main__":
    main()