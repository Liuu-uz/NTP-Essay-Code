#!/usr/bin/env python3
"""
NTP Processing序列提取器 - 添加配体版本
解决不同文件夹中相同sequence_id的覆盖问题，并添加ATP配体信息
"""

import json
import sys
import os
import pandas as pd
from typing import Dict, List, Optional
from pathlib import Path
import re

class NTPProcessingSequenceExtractor:
    def __init__(self):
        # 固定的绝对路径
        self.excel_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/NTP_Analysis_Report.xlsx"
        self.fasta_base_path = Path("/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files")
        
    def sanitize_filename(self, filename: str) -> str:
        """清理文件名中的特殊字符"""
        # 替换文件系统不允许的字符为下划线
        sanitized = re.sub(r'[<>:"/\\|?*\x00-\x1f]', '_', filename)
        
        # 移除开头和结尾的点和空格
        sanitized = sanitized.strip('. ')
        
        # 确保不为空
        if not sanitized:
            sanitized = "unnamed"
        
        # 限制长度
        if len(sanitized) > 200:
            sanitized = sanitized[:200]
        
        return sanitized
        
    def validate_sequence(self, sequence: str) -> bool:
        """验证氨基酸序列是否有效 - 包含扩展IUPAC代码"""
        valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWYXBZJUO*')
        cleaned_sequence = sequence.replace(' ', '').replace('\n', '').upper()
        return all(aa in valid_amino_acids for aa in cleaned_sequence)
    
    def clean_sequence(self, sequence: str) -> str:
        """清理序列字符串"""
        return sequence.replace(' ', '').replace('\n', '').replace('\r', '').upper().strip()
    
    def read_ntp_processing_records(self) -> List[Dict]:
        """从Excel文件读取NTP processing SF found记录"""
        try:
            print(f"读取Excel文件: {self.excel_path}")
            
            try:
                df = pd.read_excel(self.excel_path, sheet_name='Detailed_Results')
            except ValueError:
                df = pd.read_excel(self.excel_path, sheet_name=0)
                print("使用第一个工作表")
            
            print(f"Excel文件总共有 {len(df)} 行数据")
            print(f"列名: {list(df.columns)}")
            
            # 筛选Result_Category为'NTP processing SF found'的记录
            ntp_records = df[df['Result_Category'] == 'NTP processing SF found']
            print(f"找到 {len(ntp_records)} 条'NTP processing SF found'记录")
            
            if len(ntp_records) == 0:
                unique_categories = df['Result_Category'].unique()
                print("可用的Result_Category值:")
                for cat in unique_categories:
                    count = len(df[df['Result_Category'] == cat])
                    print(f"  - '{cat}': {count} 条记录")
                
                print("\n尝试模糊匹配包含'NTP'的记录...")
                ntp_like_records = df[df['Result_Category'].str.contains('NTP', case=False, na=False)]
                print(f"找到 {len(ntp_like_records)} 条包含'NTP'的记录")
                
                if len(ntp_like_records) > 0:
                    ntp_records = ntp_like_records
                    print("使用包含'NTP'的记录")
            
            # 检查重复的sequence_id
            sequence_ids = ntp_records['Sequence_ID'].astype(str).tolist()
            folders = ntp_records['Folder'].astype(str).tolist()
            
            print(f"唯一sequence_id数: {len(set(sequence_ids))}")
            print(f"总记录数: {len(sequence_ids)}")
            
            # 统计重复情况
            from collections import Counter
            id_counts = Counter(sequence_ids)
            duplicates = {seq_id: count for seq_id, count in id_counts.items() if count > 1}
            
            if duplicates:
                print(f"发现 {len(duplicates)} 个重复的sequence_id:")
                for seq_id, count in list(duplicates.items())[:5]:  # 只显示前5个
                    print(f"  - {seq_id}: {count} 次")
                    # 显示这些重复ID来自哪些文件夹
                    related_folders = [folders[i] for i, sid in enumerate(sequence_ids) if sid == seq_id]
                    print(f"    来自文件夹: {set(related_folders)}")
                if len(duplicates) > 5:
                    print(f"  ... 还有 {len(duplicates) - 5} 个重复ID")
            
            records = []
            for _, row in ntp_records.iterrows():
                records.append({
                    'folder': str(row['Folder']),
                    'sequence_id': str(row['Sequence_ID']),
                    'html_file': str(row['HTML_File']),
                    'found_superfamily': str(row['Found_Superfamily']) if pd.notna(row['Found_Superfamily']) else 'N/A',
                    'result_category': str(row['Result_Category'])
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
    
    def create_unique_name(self, folder: str, sequence_id: str) -> str:
        """创建唯一的序列名称，结合序列ID和文件夹"""
        # 方案: sequenceID_folder
        unique_name = f"{sequence_id}_{folder}"
        
        # 清理特殊字符
        unique_name = self.sanitize_filename(unique_name)
        
        return unique_name
    
    def extract_ntp_processing_sequences(self) -> List[Dict]:
        """提取所有NTP processing SF found的序列，包含ATP配体信息"""
        ntp_records = self.read_ntp_processing_records()
        
        if not ntp_records:
            print("没有找到NTP processing SF found记录")
            return []
        
        results = []
        success_count = 0
        seen_combinations = set()  # 跟踪已处理的 (folder, sequence_id) 组合
        duplicate_combinations = 0
        
        print(f"\n开始处理 {len(ntp_records)} 个NTP processing记录...")
        
        for i, record in enumerate(ntp_records, 1):
            folder = record['folder']
            sequence_id = record['sequence_id']
            found_superfamily = record['found_superfamily']
            result_category = record['result_category']
            
            # 创建唯一标识符
            combination_key = (folder, sequence_id)
            
            print(f"[{i}/{len(ntp_records)}] 处理: {sequence_id} (文件夹: {folder}, 类别: {result_category})")
            
            # 检查是否已经处理过这个组合
            if combination_key in seen_combinations:
                duplicate_combinations += 1
                print(f"  -> 跳过重复组合: {folder}/{sequence_id}")
                continue
            
            seen_combinations.add(combination_key)
            
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
            
            # 创建唯一名称
            unique_name = self.create_unique_name(folder, sequence_id)
            
            # 创建JSON条目，包含配体信息
            json_entry = {
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": sequence,
                            "count": 1
                          }
                    },
                    {
                        "ligand": {
                            "ligand": "CCD_ATP",
                            "count": 1
                        }
                    }
                ],
                "name": unique_name,  # 使用唯一名称
            }
            
            results.append(json_entry)
            success_count += 1
            print(f"  -> 成功 (唯一名称: {unique_name}, 长度: {len(sequence)}, 超家族: {found_superfamily}, 配体: ATP)")
        
        print(f"\n处理完成:")
        print(f"- 总记录数: {len(ntp_records)}")
        print(f"- 成功处理: {success_count}")
        print(f"- 跳过的重复组合: {duplicate_combinations}")
        print(f"- 唯一的 (folder, sequence_id) 组合: {len(seen_combinations)}")
        print(f"- 所有序列都包含ATP配体信息")
        
        return results
    
    def save_sequences(self, sequences: List[Dict], output_path: str, single_file: bool = False):
        """保存序列到JSON文件 - 带详细跟踪"""
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
            failed_count = 0
            failed_files = []
            seen_filenames = set()
            
            print(f"开始保存 {len(sequences)} 个序列到 {output_dir}")
            
            for i, sequence in enumerate(sequences):
                unique_name = sequence['name']
                original_id = sequence.get('original_sequence_id', 'unknown')
                folder = sequence.get('folder', 'unknown')
                
                # 确保文件名唯一
                base_filename = self.sanitize_filename(unique_name)
                filename = f"{base_filename}.json"
                
                # 如果文件名重复，添加序号
                counter = 1
                while filename in seen_filenames:
                    filename = f"{base_filename}_{counter}.json"
                    counter += 1
                
                seen_filenames.add(filename)
                filepath = output_dir / filename
                
                try:
                    with open(filepath, 'w', encoding='utf-8') as f:
                        json.dump([sequence], f, indent=2, ensure_ascii=False)
                    saved_count += 1
                    
                    # 每100个文件报告进度
                    if (i + 1) % 100 == 0:
                        print(f"进度: {i + 1}/{len(sequences)} 已保存")
                        
                except Exception as e:
                    failed_count += 1
                    failed_files.append({
                        'unique_name': unique_name, 
                        'original_id': original_id,
                        'folder': folder,
                        'error': str(e)
                    })
                    print(f"保存文件出错 {filename}: {e}")
            
            # 最终报告
            print(f"\n保存完成统计:")
            print(f"- 成功保存: {saved_count} 个文件")
            print(f"- 保存失败: {failed_count} 个文件")
            print(f"- 保存位置: {output_dir}")
            print(f"- 每个序列都包含ATP配体信息")
            
            # 显示失败的文件
            if failed_files:
                print(f"\n失败的文件详情:")
                for fail_info in failed_files[:5]:  # 只显示前5个
                    print(f"  - {fail_info['folder']}/{fail_info['original_id']}: {fail_info['error']}")
                if len(failed_files) > 5:
                    print(f"  ... 还有 {len(failed_files) - 5} 个失败")
            
            # 验证实际文件数
            actual_files = list(output_dir.glob("*.json"))
            print(f"\n验证: 目录中实际有 {len(actual_files)} 个JSON文件")
            
            if len(actual_files) != saved_count:
                print(f"警告: 报告的成功数量({saved_count})与实际文件数({len(actual_files)})不符!")
            else:
                print("✅ 文件保存验证通过!")

def main():
    # 固定的文件路径
    excel_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/NTP_Analysis_Report.xlsx"
    fasta_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files"
    output_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/ntp_processing_sequences_with_ligand"
    
    # 检查文件是否存在
    if not os.path.exists(excel_path):
        print(f"错误: Excel文件不存在 - {excel_path}")
        sys.exit(1)
    
    if not os.path.exists(fasta_path):
        print(f"错误: FASTA基础路径不存在 - {fasta_path}")
        sys.exit(1)
    
    print("=" * 60)
    print("NTP Processing 序列提取器 - 添加配体版本")
    print("=" * 60)
    print(f"Excel文件: {excel_path}")
    print(f"FASTA路径: {fasta_path}")
    print(f"输出路径: {output_path}")
    print(f"配体信息: CCD_ATP (所有序列)")
    print("=" * 60)
    
    # 创建提取器并处理
    extractor = NTPProcessingSequenceExtractor()
    sequences = extractor.extract_ntp_processing_sequences()
    
    if sequences:
        extractor.save_sequences(sequences, output_path, single_file=False)
        print(f"\n✅ 处理完成! 成功提取 {len(sequences)} 个NTP processing序列")
        print(f"JSON文件已保存到: {output_path}")
        print(f"文件命名格式: sequenceID_folder.json")
        print(f"每个序列都包含ATP配体信息")
    else:
        print("\n❌ 没有成功提取任何序列")

if __name__ == "__main__":
    main()