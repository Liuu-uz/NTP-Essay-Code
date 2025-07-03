#!/usr/bin/env python3
"""
序列转JSON生成脚本
从sequence_list.txt文件读取序列信息，生成JSON文件
"""

import json
import argparse
import sys
import os
from typing import Dict, List

class SequenceToJSONConverter:
    def __init__(self):
        pass
    
    def validate_sequence(self, sequence: str) -> bool:
        valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
        cleaned_sequence = sequence.replace(' ', '').replace('\n', '').upper()
        return all(aa in valid_amino_acids for aa in cleaned_sequence)
    
    def clean_sequence(self, sequence: str) -> str:
        return sequence.replace(' ', '').replace('\n', '').replace('\r', '').upper().strip()
    
    def process_sequence_list_file(self, txt_file: str) -> List[Dict]:
        try:
            with open(txt_file, 'r', encoding='utf-8') as f:
                content = f.read().strip()
            if not content:
                print("错误: 文件为空")
                return []
            results = []
            if content.startswith('>'):
                results = self._process_fasta_format(content)
            elif '|' in content or ':' in content:
                results = self._process_named_sequences(content)
            else:
                results = self._process_simple_sequences(content)
            return results
        except FileNotFoundError:
            print(f"错误: 文件 '{txt_file}' 不存在")
            return []
        except Exception as e:
            print(f"处理文件时出错: {e}")
            return []
    
    def _process_fasta_format(self, content: str) -> List[Dict]:
        results = []
        lines = content.split('\n')
        current_name = None
        current_sequence = ""
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if current_name and current_sequence:
                    results.append(self._create_sequence_entry(current_sequence, current_name))
                current_name = line[1:].strip()
                if not current_name:
                    current_name = f"sequence_{len(results) + 1}"
                current_sequence = ""
            elif line and not line.startswith('>'):
                current_sequence += line
        if current_name and current_sequence:
            results.append(self._create_sequence_entry(current_sequence, current_name))
        return results
    
    def _process_named_sequences(self, content: str) -> List[Dict]:
        results = []
        lines = content.split('\n')
        for i, line in enumerate(lines, 1):
            line = line.strip()
            if not line:
                continue
            name = None
            sequence = None
            if '|' in line:
                parts = line.split('|', 1)
                if len(parts) == 2:
                    name, sequence = parts
            elif ':' in line:
                parts = line.split(':', 1)
                if len(parts) == 2:
                    name, sequence = parts
            if name and sequence:
                name = name.strip()
                sequence = sequence.strip()
                if not name:
                    name = f"sequence_{i}"
                results.append(self._create_sequence_entry(sequence, name))
            else:
                print(f"警告: 第{i}行格式不正确，跳过: {line[:50]}...")
        return results
    
    def _process_simple_sequences(self, content: str) -> List[Dict]:
        results = []
        lines = content.split('\n')
        for i, line in enumerate(lines, 1):
            line = line.strip()
            if not line:
                continue
            name = f"sequence_{i}"
            results.append(self._create_sequence_entry(line, name))
        return results
    
    def _create_sequence_entry(self, sequence: str, name: str) -> Dict:
        cleaned_sequence = self.clean_sequence(sequence)
        if not cleaned_sequence:
            print(f"警告: 序列 '{name}' 为空，跳过")
            return None
        if not self.validate_sequence(cleaned_sequence):
            print(f"警告: 序列 '{name}' 包含无效字符，跳过")
            return None
        print(f"处理序列: {name} (长度: {len(cleaned_sequence)})")
        return {
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": cleaned_sequence,
                        "count": 1
                    }
                }
            ],
            "name": name.replace(' ', '_').replace('/', '_')
        }
    
    def process_single_sequence(self, sequence: str, name: str = "user_sequence") -> Dict:
        return self._create_sequence_entry(sequence, name)
    
    def save_to_json_folder(self, data: List[Dict], output_folder: str = "sequences_json"):
        try:
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
                print(f"创建文件夹: {output_folder}")
            saved_count = 0
            for entry in data:
                if entry is None:
                    continue
                filename = f"{entry['name']}.json"
                filepath = os.path.join(output_folder, filename)
                entry_data = [entry]
                with open(filepath, 'w', encoding='utf-8') as f:
                    json.dump(entry_data, f, indent=2, ensure_ascii=False)
                print(f"保存文件: {filepath}")
                saved_count += 1
            print(f"成功保存{saved_count}个JSON文件到文件夹: {output_folder}")
        except Exception as e:
            print(f"保存JSON文件时出错: {e}")
    
    def save_to_json(self, data: List[Dict], output_file: str):
        try:
            valid_data = [entry for entry in data if entry is not None]
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(valid_data, f, indent=2, ensure_ascii=False)
            print(f"成功保存到文件: {output_file}")
        except Exception as e:
            print(f"保存JSON文件时出错: {e}")

def main():
    parser = argparse.ArgumentParser(description='从序列文件生成JSON文件')
    parser.add_argument('input_file', nargs='?', default='sequence_list.txt',
                       help='输入文件路径 (默认: sequence_list.txt)')
    parser.add_argument('-f', '--folder', default='sequences_json', 
                       help='输出JSON文件夹路径 (默认: sequences_json)')
    parser.add_argument('-s', '--sequence', type=str,
                       help='直接提供一个序列字符串')
    parser.add_argument('-n', '--name', type=str, default='user_sequence',
                       help='序列名称 (与-s一起使用)')
    parser.add_argument('--single-file', action='store_true',
                       help='保存为单个JSON文件而不是分别保存到文件夹')
    
    args = parser.parse_args()
    
    converter = SequenceToJSONConverter()
    
    if args.sequence:
        print(f"处理单个序列: {args.name}")
        sequence_data = converter.process_single_sequence(args.sequence, args.name)
        if sequence_data:
            if args.single_file:
                output_file = f"{args.folder}.json" if not args.folder.endswith('.json') else args.folder
                converter.save_to_json([sequence_data], output_file)
            else:
                converter.save_to_json_folder([sequence_data], args.folder)
            print("处理完成！")
        else:
            print("序列处理失败")
    else:
        if not os.path.exists(args.input_file):
            print(f"错误: 文件 '{args.input_file}' 不存在")
            print("\n支持的文件格式:")
            print("1. 每行一个序列:")
            print("   MKLAVIFG...")
            print("   MKLVQTSR...")
            print("\n2. 名称|序列 或 名称:序列:")
            print("   protein1|MKLAVIFG...")
            print("   protein2:MKLVQTSR...")
            print("\n3. FASTA格式:")
            print("   >protein1")
            print("   MKLAVIFG...")
            print("   >protein2") 
            print("   MKLVQTSR...")
            sys.exit(1)
        
        print(f"开始处理文件: {args.input_file}")
        sequences_data = converter.process_sequence_list_file(args.input_file)
        
        if sequences_data:
            if args.single_file:
                output_file = f"{args.folder}.json" if not args.folder.endswith('.json') else args.folder
                converter.save_to_json(sequences_data, output_file)
                print(f"处理完成！成功处理了{len(sequences_data)}个序列，保存到单个文件")
            else:
                converter.save_to_json_folder(sequences_data, args.folder)
                print(f"处理完成！成功处理了{len(sequences_data)}个序列，分别保存到文件夹")
        else:
            print("未处理到任何有效序列")

if __name__ == "__main__":
    main()