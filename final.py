#!/usr/bin/env python3
"""
修改的脚本：根据Excel文件过滤CSV文件
1. 从CSV的filename列提取蛋白质ID
2. 只保留Excel中存在的蛋白质ID对应的行
3. 生成过滤后的新CSV文件
4. 报告空文件情况
"""

import pandas as pd
import os
import glob
import re
from pathlib import Path

def extract_protein_id_from_filename_entry(filename_entry):
    """
    从filename条目中提取蛋白质ID
    Example: 'a5ulA_vs_3amtA.txt' -> '3amtA'
    """
    filename = str(filename_entry)
    
    # 移除.txt扩展名
    if filename.endswith('.txt'):
        filename = filename[:-4]
    
    # 使用正则表达式提取蛋白质ID
    patterns = [
        r'_vs_([a-zA-Z0-9]{5})',  # 匹配_vs_后面的5个字符
        r'_vs_([a-zA-Z0-9]{4}[a-zA-Z])',  # 4个字符+1个字母
        r'_vs_([0-9][a-zA-Z0-9]{3}[a-zA-Z])',  # 数字开头的5字符模式
        r'_vs_([a-zA-Z0-9]{6})',  # 6个字符的情况
        r'_vs_([a-zA-Z0-9]+)',  # 任意长度的字符
    ]
    
    for pattern in patterns:
        match = re.search(pattern, filename)
        if match:
            return match.group(1)
    
    return None

def load_excel_protein_ids(excel_path):
    """
    从Excel文件的第一列加载蛋白质ID
    """
    try:
        print(f"正在读取Excel文件: {excel_path}")
        df = pd.read_excel(excel_path)
        print(f"Excel文件形状: {df.shape}")
        print(f"Excel列名: {df.columns.tolist()}")
        
        # 显示前几行数据
        print("Excel文件前5行:")
        print(df.head())
        
        # 获取第一列的值（Representative列），包含标题行后的所有数据，移除NaN值
        protein_ids = df.iloc[:, 0].dropna().tolist()
        
        # 清理数据：去除空格，转换为字符串，跳过列标题
        cleaned_ids = set()
        for i, pid in enumerate(protein_ids):
            # 跳过第一行的列标题 "Representative"
            if i == 0 and str(pid).strip().lower() in ['representative', 'protein_id', 'id']:
                continue
            
            cleaned_id = str(pid).strip()
            if cleaned_id:  # 只添加非空字符串
                cleaned_ids.add(cleaned_id)
        
        print(f"从Excel加载了 {len(cleaned_ids)} 个蛋白质ID")
        print(f"前10个蛋白质ID样本: {list(cleaned_ids)[:10]}")
        
        return cleaned_ids
    except Exception as e:
        print(f"读取Excel文件时出错: {e}")
        return set()

def process_single_csv(file_path, valid_protein_ids, output_directory):
    """
    处理单个CSV文件
    返回：(处理状态, 原始行数, 过滤后行数, 匹配的蛋白质ID列表)
    """
    filename = os.path.basename(file_path)
    
    try:
        # 读取CSV文件
        df = pd.read_csv(file_path)
        original_rows = len(df)
        
        # 检查是否有filename列
        if 'filename' not in df.columns:
            return 'no_filename_column', original_rows, 0, []
        
        # 如果原始文件就是空的
        if original_rows == 0:
            return 'originally_empty', 0, 0, []
        
        # 转换为小写便于比较
        valid_ids_lower = {vid.lower() for vid in valid_protein_ids}
        
        # 为每行提取蛋白质ID并检查是否匹配
        matched_indices = []
        matched_protein_ids = []
        
        for idx, row in df.iterrows():
            protein_id = extract_protein_id_from_filename_entry(row['filename'])
            if protein_id and protein_id.lower() in valid_ids_lower:
                matched_indices.append(idx)
                matched_protein_ids.append(protein_id)
        
        # 创建过滤后的DataFrame
        if matched_indices:
            filtered_df = df.iloc[matched_indices]
            filtered_rows = len(filtered_df)
            
            # 保存过滤后的文件
            output_path = os.path.join(output_directory, filename)
            filtered_df.to_csv(output_path, index=False)
            
            return 'filtered_success', original_rows, filtered_rows, matched_protein_ids
        else:
            # 没有匹配的行，创建空文件
            empty_df = df.iloc[0:0]  # 保持列结构但无数据
            output_path = os.path.join(output_directory, filename)
            empty_df.to_csv(output_path, index=False)
            
            return 'filtered_empty', original_rows, 0, []
    
    except Exception as e:
        print(f"处理文件 {filename} 时出错: {e}")
        return 'error', 0, 0, []

def process_csv_files(csv_directory, excel_path, output_directory=None):
    """
    处理目录中的所有CSV文件
    """
    # 从Excel加载有效的蛋白质ID
    valid_protein_ids = load_excel_protein_ids(excel_path)
    if not valid_protein_ids:
        print("错误：没有从Excel文件中加载到有效的蛋白质ID")
        return
    
    # 创建输出目录
    if output_directory:
        os.makedirs(output_directory, exist_ok=True)
    
    # 查找所有CSV和TXT文件
    file_patterns = ["*.csv", "*.txt"]
    all_files = []
    for pattern in file_patterns:
        files = glob.glob(os.path.join(csv_directory, pattern))
        all_files.extend(files)
    
    print(f"\n找到 {len(all_files)} 个文件需要处理")
    
    # 统计变量
    results = {
        'originally_empty': [],      # 原本就是空的
        'filtered_empty': [],        # 过滤后变空的
        'filtered_success': [],      # 成功过滤且有数据
        'no_filename_column': [],    # 没有filename列
        'error': []                  # 处理出错
    }
    
    # 处理每个文件
    for file_path in all_files:
        filename = os.path.basename(file_path)
        print(f"\n处理文件: {filename}")
        
        status, original_rows, filtered_rows, matched_ids = process_single_csv(
            file_path, valid_protein_ids, output_directory
        )
        
        file_info = {
            'filename': filename,
            'original_rows': original_rows,
            'filtered_rows': filtered_rows,
            'matched_ids': matched_ids
        }
        
        results[status].append(file_info)
        
        # 输出处理结果
        if status == 'originally_empty':
            print(f"  ⚠️  原本就是空文件")
        elif status == 'filtered_empty':
            print(f"  ❌ 过滤后变空 (原有 {original_rows} 行)")
        elif status == 'filtered_success':
            print(f"  ✅ 成功过滤: {original_rows} -> {filtered_rows} 行")
            print(f"     匹配的蛋白质ID: {matched_ids}")
        elif status == 'no_filename_column':
            print(f"  ❌ 没有filename列")
        elif status == 'error':
            print(f"  ❌ 处理出错")
    
    # 输出详细总结
    print(f"\n" + "="*80)
    print(f"处理完成总结")
    print(f"="*80)
    
    total_files = len(all_files)
    print(f"总文件数: {total_files}")
    print(f"成功过滤且有数据: {len(results['filtered_success'])}")
    print(f"过滤后变空: {len(results['filtered_empty'])}")
    print(f"原本就是空: {len(results['originally_empty'])}")
    print(f"没有filename列: {len(results['no_filename_column'])}")
    print(f"处理出错: {len(results['error'])}")
    
    # 保存详细报告
    if output_directory:
        report_file = os.path.join(output_directory, "processing_report.txt")
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("CSV文件处理报告\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"总文件数: {total_files}\n")
            f.write(f"成功过滤且有数据: {len(results['filtered_success'])}\n")
            f.write(f"过滤后变空: {len(results['filtered_empty'])}\n")
            f.write(f"原本就是空: {len(results['originally_empty'])}\n")
            f.write(f"没有filename列: {len(results['no_filename_column'])}\n")
            f.write(f"处理出错: {len(results['error'])}\n\n")
            
            # 详细列出需要关注的文件
            if results['filtered_empty']:
                f.write("⚠️ 过滤后变空的文件:\n")
                f.write("-" * 40 + "\n")
                for info in results['filtered_empty']:
                    f.write(f"  {info['filename']} (原有 {info['original_rows']} 行)\n")
                f.write("\n")
            
            if results['originally_empty']:
                f.write("⚠️ 原本就是空的文件:\n")
                f.write("-" * 40 + "\n")
                for info in results['originally_empty']:
                    f.write(f"  {info['filename']}\n")
                f.write("\n")
            
            if results['filtered_success']:
                f.write("✅ 成功过滤的文件:\n")
                f.write("-" * 40 + "\n")
                for info in results['filtered_success']:
                    f.write(f"  {info['filename']}: {info['original_rows']} -> {info['filtered_rows']} 行\n")
                    f.write(f"     匹配ID: {', '.join(info['matched_ids'])}\n")
                f.write("\n")
            
            if results['no_filename_column'] or results['error']:
                f.write("❌ 问题文件:\n")
                f.write("-" * 40 + "\n")
                for info in results['no_filename_column']:
                    f.write(f"  {info['filename']} (没有filename列)\n")
                for info in results['error']:
                    f.write(f"  {info['filename']} (处理出错)\n")
        
        print(f"\n详细报告已保存到: {report_file}")
    
    # 特别提醒空文件情况
    empty_files_count = len(results['filtered_empty']) + len(results['originally_empty'])
    if empty_files_count > 0:
        print(f"\n🚨 注意：有 {empty_files_count} 个文件是空的或过滤后变空！")
        print("详情请查看处理报告文件。")

def main():
    # 配置 - 更新这些路径
    csv_directory = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/batch_results"
    excel_path = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/SUPFAM based NTP processing.xlsx"
    
    # 输出目录以保持原始文件不变
    output_directory = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/filtered_results"
    
    # 验证路径是否存在
    if not os.path.exists(csv_directory):
        print(f"错误: CSV目录未找到: {csv_directory}")
        return
    
    if not os.path.exists(excel_path):
        print(f"错误: Excel文件未找到: {excel_path}")
        return
    
    print(f"CSV目录: {csv_directory}")
    print(f"Excel文件: {excel_path}")
    print(f"输出目录: {output_directory}")
    print("将根据Excel文件过滤CSV内容，并报告空文件情况")
    
    # 处理文件
    process_csv_files(csv_directory, excel_path, output_directory)

if __name__ == "__main__":
    main()