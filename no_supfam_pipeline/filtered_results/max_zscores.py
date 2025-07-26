import os
import csv
import glob

def extract_first_data_row_from_csv(csv_file_path):
    """从CSV文件中提取第一行数据（跳过标题行）"""
    try:
        with open(csv_file_path, 'r', encoding='utf-8') as file:
            reader = csv.reader(file)
            
            # 跳过标题行
            try:
                header = next(reader)
                print(f"  标题行: {header}")
            except StopIteration:
                print(f"  文件为空: {csv_file_path}")
                return None, None
            
            # 读取第一行数据
            try:
                first_data_row = next(reader)
                return header, first_data_row
            except StopIteration:
                print(f"  只有标题行，没有数据: {csv_file_path}")
                return header, None
                
    except Exception as e:
        print(f"读取文件 {csv_file_path} 时出错: {e}")
        return None, None

def extract_first_rows_batch(input_directory, output_file):
    """批量提取所有CSV文件的第一行数据并合并"""
    
    # 查找所有CSV文件
    csv_pattern = os.path.join(input_directory, "*.csv")
    csv_files = glob.glob(csv_pattern)
    
    if not csv_files:
        print(f"在目录 {input_directory} 中未找到CSV文件")
        return
    
    print(f"找到 {len(csv_files)} 个CSV文件")
    
    all_headers = set()
    successful_extractions = []
    failed_files = []
    
    # 首先收集所有可能的列标题
    for csv_file in csv_files:
        header, _ = extract_first_data_row_from_csv(csv_file)
        if header:
            all_headers.update(header)
    
    # 创建统一的列标题列表
    unified_headers = ['Source_File'] + sorted(list(all_headers))
    
    # 处理每个文件
    for csv_file in csv_files:
        print(f"\n处理文件: {os.path.basename(csv_file)}")
        
        header, first_data_row = extract_first_data_row_from_csv(csv_file)
        
        if first_data_row:
            # 创建一个字典来存储这个文件的数据
            row_data = {'Source_File': os.path.basename(csv_file)}
            
            # 将数据按列标题对应
            for i, value in enumerate(first_data_row):
                if i < len(header):
                    row_data[header[i]] = value.strip() if value else ''
            
            successful_extractions.append(row_data)
            print(f"  成功提取第一行数据: {first_data_row}")
        else:
            failed_files.append(os.path.basename(csv_file))
    
    # 保存结果
    if successful_extractions:
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=unified_headers)
            
            # 写入标题行
            writer.writeheader()
            
            # 写入所有数据
            for row_data in successful_extractions:
                writer.writerow(row_data)
        
        print(f"\n提取完成！")
        print(f"成功处理了 {len(successful_extractions)} 个文件")
        print(f"失败的文件: {len(failed_files)}")
        if failed_files:
            print(f"失败文件列表: {', '.join(failed_files)}")
        print(f"结果保存在: {output_file}")
        
    else:
        print("没有成功提取到任何数据")

def simple_extract_first_rows(input_directory, output_file):
    """简化版本：直接提取每个文件的第一行数据"""
    
    csv_files = glob.glob(os.path.join(input_directory, "*.csv"))
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        
        # 写入表头
        first_file = True
        
        for csv_file in csv_files:
            file_name = os.path.basename(csv_file)
            print(f"处理文件: {file_name}")
            
            try:
                with open(csv_file, 'r', encoding='utf-8') as infile:
                    reader = csv.reader(infile)
                    
                    # 跳过标题行
                    try:
                        header = next(reader)
                        if first_file:
                            # 第一个文件时写入列标题
                            writer.writerow(['Source_File'] + header)
                            first_file = False
                    except StopIteration:
                        print(f"  文件为空: {file_name}")
                        continue
                    
                    # 读取第一行数据
                    try:
                        first_data_row = next(reader)
                        writer.writerow([file_name] + first_data_row)
                        print(f"  提取数据: {first_data_row}")
                    except StopIteration:
                        print(f"  只有标题，无数据: {file_name}")
                        
            except Exception as e:
                print(f"  处理 {file_name} 时出错: {e}")
    
    print(f"\n提取完成，结果保存在: {output_file}")

def extract_max_score_rows(input_directory, output_file, score_column_index=1):
    """提取每个文件中分数最高的行（假设第二列是分数）"""
    
    csv_files = glob.glob(os.path.join(input_directory, "*.csv"))
    
    results = []
    
    for csv_file in csv_files:
        file_name = os.path.basename(csv_file)
        print(f"处理文件: {file_name}")
        
        try:
            with open(csv_file, 'r', encoding='utf-8') as infile:
                reader = csv.reader(infile)
                
                # 读取标题行
                try:
                    header = next(reader)
                except StopIteration:
                    print(f"  文件为空: {file_name}")
                    continue
                
                # 找到分数最高的行
                max_score = float('-inf')
                max_row = None
                
                for row in reader:
                    if len(row) > score_column_index:
                        try:
                            score = float(row[score_column_index])
                            if score > max_score:
                                max_score = score
                                max_row = row
                        except (ValueError, IndexError):
                            continue
                
                if max_row:
                    result = {
                        'Source_File': file_name,
                        'Max_Score': max_score
                    }
                    
                    # 添加所有列的数据
                    for i, value in enumerate(max_row):
                        if i < len(header):
                            result[header[i]] = value.strip() if value else ''
                    
                    results.append(result)
                    print(f"  最高分数: {max_score}, 数据: {max_row}")
                else:
                    print(f"  未找到有效的分数数据")
                    
        except Exception as e:
            print(f"  处理 {file_name} 时出错: {e}")
    
    # 保存结果
    if results:
        # 获取所有可能的列名
        all_columns = set()
        for result in results:
            all_columns.update(result.keys())
        
        fieldnames = ['Source_File', 'Max_Score'] + sorted([col for col in all_columns if col not in ['Source_File', 'Max_Score']])
        
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        
        print(f"\n提取完成！找到 {len(results)} 个最高分数行")
        print(f"结果保存在: {output_file}")
    else:
        print("没有找到任何有效数据")

def main():
    """主函数"""
    # 设置输入和输出路径
    input_directory = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/no_supfam_pipeline/filtered_results"
    output_file = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/no_supfam_pipeline/extracted_first_rows.csv"
    
    # 检查输入目录是否存在
    if not os.path.exists(input_directory):
        print(f"错误: 输入目录不存在: {input_directory}")
        return
    
    print("提取CSV文件的第一行数据（跳过标题行）...")
    print(f"输入目录: {input_directory}")
    print(f"输出文件: {output_file}")
    
    # 选择提取方式
    print("\n选择提取方式:")
    print("1. 提取每个文件的第一行数据（跳过标题）")
    print("2. 提取每个文件分数最高的行（假设第二列是分数）")
    print("3. 简化版提取（快速）")
    
    choice = input("请输入选择 (1/2/3): ").strip()
    
    if choice == '1':
        extract_first_rows_batch(input_directory, output_file)
    elif choice == '2':
        extract_max_score_rows(input_directory, output_file)
    elif choice == '3':
        simple_extract_first_rows(input_directory, output_file)
    else:
        print("无效选择，默认执行方式1")
        extract_first_rows_batch(input_directory, output_file)

if __name__ == "__main__":
    main()