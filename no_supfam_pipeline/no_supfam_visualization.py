import os
import csv
import pymol
from pymol import cmd
import glob

def parse_extracted_csv(csv_file_path):
    """解析extracted_first_rows.csv文件"""
    results = []
    try:
        with open(csv_file_path, 'r', encoding='utf-8') as file:
            reader = csv.DictReader(file)
            for row in reader:
                # 从Source_File提取查询ID (例如: Q67ZM7_zscores.csv -> Q67ZM7)
                query_id = row['Source_File'].replace('_zscores.csv', '')
                
                # 从filename的最后5个字母取前4个作为模板ID
                filename = row['filename']
                # 去掉.txt后缀
                filename_no_ext = filename.replace('.txt', '')
                
                if len(filename_no_ext) >= 5:
                    # 取最后5个字符，然后取前4个
                    last_5_chars = filename_no_ext[-5:]
                    template_id = last_5_chars[:4]
                    
                    print(f"  文件名: {filename} -> 最后5个字符: {last_5_chars} -> 模板ID: {template_id}")
                    
                    results.append({
                        'query_id': query_id,
                        'template_id': template_id,
                        'z_score': float(row['Z-score']),
                        'filename': filename
                    })
                else:
                    print(f"  警告: 文件名 {filename} 长度不足5个字符")
                    
        return results
    except Exception as e:
        print(f"解析CSV文件时出错: {e}")
        return []

def find_query_pdb(predicted_dir, query_id):
    """在predicted_structures目录中查找查询结构PDB文件"""
    # 尝试多种可能的文件名格式
    possible_patterns = [
        f"{query_id.lower()}*.pdb",
        f"{query_id[:4].lower()}*.pdb",
        f"{query_id}*.pdb"
    ]
    
    for pattern in possible_patterns:
        query_pattern = os.path.join(predicted_dir, pattern)
        query_files = glob.glob(query_pattern)
        if query_files:
            return query_files[0]
    
    return None

def find_template_pdb(template_dir, template_id):
    """在template目录中查找模板PDB文件"""
    if not os.path.exists(template_dir):
        print(f"    模板目录不存在: {template_dir}")
        return None
    
    # 尝试多种文件名格式
    possible_patterns = [
        f"{template_id.lower()}*.pdb",
        f"{template_id.upper()}*.pdb", 
        f"{template_id}*.pdb"
    ]
    
    for pattern in possible_patterns:
        template_pattern = os.path.join(template_dir, pattern)
        template_files = glob.glob(template_pattern)
        if template_files:
            print(f"    找到模板文件: {template_files[0]} (使用模式: {pattern})")
            return template_files[0]
    
    print(f"    未找到模板文件，尝试的模式: {possible_patterns}")
    return None

def setup_pymol_session():
    """初始化PyMOL会话"""
    pymol.finish_launching(['pymol', '-qc'])
    cmd.reinitialize()

def visualize_structures_with_atp(query_pdb, template_pdb, output_prefix, query_id, template_id, z_score):
    """使用PyMOL可视化结构，包括ATP和活性位点，无文字标签"""
    
    # 清除当前会话
    cmd.delete('all')
    
    # 加载结构文件
    query_name = 'query_structure'
    template_name = 'template_structure'
    
    try:
        cmd.load(query_pdb, query_name)
        cmd.load(template_pdb, template_name)
    except Exception as e:
        print(f"    加载PDB文件失败: {e}")
        return None
    
    # 结构叠合 - 使用CA原子进行对齐
    try:
        alignment_result = cmd.align(query_name, template_name)
        print(f"    结构对齐结果: RMSD = {alignment_result[0]:.3f} Å, 对齐原子数 = {alignment_result[1]}")
    except Exception as e:
        print(f"    结构对齐失败: {e}")
        return None
    
    # 设置结构显示样式
    cmd.hide('everything')
    cmd.show('cartoon', 'all')
    cmd.color('cyan', query_name)
    cmd.color('orange', template_name)
    
    # 查找并显示ATP分子
    atp_selection = 'resn ATP or resn ADP or resn AMP'
    nucleotide_found = False
    
    if cmd.count_atoms(atp_selection) > 0:
        cmd.show('sticks', atp_selection)
        cmd.color('red', atp_selection)
        print("    发现ATP/ADP/AMP分子")
        nucleotide_found = True
        
        # 显示ATP周围的活性位点残基（5Å范围内）
        active_site = f'byres ({atp_selection} around 5)'
        cmd.show('sticks', active_site)
        cmd.color('yellow', f'{active_site} and not {atp_selection}')
        print("    显示活性位点残基")
    else:
        # 如果没有ATP，尝试显示其他核苷酸
        other_nucleotides = 'resn GTP or resn GDP or resn GMP or resn CTP or resn CDP or resn CMP or resn UTP or resn UDP or resn UMP'
        if cmd.count_atoms(other_nucleotides) > 0:
            cmd.show('sticks', other_nucleotides)
            cmd.color('red', other_nucleotides)
            active_site = f'byres ({other_nucleotides} around 5)'
            cmd.show('sticks', active_site)
            cmd.color('yellow', f'{active_site} and not {other_nucleotides}')
            print("    发现其他核苷酸分子")
            nucleotide_found = True
        else:
            print("    未发现核苷酸分子")
    
    # 设置视图
    cmd.zoom('all')
    cmd.orient()
    
    # 设置高质量渲染参数
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_shadows', 1)  # 启用阴影增加立体感
    cmd.set('ambient', 0.2)
    cmd.set('direct', 0.8)
    cmd.set('reflect', 0.5)
    cmd.set('shininess', 50)
    cmd.set('spec_reflect', 0.8)
    
    # 保存白色背景版本
    cmd.bg_color('white')
    png_file_white = f'{output_prefix}_clean_white.png'
    cmd.png(png_file_white, width=1800, height=1400, dpi=300, ray=1)
    
    
    # 保存叠合后的结构
    pdb_file = f'{output_prefix}_aligned.pdb'
    cmd.save(pdb_file, 'all')
    
    # 添加标题信息到终端输出（不显示在图像中）
    title_text = f"Query: {query_id} vs Template: {template_id} (Z-score: {z_score:.2f})"
    print(f"    {title_text}")
    
    return {
        'rmsd': alignment_result[0],
        'aligned_atoms': alignment_result[1],
        'nucleotide_found': nucleotide_found,
        'png_file_white': png_file_white,
        'pdb_file': pdb_file
    }

def process_csv_data(csv_file_path, predicted_dir, template_dir, output_dir, max_structures=None):
    """基于CSV数据处理结构可视化"""
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 初始化PyMOL
    setup_pymol_session()
    
    # 解析CSV文件
    csv_data = parse_extracted_csv(csv_file_path)
    if not csv_data:
        print("无法解析CSV文件或文件为空")
        return []
    
    print(f"从CSV文件中解析到 {len(csv_data)} 条记录")
    
    # 如果设置了最大处理数量，则只处理前N个
    if max_structures:
        csv_data = csv_data[:max_structures]
        print(f"限制处理前 {max_structures} 个结构")
    
    results = []
    successful_count = 0
    
    for i, data in enumerate(csv_data, 1):
        query_id = data['query_id']
        template_id = data['template_id']
        z_score = data['z_score']
        filename = data['filename']
        
        print(f"\n[{i}/{len(csv_data)}] 处理: {query_id} vs {template_id} (Z-score: {z_score:.2f})")
        print(f"  原始文件名: {filename}")
        
        # 查找查询结构PDB文件
        query_pdb = find_query_pdb(predicted_dir, query_id)
        if not query_pdb:
            print(f"    未找到查询结构: {query_id}")
            continue
        
        # 查找模板结构PDB文件
        template_pdb = find_template_pdb(template_dir, template_id)
        if not template_pdb:
            print(f"    未找到模板结构: {template_id}")
            continue
        
        print(f"    查询结构: {query_pdb}")
        print(f"    模板结构: {template_pdb}")
        
        # 创建输出文件前缀
        output_prefix = os.path.join(output_dir, f"{query_id}_vs_{template_id}_zscore{z_score:.1f}")
        
        try:
            # 进行结构可视化
            vis_result = visualize_structures_with_atp(
                query_pdb, template_pdb, output_prefix, 
                query_id, template_id, z_score
            )
            
            if vis_result:
                # 记录结果
                result = {
                    'query_id': query_id,
                    'template_id': template_id,
                    'z_score': z_score,
                    'filename': filename,
                    'query_pdb': query_pdb,
                    'template_pdb': template_pdb,
                    'rmsd': vis_result['rmsd'],
                    'aligned_atoms': vis_result['aligned_atoms'],
                    'nucleotide_found': vis_result['nucleotide_found'],
                    'png_file_white': vis_result['png_file_white'],
                    'pdb_file': vis_result['pdb_file']
                }
                results.append(result)
                successful_count += 1
                print(f"    ✓ 成功处理 ({successful_count}/{i})")
            
        except Exception as e:
            print(f"    ✗ 处理失败: {e}")
            continue
    
    # 保存处理结果摘要
    summary_file = os.path.join(output_dir, 'visualization_summary.csv')
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Query_ID', 'Template_ID', 'Z_Score', 'Original_Filename', 'Query_PDB', 'Template_PDB', 
            'RMSD', 'Aligned_Atoms', 'Nucleotide_Found', 'PNG_Clean_White', 'PDB_File'
        ])
        
        for result in results:
            writer.writerow([
                result['query_id'],
                result['template_id'],
                f"{result['z_score']:.2f}",
                result['filename'],
                result['query_pdb'],
                result['template_pdb'],
                f"{result['rmsd']:.3f}",
                result['aligned_atoms'],
                result['nucleotide_found'],
                os.path.basename(result['png_file_white']),
                os.path.basename(result['pdb_file'])
            ])
    
    print(f"\n🎉 处理完成！")
    print(f"总记录数: {len(csv_data)}")
    print(f"成功处理: {successful_count}")
    print(f"失败数量: {len(csv_data) - successful_count}")
    print(f"结果摘要保存在: {summary_file}")
    
    return results

def main():
    """主函数"""
    # 设置路径
    base_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code"
    pipeline_dir = os.path.join(base_dir, "no_supfam_pipeline")
    
    csv_file_path = os.path.join(pipeline_dir, "extracted_first_rows.csv")
    predicted_dir = os.path.join(pipeline_dir, "predicted_structures")
    template_dir = os.path.join(pipeline_dir, "template")  # 明确指定template目录
    output_dir = os.path.join(pipeline_dir, "pymol_visualization_output")
    
    # 检查文件和目录是否存在
    if not os.path.exists(csv_file_path):
        print(f"错误: CSV文件不存在: {csv_file_path}")
        return
    
    if not os.path.exists(predicted_dir):
        print(f"错误: 预测结构目录不存在: {predicted_dir}")
        return
    
    if not os.path.exists(template_dir):
        print(f"错误: 模板目录不存在: {template_dir}")
        return
    
    print("基于CSV文件进行PyMOL结构可视化 (使用filename最后5个字符的前4个作为模板ID)")
    print(f"CSV文件: {csv_file_path}")
    print(f"预测结构目录: {predicted_dir}")
    print(f"模板目录: {template_dir}")
    print(f"输出目录: {output_dir}")
    
    # 询问是否限制处理数量（用于测试）
    response = input("\n是否限制处理数量？(输入数字限制，回车处理所有): ").strip()
    max_structures = None
    if response.isdigit():
        max_structures = int(response)
        print(f"将只处理前 {max_structures} 个结构")
    
    # 处理文件
    results = process_csv_data(csv_file_path, predicted_dir, template_dir, output_dir, max_structures)
    
    # 显示统计信息
    if results:
        z_scores = [r['z_score'] for r in results]
        rmsds = [r['rmsd'] for r in results]
        nucleotide_count = sum(1 for r in results if r['nucleotide_found'])
        
        print(f"\n📊 统计信息:")
        print(f"Z-score范围: {min(z_scores):.2f} - {max(z_scores):.2f}")
        print(f"平均Z-score: {sum(z_scores)/len(z_scores):.2f}")
        print(f"RMSD范围: {min(rmsds):.3f} - {max(rmsds):.3f} Å")
        print(f"平均RMSD: {sum(rmsds)/len(rmsds):.3f} Å")
        print(f"包含核苷酸的结构: {nucleotide_count}/{len(results)} ({nucleotide_count/len(results)*100:.1f}%)")
    
    # 退出PyMOL
    cmd.quit()

if __name__ == "__main__":
    main()
