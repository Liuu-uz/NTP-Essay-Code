import os
import csv
import pymol
from pymol import cmd
import glob
import pandas as pd

def parse_xlsx_file(xlsx_file_path):
    """解析xlsx文件，提取EC number、sequenceID和PDBID+链信息"""
    results = []
    try:
        # 读取xlsx文件
        df = pd.read_excel(xlsx_file_path)
        
        print(f"读取到 {len(df)} 行数据")
        print(f"列名: {list(df.columns)}")
        
        for index, row in df.iterrows():
            try:
                # 第一列：EC number
                ec_number = str(row.iloc[0]).strip()
                
                # 第三列：sequenceID  
                sequence_id = str(row.iloc[2]).strip()
                
                # 第七列：匹配的PDBID+链
                pdb_chain = str(row.iloc[6]).strip()
                
                # 验证数据有效性
                if ec_number and sequence_id and pdb_chain and pdb_chain != 'nan':
                    # 提取PDB ID（前4个字符）和链ID
                    if len(pdb_chain) >= 4:
                        pdb_id = pdb_chain[:4].lower()
                        chain_id = pdb_chain[4:] if len(pdb_chain) > 4 else 'A'
                        
                        results.append({
                            'ec_number': ec_number,
                            'sequence_id': sequence_id,
                            'pdb_chain': pdb_chain,
                            'pdb_id': pdb_id,
                            'chain_id': chain_id
                        })
                        
                        print(f"  解析: EC={ec_number}, SeqID={sequence_id}, PDB={pdb_id}, Chain={chain_id}")
                
            except Exception as e:
                print(f"  跳过第 {index+1} 行，解析错误: {e}")
                continue
                
        return results
        
    except Exception as e:
        print(f"解析xlsx文件时出错: {e}")
        return []

def find_prediction_pdb(predicted_dir, sequence_id, ec_number):
    """在supfam_pdb_files目录中查找预测结构PDB文件"""
    # 文件命名规则：sequenceID_ECnumber.pdb
    expected_filename = f"{sequence_id}_{ec_number}.pdb"
    expected_path = os.path.join(predicted_dir, expected_filename)
    
    if os.path.exists(expected_path):
        print(f"    找到预测文件: {expected_filename}")
        return expected_path
    
    # 尝试其他可能的格式
    possible_patterns = [
        f"{sequence_id}_{ec_number}*.pdb",
        f"{sequence_id}_*.pdb",
        f"*{sequence_id}*.pdb"
    ]
    
    for pattern in possible_patterns:
        search_pattern = os.path.join(predicted_dir, pattern)
        files = glob.glob(search_pattern)
        if files:
            print(f"    找到预测文件: {os.path.basename(files[0])} (使用模式: {pattern})")
            return files[0]
    
    print(f"    未找到预测文件，期望: {expected_filename}")
    return None

def find_template_pdb(template_dir, pdb_id, chain_id):
    """在pdb_files目录中查找模板PDB文件"""
    if not os.path.exists(template_dir):
        print(f"    模板目录不存在: {template_dir}")
        return None
    
    # 尝试多种文件名格式
    possible_patterns = [
        f"{pdb_id}.pdb",
        f"{pdb_id.upper()}.pdb",
        f"pdb{pdb_id}.ent",
        f"{pdb_id}*.pdb"
    ]
    
    for pattern in possible_patterns:
        template_pattern = os.path.join(template_dir, pattern)
        template_files = glob.glob(template_pattern)
        if template_files:
            print(f"    找到模板文件: {template_files[0]} (使用模式: {pattern})")
            return template_files[0]
    
    print(f"    未找到模板文件，尝试的PDB ID: {pdb_id}")
    return None

def setup_pymol_session():
    """初始化PyMOL会话"""
    pymol.finish_launching(['pymol', '-qc'])
    cmd.reinitialize()

def extract_chain_from_pdb(pdb_file, chain_id, output_path):
    """从PDB文件中提取指定链并保存"""
    try:
        # 加载PDB文件
        temp_name = 'temp_structure'
        cmd.load(pdb_file, temp_name)
        
        # 检查指定链是否存在
        chains = cmd.get_chains(temp_name)
        print(f"    模板文件中的链: {chains}")
        
        if chain_id.upper() in [c.upper() for c in chains]:
            # 选择指定链
            chain_selection = f'{temp_name} and chain {chain_id.upper()}'
            cmd.create('extracted_chain', chain_selection)
            
            # 保存提取的链
            cmd.save(output_path, 'extracted_chain')
            
            # 清理
            cmd.delete(temp_name)
            cmd.delete('extracted_chain')
            
            print(f"    成功提取链 {chain_id.upper()} 到: {output_path}")
            return True
        else:
            print(f"    警告: 链 {chain_id.upper()} 不存在，使用完整结构")
            cmd.save(output_path, temp_name)
            cmd.delete(temp_name)
            return True
            
    except Exception as e:
        print(f"    提取链失败: {e}")
        return False

def visualize_structures_with_atp(query_pdb, template_pdb, output_prefix, sequence_id, pdb_id, chain_id, ec_number):
    """使用PyMOL可视化结构，包括ATP和活性位点"""
    
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
    
    # 关闭cartoon的箭头/方向性显示
    cmd.set('cartoon_fancy_helices', 0)  # 关闭螺旋的花式显示
    cmd.set('cartoon_fancy_sheets', 0)   # 关闭片层的花式显示
    cmd.set('cartoon_flat_sheets', 1)    # 使用平坦的片层
    cmd.set('cartoon_smooth_loops', 1)   # 平滑的loop区域
    cmd.set('cartoon_highlight_color', 'grey50')  # 设置高亮颜色
    cmd.set('cartoon_putty_transform', 0)  # 关闭putty变换
    cmd.set('cartoon_ring_mode', 0)      # 关闭环形模式
    cmd.set('cartoon_tube_radius', 0.4)  # 设置管状半径
    cmd.set('cartoon_oval_length', 1.2)  # 设置椭圆长度
    cmd.set('cartoon_oval_width', 0.3)   # 设置椭圆宽度
    
    # 查找并显示ATP分子
    atp_selection = 'resn ATP or resn ADP or resn AMP'
    nucleotide_found = False
    
    if cmd.count_atoms(atp_selection) > 0:
        cmd.show('sticks', atp_selection)
        cmd.color('red', atp_selection)
        cmd.show('spheres', f'{atp_selection} and name P*')
        cmd.set('sphere_scale', 0.3, f'{atp_selection} and name P*')
        print("    发现ATP/ADP/AMP分子")
        nucleotide_found = True
        
        # 显示ATP周围的活性位点残基（4Å范围内）
        active_site = f'byres ({atp_selection} around 4)'
        cmd.show('sticks', f'{active_site} and sidechain')
        cmd.color('yellow', f'{active_site} and sidechain and not {atp_selection}')
        print("    显示活性位点残基")
    else:
        # 如果没有ATP，尝试显示其他核苷酸
        other_nucleotides = 'resn GTP or resn GDP or resn GMP or resn CTP or resn CDP or resn CMP or resn UTP or resn UDP or resn UMP'
        if cmd.count_atoms(other_nucleotides) > 0:
            cmd.show('sticks', other_nucleotides)
            cmd.color('red', other_nucleotides)
            cmd.show('spheres', f'{other_nucleotides} and name P*')
            cmd.set('sphere_scale', 0.3, f'{other_nucleotides} and name P*')
            active_site = f'byres ({other_nucleotides} around 4)'
            cmd.show('sticks', f'{active_site} and sidechain')
            cmd.color('yellow', f'{active_site} and sidechain and not {other_nucleotides}')
            print("    发现其他核苷酸分子")
            nucleotide_found = True
        else:
            print("    未发现核苷酸分子")
    
    # 设置视图
    cmd.zoom('all')
    cmd.orient()
    
    # 关闭所有可能产生箭头的设置
    cmd.set('cgo_line_width', 1.0)
    cmd.set('dash_gap', 0.0)
    cmd.set('dash_length', 0.25)
    cmd.set('dash_round_ends', 1)
    cmd.hide('labels')  # 隐藏所有标签
    cmd.hide('nonbonded')  # 隐藏非键合原子
    
    # 设置高质量渲染参数，避免箭头显示
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_shadows', 1)
    cmd.set('ambient', 0.2)
    cmd.set('direct', 0.8)
    cmd.set('reflect', 0.5)
    cmd.set('shininess', 50)
    cmd.set('spec_reflect', 0.8)
    cmd.set('ray_opaque_background', 1)  # 确保背景不透明
    
    # 保存图像
    cmd.bg_color('white')
    png_file = f'{output_prefix}_visualization.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    
    # 保存叠合后的结构
    pdb_file = f'{output_prefix}_aligned.pdb'
    cmd.save(pdb_file, 'all')
    
    # 输出信息
    title_text = f"Query: {sequence_id} (EC: {ec_number}) vs Template: {pdb_id}_{chain_id}"
    print(f"    {title_text}")
    
    return {
        'rmsd': alignment_result[0],
        'aligned_atoms': alignment_result[1],
        'nucleotide_found': nucleotide_found,
        'png_file': png_file,
        'pdb_file': pdb_file
    }

def process_xlsx_data(xlsx_file_path, predicted_dir, template_dir, output_dir, max_structures=None):
    """基于xlsx数据处理结构可视化"""
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 初始化PyMOL
    setup_pymol_session()
    
    # 解析xlsx文件
    xlsx_data = parse_xlsx_file(xlsx_file_path)
    if not xlsx_data:
        print("无法解析xlsx文件或文件为空")
        return []
    
    print(f"从xlsx文件中解析到 {len(xlsx_data)} 条记录")
    
    # 如果设置了最大处理数量，则只处理前N个
    if max_structures:
        xlsx_data = xlsx_data[:max_structures]
        print(f"限制处理前 {max_structures} 个结构")
    
    results = []
    successful_count = 0
    
    for i, data in enumerate(xlsx_data, 1):
        sequence_id = data['sequence_id']
        ec_number = data['ec_number']
        pdb_id = data['pdb_id']
        chain_id = data['chain_id']
        pdb_chain = data['pdb_chain']
        
        print(f"\n[{i}/{len(xlsx_data)}] 处理: {sequence_id} (EC: {ec_number}) vs {pdb_id}_{chain_id}")
        
        # 查找预测结构PDB文件
        query_pdb = find_prediction_pdb(predicted_dir, sequence_id, ec_number)
        if not query_pdb:
            print(f"    未找到预测结构: {sequence_id}_{ec_number}")
            continue
        
        # 查找模板结构PDB文件
        template_pdb_raw = find_template_pdb(template_dir, pdb_id, chain_id)
        if not template_pdb_raw:
            print(f"    未找到模板结构: {pdb_id}")
            continue
        
        # 从模板PDB中提取指定链
        template_pdb = os.path.join(output_dir, f"temp_{pdb_id}_{chain_id}.pdb")
        if not extract_chain_from_pdb(template_pdb_raw, chain_id, template_pdb):
            print(f"    提取模板链失败: {pdb_id}_{chain_id}")
            continue
        
        print(f"    预测结构: {query_pdb}")
        print(f"    模板结构: {template_pdb}")
        
        # 创建输出文件前缀
        output_prefix = os.path.join(output_dir, f"{sequence_id}_{ec_number}_vs_{pdb_id}_{chain_id}")
        
        try:
            # 进行结构可视化
            vis_result = visualize_structures_with_atp(
                query_pdb, template_pdb, output_prefix, 
                sequence_id, pdb_id, chain_id, ec_number
            )
            
            if vis_result:
                # 记录结果
                result = {
                    'sequence_id': sequence_id,
                    'ec_number': ec_number,
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'pdb_chain': pdb_chain,
                    'query_pdb': query_pdb,
                    'template_pdb': template_pdb,
                    'rmsd': vis_result['rmsd'],
                    'aligned_atoms': vis_result['aligned_atoms'],
                    'nucleotide_found': vis_result['nucleotide_found'],
                    'png_file': vis_result['png_file'],
                    'pdb_file': vis_result['pdb_file']
                }
                results.append(result)
                successful_count += 1
                print(f"    ✓ 成功处理 ({successful_count}/{i})")
            
        except Exception as e:
            print(f"    ✗ 处理失败: {e}")
            continue
        
        finally:
            # 清理临时文件
            if os.path.exists(template_pdb):
                try:
                    os.remove(template_pdb)
                except:
                    pass
    
    # 保存处理结果摘要
    summary_file = os.path.join(output_dir, 'visualization_summary.csv')
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Sequence_ID', 'EC_Number', 'PDB_ID', 'Chain_ID', 'PDB_Chain', 'Query_PDB', 'Template_PDB',
            'RMSD', 'Aligned_Atoms', 'Nucleotide_Found', 'PNG_File', 'PDB_File'
        ])
        
        for result in results:
            writer.writerow([
                result['sequence_id'],
                result['ec_number'],
                result['pdb_id'],
                result['chain_id'],
                result['pdb_chain'],
                result['query_pdb'],
                result['template_pdb'],
                f"{result['rmsd']:.3f}",
                result['aligned_atoms'],
                result['nucleotide_found'],
                os.path.basename(result['png_file']),
                os.path.basename(result['pdb_file'])
            ])
    
    print(f"\n🎉 处理完成！")
    print(f"总记录数: {len(xlsx_data)}")
    print(f"成功处理: {successful_count}")
    print(f"失败数量: {len(xlsx_data) - successful_count}")
    print(f"结果摘要保存在: {summary_file}")
    
    return results

def main():
    """主函数"""
    # 设置路径
    base_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code"
    
    xlsx_file_path = os.path.join(base_dir, "seq-supfam", "NTP_Analysis_Report.xlsx")
    predicted_dir = os.path.join(base_dir, "supfam_pdb_files")
    template_dir = os.path.join(base_dir, "pdb_files")
    output_dir = os.path.join(base_dir, "structure_visualization_output")
    
    # 检查文件和目录是否存在
    if not os.path.exists(xlsx_file_path):
        print(f"错误: xlsx文件不存在: {xlsx_file_path}")
        return
    
    if not os.path.exists(predicted_dir):
        print(f"错误: 预测结构目录不存在: {predicted_dir}")
        return
    
    if not os.path.exists(template_dir):
        print(f"错误: 模板目录不存在: {template_dir}")
        return
    
    print("基于xlsx文件进行PyMOL结构可视化")
    print(f"xlsx文件: {xlsx_file_path}")
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
    results = process_xlsx_data(xlsx_file_path, predicted_dir, template_dir, output_dir, max_structures)
    
    # 显示统计信息
    if results:
        rmsds = [r['rmsd'] for r in results]
        nucleotide_count = sum(1 for r in results if r['nucleotide_found'])
        
        print(f"\n📊 统计信息:")
        print(f"RMSD范围: {min(rmsds):.3f} - {max(rmsds):.3f} Å")
        print(f"平均RMSD: {sum(rmsds)/len(rmsds):.3f} Å")
        print(f"包含核苷酸的结构: {nucleotide_count}/{len(results)} ({nucleotide_count/len(results)*100:.1f}%)")
    
    # 退出PyMOL
    cmd.quit()

if __name__ == "__main__":
    main()