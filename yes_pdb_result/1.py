import os
import csv
import pymol
from pymol import cmd
import glob
import pandas as pd
import numpy as np
import math

def load_xlsx_data(xlsx_file_path):
    """加载xlsx文件数据并建立索引"""
    try:
        df = pd.read_excel(xlsx_file_path)
        print(f"读取到 {len(df)} 行数据")
        print(f"列名: {list(df.columns)}")
        
        # 建立查找索引：{(sequence_id, ec_number): template_info}
        template_index = {}
        
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
                        
                        # 使用(sequence_id, ec_number)作为键
                        key = (sequence_id, ec_number)
                        template_index[key] = {
                            'ec_number': ec_number,
                            'sequence_id': sequence_id,
                            'pdb_chain': pdb_chain,
                            'pdb_id': pdb_id,
                            'chain_id': chain_id
                        }
                        
                        print(f"  索引: {key} -> PDB={pdb_id}, Chain={chain_id}")
                
            except Exception as e:
                print(f"  跳过第 {index+1} 行，解析错误: {e}")
                continue
                
        print(f"建立了 {len(template_index)} 个模板索引")
        return template_index
        
    except Exception as e:
        print(f"加载xlsx文件时出错: {e}")
        return {}

def parse_prediction_filename(filename):
    """解析预测PDB文件名，提取sequence_id和ec_number"""
    # 文件命名规则：sequenceID_ECnumber.pdb
    base_name = os.path.splitext(filename)[0]  # 去掉.pdb后缀
    
    # 尝试按下划线分割
    parts = base_name.split('_')
    
    if len(parts) >= 2:
        # 假设最后一部分是EC number，前面的部分组成sequence_id
        ec_number = parts[-1]
        sequence_id = '_'.join(parts[:-1])
        return sequence_id, ec_number
    else:
        print(f"    警告: 无法解析文件名格式: {filename}")
        return None, None

def scan_prediction_files(predicted_dir):
    """扫描预测结构目录，获取所有PDB文件信息"""
    if not os.path.exists(predicted_dir):
        print(f"错误: 预测结构目录不存在: {predicted_dir}")
        return []
    
    # 查找所有.pdb文件
    pdb_pattern = os.path.join(predicted_dir, "*.pdb")
    pdb_files = glob.glob(pdb_pattern)
    
    print(f"在 {predicted_dir} 中找到 {len(pdb_files)} 个PDB文件")
    
    prediction_data = []
    for pdb_file in pdb_files:
        filename = os.path.basename(pdb_file)
        sequence_id, ec_number = parse_prediction_filename(filename)
        
        if sequence_id and ec_number:
            prediction_data.append({
                'pdb_file': pdb_file,
                'filename': filename,
                'sequence_id': sequence_id,
                'ec_number': ec_number
            })
            print(f"  解析: {filename} -> SeqID={sequence_id}, EC={ec_number}")
        else:
            print(f"  跳过: {filename}")
    
    print(f"成功解析 {len(prediction_data)} 个预测文件")
    return prediction_data

def find_template_info(prediction_data, template_index):
    """根据预测文件信息查找对应的模板信息"""
    matched_pairs = []
    unmatched_predictions = []
    
    for pred in prediction_data:
        key = (pred['sequence_id'], pred['ec_number'])
        
        if key in template_index:
            template_info = template_index[key]
            matched_pair = {
                'prediction': pred,
                'template': template_info
            }
            matched_pairs.append(matched_pair)
            print(f"  匹配: {pred['filename']} -> {template_info['pdb_id']}_{template_info['chain_id']}")
        else:
            unmatched_predictions.append(pred)
            print(f"  未匹配: {pred['filename']} (SeqID={pred['sequence_id']}, EC={pred['ec_number']})")
    
    print(f"\n匹配结果:")
    print(f"成功匹配: {len(matched_pairs)} 对")
    print(f"未匹配的预测文件: {len(unmatched_predictions)} 个")
    
    return matched_pairs, unmatched_predictions

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

def set_high_quality_render():
    """设置高质量渲染参数"""
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_shadows', 1)
    cmd.set('ambient', 0.2)
    cmd.set('direct', 0.8)
    cmd.set('reflect', 0.5)
    cmd.set('shininess', 50)
    cmd.set('spec_reflect', 0.8)
    cmd.set('antialias', 2)
    cmd.set('ray_opaque_background', 1)
    
    # 关闭可能产生箭头的设置
    cmd.set('cartoon_fancy_helices', 0)
    cmd.set('cartoon_fancy_sheets', 0)
    cmd.set('cartoon_flat_sheets', 1)
    cmd.set('cartoon_smooth_loops', 1)
    cmd.set('cartoon_highlight_color', 'grey50')
    cmd.set('cartoon_putty_transform', 0)
    cmd.set('cartoon_ring_mode', 0)
    cmd.set('cartoon_tube_radius', 0.4)
    cmd.set('cartoon_oval_length', 1.2)
    cmd.set('cartoon_oval_width', 0.3)

def create_view1_overlap_comparison(query_name, template_name, output_prefix):
    """视图1: 结构重合对比图"""
    print("    生成视图1: 结构重合对比")
    
    # 清除所有显示
    cmd.hide('everything')
    
    # 显示蛋白质主链
    cmd.show('cartoon', 'all')
    cmd.color('cyan', query_name)
    cmd.color('orange', template_name)
    
    # 设置视图
    cmd.zoom('all')
    cmd.orient()
    
    # 设置渲染质量
    set_high_quality_render()
    cmd.bg_color('white')
    
    # 保存图像
    png_file = f'{output_prefix}_view1_overlap.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    print(f"    保存: {png_file}")
    
    return png_file

def create_view2_with_atp(query_name, template_name, output_prefix):
    """视图2: 蛋白质链 + ATP分子（完全按照原始代码逻辑）"""
    print("    生成视图2: 蛋白质链 + ATP")
    
    # 清除所有显示
    cmd.hide('everything')
    
    # 显示蛋白质主链
    cmd.show('cartoon', 'all')
    cmd.color('cyan', query_name)
    cmd.color('orange', template_name)
    
    # 查找并显示ATP分子（完全按照原始代码逻辑）
    atp_selection = 'resn ATP or resn ADP or resn AMP'
    nucleotide_found = False
    
    if cmd.count_atoms(atp_selection) > 0:
        cmd.show('sticks', atp_selection)
        cmd.color('red', atp_selection)
        cmd.show('spheres', f'{atp_selection} and name P*')
        cmd.set('sphere_scale', 0.3, f'{atp_selection} and name P*')
        print("    发现ATP/ADP/AMP分子")
        nucleotide_found = True
        
        # 显示ATP周围的活性位点残基（4Å范围内），按链分色
        active_site = f'byres ({atp_selection} around 4.0) and polymer.protein'

        query_sidechains    = f'({active_site} and sidechain and {query_name}) and not ({atp_selection})'
        template_sidechains = f'({active_site} and sidechain and {template_name}) and not ({atp_selection})'

        cmd.show('sticks', query_sidechains)
        cmd.show('sticks', template_sidechains)
        cmd.set('stick_radius', 0.2, f'{active_site} and sidechain')

        cmd.color('forest', query_sidechains)     # Query 链侧链
        cmd.color('gray50', template_sidechains)  # Template 链侧链
        
        print("    显示活性位点残基")
    else:
        # 如果没有ATP，尝试显示其他核苷酸
        other_nucleotides = 'resn GTP or resn GDP or resn GMP or resn CTP or resn CDP or resn CMP or resn UTP or resn UDP or resn UMP'
        if cmd.count_atoms(other_nucleotides) > 0:
            cmd.show('sticks', other_nucleotides)
            cmd.color('red', other_nucleotides)
            cmd.show('spheres', f'{other_nucleotides} and name P*')
            cmd.set('sphere_scale', 0.3, f'{other_nucleotides} and name P*')
            atp_selection = other_nucleotides
            nucleotide_found = True
            
            # 显示其他核苷酸周围的活性位点残基
            active_site = f'byres ({other_nucleotides} around 4)'
            cmd.show('sticks', f'{active_site} and sidechain')
            cmd.color('yellow', f'{active_site} and sidechain and not {other_nucleotides}')
            print("    发现其他核苷酸分子")
        else:
            print("    未发现核苷酸分子")
    
    # 设置视图
    cmd.zoom('all')
    cmd.orient()
    
    # 设置渲染质量
    set_high_quality_render()
    cmd.bg_color('white')
    
    # 保存图像
    png_file = f'{output_prefix}_view2_with_atp.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    print(f"    保存: {png_file}")
    
    return png_file, nucleotide_found, atp_selection if nucleotide_found else None

def create_view3_atp_with_4a_residues(query_name, template_name, output_prefix, atp_selection):
    """视图3: ATP + 蛋白质链 + ATP周围4Å氨基酸（完全按照原始代码逻辑）+ 保存PSE文件"""
    print("    生成视图3: ATP + 蛋白质链 + 4Å氨基酸 + PSE会话文件")
    
    if not atp_selection:
        print("    无ATP分子，无法生成视图3")
        return None
    
    # 不清除显示，在视图2基础上添加（保持与原始代码一致）
    # cmd.hide('everything')  # 注释掉，保持视图2的基础显示
    
    # 确保蛋白质链和ATP都显示（按原始代码逻辑）
    cmd.show('cartoon', 'all')
    cmd.color('cyan', query_name)
    cmd.color('orange', template_name)
    
    # 确保ATP显示（与原始代码完全一致）
    cmd.show('sticks', atp_selection)
    cmd.color('red', atp_selection)
    cmd.show('spheres', f'{atp_selection} and name P*')
    cmd.set('sphere_scale', 0.3, f'{atp_selection} and name P*')
    
    # 显示ATP周围4Å的氨基酸侧链（完全按照原始代码）
    active_site_4a = f'byres ({atp_selection} around 4.0) and polymer.protein'
    
    if cmd.count_atoms(active_site_4a) > 0:
        # 4Å范围内的蛋白侧链，按链分色
        query_sidechains    = f'({active_site_4a} and sidechain and {query_name}) and not ({atp_selection})'
        template_sidechains = f'({active_site_4a} and sidechain and {template_name}) and not ({atp_selection})'

        cmd.show('sticks', query_sidechains)
        cmd.show('sticks', template_sidechains)
        cmd.set('stick_radius', 0.2, f'{active_site_4a} and sidechain')

        cmd.color('forest',  query_sidechains)     # Query 链侧链
        cmd.color('gray50',  template_sidechains)  # Template 链侧链
        
        active_site_residues = cmd.count_atoms(f'{active_site_4a} and name CA')
        query_active_residues = cmd.count_atoms(f'{active_site_4a} and {query_name} and name CA')
        template_active_residues = cmd.count_atoms(f'{active_site_4a} and {template_name} and name CA')
        
        print(f"    添加4Å氨基酸侧链: {active_site_residues} 个残基")
        print(f"      Query链: {query_active_residues} 个残基")
        print(f"      Template链: {template_active_residues} 个残基")
    else:
        print("    ATP周围4Å内未找到氨基酸")
        active_site_residues = 0
        query_active_residues = 0
        template_active_residues = 0
    
    # 按照原始代码：设置视图聚焦ATP和活性位点
    focus_selection = f'{atp_selection} or {active_site_4a}'
    
    if cmd.count_atoms(focus_selection) > 0:
        cmd.zoom(focus_selection, buffer=2.0)  # 按原始代码的视图设置
        cmd.orient(focus_selection)
        
        # 微调视角（按原始代码）
        cmd.turn('x', 15)   # 向下倾斜
        cmd.turn('y', 25)   # 侧面角度
        cmd.turn('z', -10)  # 微调旋转
        
        print("    聚焦ATP和4Å氨基酸区域")
    else:
        cmd.zoom('all')
        cmd.orient()
    
    # 设置渲染质量
    set_high_quality_render()
    cmd.bg_color('white')
    
    # 保存PNG图像
    png_file = f'{output_prefix}_view3_atp_4a_residues.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    print(f"    保存PNG: {png_file}")
    
    # 保存PSE会话文件
    pse_file = f'{output_prefix}_view3_atp_4a_residues.pse'
    cmd.save(pse_file)
    print(f"    保存PSE会话文件: {pse_file}")
    
    return {
        'png_file': png_file,
        'pse_file': pse_file,
        'active_site_residues': active_site_residues,
        'query_active_residues': query_active_residues,
        'template_active_residues': template_active_residues,
        'hydrogen_bonds': False,  # 原始代码没有氢键，简化处理
        'atp_atoms': cmd.count_atoms(atp_selection)
    }

def calculate_overlap_metrics(query_name, template_name):
    """计算结构重叠程度的多种指标"""
    metrics = {}
    
    try:
        # 1. 基本对齐信息
        alignment_result = cmd.align(query_name, template_name)
        metrics['rmsd'] = alignment_result[0]
        metrics['aligned_atoms'] = alignment_result[1]
        metrics['alignment_cycles'] = alignment_result[2]
        
        # 2. 获取结构信息
        query_atoms = cmd.count_atoms(f'{query_name} and name CA')
        template_atoms = cmd.count_atoms(f'{template_name} and name CA')
        
        metrics['query_ca_atoms'] = query_atoms
        metrics['template_ca_atoms'] = template_atoms
        
        # 3. 计算序列覆盖度
        metrics['sequence_coverage_query'] = (metrics['aligned_atoms'] / query_atoms) * 100 if query_atoms > 0 else 0
        metrics['sequence_coverage_template'] = (metrics['aligned_atoms'] / template_atoms) * 100 if template_atoms > 0 else 0
        metrics['avg_sequence_coverage'] = (metrics['sequence_coverage_query'] + metrics['sequence_coverage_template']) / 2
        
        # 4. 计算TM-score (简化版本)
        Lq = query_atoms
        if Lq > 15:
            d0 = 1.24 * ((Lq - 15) ** (1/3)) - 1.8
        else:
            d0 = 0.5
        
        avg_distance = metrics['rmsd']
        tm_score_term = 1 / (1 + (avg_distance / d0) ** 2)
        metrics['tm_score_estimate'] = tm_score_term
        
        # 5. 计算GDT-TS估计
        rmsd = metrics['rmsd']
        
        def estimate_percentage_under_threshold(rmsd, threshold):
            if rmsd <= threshold:
                return 95
            elif rmsd <= threshold * 2:
                return max(50 - (rmsd - threshold) * 30, 10)
            else:
                return max(20 - rmsd * 2, 0)
        
        p1 = estimate_percentage_under_threshold(rmsd, 1.0)
        p2 = estimate_percentage_under_threshold(rmsd, 2.0)
        p4 = estimate_percentage_under_threshold(rmsd, 4.0)
        p8 = estimate_percentage_under_threshold(rmsd, 8.0)
        
        metrics['gdt_ts_estimate'] = (p1 + p2 + p4 + p8) / 4
        
        # 6. 计算结构相似性等级
        def classify_similarity(rmsd, coverage):
            if rmsd <= 1.0 and coverage >= 80:
                return "Very High"
            elif rmsd <= 2.0 and coverage >= 70:
                return "High"
            elif rmsd <= 3.0 and coverage >= 60:
                return "Medium"
            elif rmsd <= 5.0 and coverage >= 40:
                return "Low"
            else:
                return "Very Low"
        
        metrics['similarity_class'] = classify_similarity(rmsd, metrics['avg_sequence_coverage'])
        
        # 7. 计算重叠置信度评分 (0-100)
        rmsd_score = max(0, 100 - rmsd * 20)
        coverage_score = metrics['avg_sequence_coverage']
        alignment_score = min(100, (metrics['aligned_atoms'] / max(query_atoms, template_atoms)) * 100)
        
        metrics['overlap_confidence_score'] = (rmsd_score * 0.4 + coverage_score * 0.4 + alignment_score * 0.2)
        
        # 8. 计算结构相关性
        metrics['structural_correlation'] = math.exp(-rmsd / 2.0)
        
        print(f"    重叠指标计算完成:")
        print(f"      RMSD: {metrics['rmsd']:.3f} Å")
        print(f"      序列覆盖度: {metrics['avg_sequence_coverage']:.1f}%")
        print(f"      TM-score估计: {metrics['tm_score_estimate']:.3f}")
        print(f"      GDT-TS估计: {metrics['gdt_ts_estimate']:.1f}")
        print(f"      相似性等级: {metrics['similarity_class']}")
        print(f"      重叠置信度: {metrics['overlap_confidence_score']:.1f}/100")
        
    except Exception as e:
        print(f"    计算重叠指标时出错: {e}")
        # 返回默认值
        metrics = {
            'rmsd': 999.0, 'aligned_atoms': 0, 'alignment_cycles': 0,
            'query_ca_atoms': 0, 'template_ca_atoms': 0,
            'sequence_coverage_query': 0, 'sequence_coverage_template': 0, 'avg_sequence_coverage': 0,
            'tm_score_estimate': 0, 'gdt_ts_estimate': 0,
            'similarity_class': "Very Low", 'overlap_confidence_score': 0, 'structural_correlation': 0
        }
    
    return metrics

def visualize_three_views(query_pdb, template_pdb, output_prefix, sequence_id, pdb_id, chain_id, ec_number):
    """生成三个视图的结构可视化"""
    
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
    
    # 计算重叠指标（包含结构对齐）
    overlap_metrics = calculate_overlap_metrics(query_name, template_name)
    
    # 生成三个视图
    results = {}
    
    # 视图1: 结构重合对比
    view1_file = create_view1_overlap_comparison(query_name, template_name, output_prefix)
    results['view1_overlap'] = view1_file
    
    # 视图2: 蛋白质链 + ATP
    view2_file, nucleotide_found, atp_selection = create_view2_with_atp(query_name, template_name, output_prefix)
    results['view2_with_atp'] = view2_file
    results['nucleotide_found'] = nucleotide_found
    
    # 视图3: ATP + 蛋白质链 + 4Å氨基酸
    if nucleotide_found and atp_selection:
        view3_result = create_view3_atp_with_4a_residues(query_name, template_name, output_prefix, atp_selection)
        if view3_result:
            results['view3_atp_4a'] = view3_result['png_file']
            results['active_site_residues'] = view3_result['active_site_residues']
            results['query_active_residues'] = view3_result['query_active_residues']
            results['template_active_residues'] = view3_result['template_active_residues']
            results['hydrogen_bonds'] = view3_result['hydrogen_bonds']
            results['atp_atoms'] = view3_result['atp_atoms']
        else:
            results['view3_atp_4a'] = None
    else:
        results['view3_atp_4a'] = None
        print("    跳过视图3: 未发现核苷酸分子")
    
    # 保存叠合后的结构
    pdb_file = f'{output_prefix}_aligned.pdb'
    cmd.save(pdb_file, 'all')
    results['pdb_file'] = pdb_file
    
    # 合并重叠指标
    results.update(overlap_metrics)
    
    # 输出总结信息
    title_text = f"Query: {sequence_id} (EC: {ec_number}) vs Template: {pdb_id}_{chain_id}"
    print(f"    {title_text}")
    print(f"    生成的文件:")
    print(f"      视图1 (重合对比): {os.path.basename(results['view1_overlap'])}")
    print(f"      视图2 (含ATP): {os.path.basename(results['view2_with_atp'])}")
    if results['view3_atp_4a']:
        print(f"      视图3 (ATP+4Å氨基酸): {os.path.basename(results['view3_atp_4a'])}")
        print(f"      活性位点残基数: {results.get('active_site_residues', 'N/A')}")
    
    return results

def process_three_view_analysis(predicted_dir, xlsx_file_path, template_dir, output_dir, max_structures=None):
    """基于预测文件进行三视图结构分析的主处理函数"""
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 初始化PyMOL
    setup_pymol_session()
    
    # 1. 加载xlsx模板数据
    print("步骤1: 加载xlsx模板数据...")
    template_index = load_xlsx_data(xlsx_file_path)
    if not template_index:
        print("无法加载xlsx文件或文件为空")
        return []
    
    # 2. 扫描预测文件
    print("\n步骤2: 扫描预测PDB文件...")
    prediction_data = scan_prediction_files(predicted_dir)
    if not prediction_data:
        print("未找到预测PDB文件")
        return []
    
    # 3. 匹配预测文件和模板
    print("\n步骤3: 匹配预测文件和模板...")
    matched_pairs, unmatched_predictions = find_template_info(prediction_data, template_index)
    
    if not matched_pairs:
        print("没有找到匹配的预测-模板对")
        return []
    
    # 如果设置了最大处理数量，则只处理前N个
    if max_structures:
        matched_pairs = matched_pairs[:max_structures]
        print(f"\n限制处理前 {max_structures} 个结构对")
    
    # 4. 处理每个匹配对
    print(f"\n步骤4: 处理 {len(matched_pairs)} 个匹配对...")
    results = []
    successful_count = 0
    
    for i, pair in enumerate(matched_pairs, 1):
        pred_info = pair['prediction']
        template_info = pair['template']
        
        sequence_id = pred_info['sequence_id']
        ec_number = pred_info['ec_number']
        pdb_id = template_info['pdb_id']
        chain_id = template_info['chain_id']
        pdb_chain = template_info['pdb_chain']
        
        print(f"\n[{i}/{len(matched_pairs)}] 处理: {pred_info['filename']} vs {pdb_id}_{chain_id}")
        
        # 预测结构文件路径
        query_pdb = pred_info['pdb_file']
        print(f"    预测文件: {query_pdb}")
        
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
        
        print(f"    模板文件: {template_pdb}")
        
        # 创建输出文件前缀
        output_prefix = os.path.join(output_dir, f"{sequence_id}_{ec_number}_vs_{pdb_id}_{chain_id}")
        
        try:
            # 进行三视图结构可视化和分析
            vis_result = visualize_three_views(
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
                    'prediction_file': pred_info['filename'],
                    'query_pdb': query_pdb,
                    'template_pdb': template_pdb_raw,
                    **vis_result  # 包含所有三视图文件和重叠指标（包括PSE文件）
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
    
    # 5. 保存结果
    print(f"\n步骤5: 保存三视图分析结果...")
    summary_file = save_three_view_summary(results, output_dir)
    
    # 6. 生成统计报告
    generate_three_view_statistics(results)
    
    print(f"\n🎉 三视图处理完成！")
    print(f"扫描到的预测文件: {len(prediction_data)}")
    print(f"成功匹配的对: {len(matched_pairs)}")
    print(f"成功处理: {successful_count}")
    print(f"失败数量: {len(matched_pairs) - successful_count}")
    print(f"未匹配的预测文件: {len(unmatched_predictions)}")
    print(f"详细结果保存在: {summary_file}")
    
    return results

def save_three_view_summary(results, output_dir):
    """保存包含三视图信息的结果摘要"""
    summary_file = os.path.join(output_dir, 'three_view_analysis_summary.csv')
    
    if not results:
        return summary_file
    
    # 获取所有可能的字段名
    all_fields = set()
    for result in results:
        all_fields.update(result.keys())
    
    # 定义字段顺序
    field_order = [
        'sequence_id', 'ec_number', 'pdb_id', 'chain_id', 'pdb_chain', 'prediction_file',
        'rmsd', 'aligned_atoms', 'query_ca_atoms', 'template_ca_atoms',
        'sequence_coverage_query', 'sequence_coverage_template', 'avg_sequence_coverage',
        'tm_score_estimate', 'gdt_ts_estimate', 'similarity_class', 'overlap_confidence_score', 'structural_correlation',
        'nucleotide_found', 'view1_overlap', 'view2_with_atp', 'view3_atp_4a', 'pse_file',
        'active_site_residues', 'query_active_residues', 'template_active_residues', 
        'hydrogen_bonds', 'atp_atoms', 'pdb_file', 'query_pdb', 'template_pdb'
    ]
    
    # 确保所有字段都包含在内
    remaining_fields = all_fields - set(field_order)
    final_fields = field_order + sorted(remaining_fields)
    
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=final_fields)
        writer.writeheader()
        
        for result in results:
            # 格式化数值字段
            formatted_result = result.copy()
            for key, value in formatted_result.items():
                if isinstance(value, float):
                    if 'score' in key.lower() or 'coverage' in key.lower():
                        formatted_result[key] = f"{value:.1f}"
                    elif key == 'rmsd':
                        formatted_result[key] = f"{value:.3f}"
                    elif 'correlation' in key.lower():
                        formatted_result[key] = f"{value:.4f}"
                    else:
                        formatted_result[key] = f"{value:.2f}"
                elif key in ['view1_overlap', 'view2_with_atp', 'view3_atp_4a', 'pdb_file', 'pse_file']:
                    # 只保存文件名，不保存完整路径（包括PSE文件）
                    formatted_result[key] = os.path.basename(value) if value else ''
            
            writer.writerow(formatted_result)
    
    return summary_file

def generate_three_view_statistics(results):
    """生成三视图的统计报告"""
    if not results:
        return
    
    print(f"\n📊 三视图结构分析统计报告:")
    print(f"="*50)
    
    # 基本统计
    total_count = len(results)
    
    # RMSD统计
    rmsds = [r['rmsd'] for r in results if 'rmsd' in r and r['rmsd'] < 999]
    if rmsds:
        print(f"RMSD统计:")
        print(f"  范围: {min(rmsds):.3f} - {max(rmsds):.3f} Å")
        print(f"  平均值: {sum(rmsds)/len(rmsds):.3f} Å")
        print(f"  中位数: {sorted(rmsds)[len(rmsds)//2]:.3f} Å")
    
    # 覆盖度统计
    coverages = [r['avg_sequence_coverage'] for r in results if 'avg_sequence_coverage' in r]
    if coverages:
        print(f"\n序列覆盖度统计:")
        print(f"  范围: {min(coverages):.1f}% - {max(coverages):.1f}%")
        print(f"  平均值: {sum(coverages)/len(coverages):.1f}%")
    
    # 相似性等级分布
    similarity_classes = [r['similarity_class'] for r in results if 'similarity_class' in r]
    if similarity_classes:
        from collections import Counter
        class_counts = Counter(similarity_classes)
        print(f"\n相似性等级分布:")
        for class_name, count in class_counts.most_common():
            percentage = (count / len(similarity_classes)) * 100
            print(f"  {class_name}: {count} ({percentage:.1f}%)")
    
    # 三视图生成统计
    view1_count = sum(1 for r in results if r.get('view1_overlap'))
    view2_count = sum(1 for r in results if r.get('view2_with_atp'))
    view3_count = sum(1 for r in results if r.get('view3_atp_4a'))
    nucleotide_count = sum(1 for r in results if r.get('nucleotide_found', False))
    
    print(f"\n三视图生成统计:")
    print(f"  视图1 (重合对比): {view1_count}/{total_count} ({view1_count/total_count*100:.1f}%)")
    print(f"  视图2 (含ATP): {view2_count}/{total_count} ({view2_count/total_count*100:.1f}%)")
    print(f"  视图3 (ATP+4Å氨基酸): {view3_count}/{total_count} ({view3_count/total_count*100:.1f}%)")
    print(f"  包含核苷酸分子: {nucleotide_count}/{total_count} ({nucleotide_count/total_count*100:.1f}%)")
    
    # 活性位点统计
    active_sites = [r['active_site_residues'] for r in results if r.get('active_site_residues')]
    if active_sites:
        print(f"\n活性位点统计 (ATP周围4Å):")
        print(f"  氨基酸残基数范围: {min(active_sites)} - {max(active_sites)}")
        print(f"  平均氨基酸残基数: {sum(active_sites)/len(active_sites):.1f}")
    
    # 氢键统计
    hydrogen_bond_count = sum(1 for r in results if r.get('hydrogen_bonds', False))
    if hydrogen_bond_count > 0:
        print(f"\n氢键相互作用:")
        print(f"  检测到氢键的结构: {hydrogen_bond_count}/{view3_count} ({hydrogen_bond_count/view3_count*100:.1f}%)")
    
    # 高质量重叠统计
    high_quality = [r for r in results if r.get('overlap_confidence_score', 0) >= 70]
    print(f"\n高质量重叠 (置信度≥70): {len(high_quality)}/{total_count} ({len(high_quality)/total_count*100:.1f}%)")
    
    # TM-score和GDT-TS统计
    tm_scores = [r['tm_score_estimate'] for r in results if 'tm_score_estimate' in r]
    if tm_scores:
        print(f"\nTM-score估计统计:")
        print(f"  范围: {min(tm_scores):.3f} - {max(tm_scores):.3f}")
        print(f"  平均值: {sum(tm_scores)/len(tm_scores):.3f}")
    
    gdt_scores = [r['gdt_ts_estimate'] for r in results if 'gdt_ts_estimate' in r]
    if gdt_scores:
        print(f"\nGDT-TS估计统计:")
        print(f"  范围: {min(gdt_scores):.1f} - {max(gdt_scores):.1f}")
        print(f"  平均值: {sum(gdt_scores)/len(gdt_scores):.1f}")

def main():
    """主函数"""
    # 设置路径
    base_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code"
    
    predicted_dir = os.path.join(base_dir, "yes_pdb_result")
    xlsx_file_path = os.path.join(base_dir, "seq-supfam", "NTP_Analysis_Report.xlsx")
    template_dir = os.path.join(base_dir, "pdb_files")
    output_dir = os.path.join(base_dir, "three_view_analysis_output")
    
    # 检查文件和目录是否存在
    if not os.path.exists(predicted_dir):
        print(f"错误: 预测结构目录不存在: {predicted_dir}")
        return
    
    if not os.path.exists(xlsx_file_path):
        print(f"错误: xlsx文件不存在: {xlsx_file_path}")
        return
    
    if not os.path.exists(template_dir):
        print(f"错误: 模板目录不存在: {template_dir}")
        return
    
    print("基于预测PDB文件进行三视图结构分析")
    print("生成三个视图:")
    print("  1. 结构重合对比图（PNG）")
    print("  2. 蛋白质链 + ATP图（PNG）")
    print("  3. ATP + 蛋白质链 + ATP周围4Å氨基酸图（PNG + PSE会话文件）")
    print(f"预测结构目录: {predicted_dir}")
    print(f"xlsx模板文件: {xlsx_file_path}")
    print(f"模板结构目录: {template_dir}")
    print(f"输出目录: {output_dir}")
    
    # 询问是否限制处理数量（用于测试）
    response = input("\n是否限制处理数量？(输入数字限制，回车处理所有): ").strip()
    max_structures = None
    if response.isdigit():
        max_structures = int(response)
        print(f"将只处理前 {max_structures} 个结构对")
    
    # 处理文件
    try:
        results = process_three_view_analysis(
            predicted_dir, xlsx_file_path, template_dir, output_dir, max_structures
        )
        
        print(f"\n🎯 三视图分析完成！")
        print(f"每个成功处理的结构对都生成了3个视图:")
        print(f"  - 视图1: 结构重合对比")
        print(f"  - 视图2: 蛋白质链 + ATP")
        print(f"  - 视图3: ATP + 蛋白质链 + 4Å氨基酸")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
    
    finally:
        # 退出PyMOL
        try:
            cmd.quit()
        except:
            pass

if __name__ == "__main__":
    main()