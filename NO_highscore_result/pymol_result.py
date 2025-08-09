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

def create_view1_protein_chains(query_name, template_name, output_prefix):
    """视图1: 两条蛋白质链的对比"""
    print("    生成视图1: 蛋白质链对比")
    
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
    png_file = f'{output_prefix}_view1_chains.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    print(f"    保存: {png_file}")
    
    return png_file

def create_view2_with_atp(query_name, template_name, output_prefix):
    """视图2: 只显示两条链重合部分 + 重合区域中的ATP分子"""
    print("    生成视图2: 重合区域蛋白质链 + 重合区域ATP")
    
    # 清除所有显示
    cmd.hide('everything')
    
    # 找到两条链的重合区域 (使用8Å距离定义重合)
    overlap_query = f'{query_name} and (byres ({template_name} around 8.0))'
    overlap_template = f'{template_name} and (byres ({query_name} around 8.0))'
    overlap_region = f'({overlap_query}) or ({overlap_template})'
    
    print(f"    Query链重合区域残基数: {cmd.count_atoms(f'{overlap_query} and name CA')}")
    print(f"    Template链重合区域残基数: {cmd.count_atoms(f'{overlap_template} and name CA')}")
    
    # 显示重合区域的蛋白质链
    if cmd.count_atoms(overlap_region) > 0:
        cmd.show('cartoon', overlap_region)
        cmd.color('cyan', overlap_query)
        cmd.color('orange', overlap_template)
        print("    显示重合区域蛋白质链")
    else:
        print("    警告: 未找到重合区域，显示完整蛋白质链")
        # 如果没有重合区域，显示完整链
        cmd.show('cartoon', 'all')
        cmd.color('cyan', query_name)
        cmd.color('orange', template_name)
        overlap_region = 'all'
    
    # 查找所有ATP分子
    all_atp_selection = 'resn ATP or resn ADP or resn AMP'
    if cmd.count_atoms(all_atp_selection) == 0:
        all_atp_selection = 'resn GTP or resn GDP or resn GMP or resn CTP or resn CDP or resn CMP or resn UTP or resn UDP or resn UMP'
    
    nucleotide_found = False
    
    if cmd.count_atoms(all_atp_selection) > 0:
        # 只显示重合区域中的ATP (在重合区域6Å范围内的ATP)
        overlap_atp = f'{all_atp_selection} and (({overlap_region}) around 6.0)'
        
        if cmd.count_atoms(overlap_atp) > 0:
            # 进一步筛选：选择距离重合区域最近的1-2个ATP分子
            cmd.select('temp_overlap_atp', overlap_atp)
            atp_list = []
            cmd.iterate('temp_overlap_atp and name C1\'', 'atp_list.append((chain, resi))', space={'atp_list': atp_list})
            
            if atp_list:
                # 限制显示ATP数量，只保留前2个（如果有多个的话）
                selected_atp_residues = atp_list[:2]  # 最多2个ATP分子
                
                # 构建精确的ATP选择
                atp_parts = []
                for chain, resi in selected_atp_residues:
                    atp_parts.append(f'(chain {chain} and resi {resi})')
                
                if atp_parts:
                    final_atp = f'({" or ".join(atp_parts)}) and resn ATP+ADP+AMP+GTP+GDP+GMP+CTP+CDP+CMP+UTP+UDP+UMP'
                    
                    cmd.show('sticks', final_atp)
                    cmd.color('red', final_atp)
                    cmd.show('spheres', f'{final_atp} and name P*')
                    cmd.set('sphere_scale', 0.3, f'{final_atp} and name P*')
                    nucleotide_found = True
                    atp_selection = final_atp
                    
                    print(f"    显示筛选后的ATP: {len(selected_atp_residues)} 个分子，{cmd.count_atoms(final_atp)} 个原子")
                    
                    # 统计信息
                    total_atp_molecules = len(atp_list)
                    displayed_molecules = len(selected_atp_residues)
                    hidden_molecules = total_atp_molecules - displayed_molecules
                    print(f"    ATP分子统计: 重合区域总数 {total_atp_molecules}, 显示 {displayed_molecules}, 隐藏 {hidden_molecules}")
            
            cmd.delete('temp_overlap_atp')
        else:
            print("    重合区域6Å范围内未发现ATP分子")
            # 如果6Å内没有ATP，尝试更大范围但限制数量
            extended_atp = f'{all_atp_selection} and (({overlap_region}) around 12.0)'
            if cmd.count_atoms(extended_atp) > 0:
                cmd.select('temp_extended_atp', extended_atp)
                atp_list = []
                cmd.iterate('temp_extended_atp and name C1\'', 'atp_list.append((chain, resi))', space={'atp_list': atp_list})
                
                if atp_list:
                    # 只显示最近的1个ATP
                    selected_atp = atp_list[:1]
                    chain, resi = selected_atp[0]
                    final_atp = f'chain {chain} and resi {resi} and resn ATP+ADP+AMP+GTP+GDP+GMP+CTP+CDP+CMP+UTP+UDP+UMP'
                    
                    cmd.show('sticks', final_atp)
                    cmd.color('red', final_atp)
                    cmd.show('spheres', f'{final_atp} and name P*')
                    cmd.set('sphere_scale', 0.3, f'{final_atp} and name P*')
                    nucleotide_found = True
                    atp_selection = final_atp
                    print(f"    显示扩展范围最近ATP: 1 个分子")
                
                cmd.delete('temp_extended_atp')
            else:
                atp_selection = ""
    else:
        print("    警告: 未发现任何核苷酸分子")
        atp_selection = ""
    
    # 设置视图 - 聚焦重合区域
    if nucleotide_found:
        focus_selection = f'{overlap_region} or {atp_selection}'
    else:
        focus_selection = overlap_region
    
    if cmd.count_atoms(focus_selection) > 0:
        cmd.zoom(focus_selection)
        cmd.orient(focus_selection)
    else:
        cmd.zoom('all')
        cmd.orient()
    
    # 设置渲染质量
    set_high_quality_render()
    cmd.bg_color('white')
    
    # 保存图像
    png_file = f'{output_prefix}_view2_overlap_with_atp.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    print(f"    保存: {png_file}")
    
    return png_file, nucleotide_found

def create_view3_atp_closeup(query_name, template_name, output_prefix):
    """视图3: 基于视图2的ATP区域直接放大 - 保持相同颜色和风格"""
    print("    生成视图3: 基于视图2的ATP区域放大")
    
    # 不要清除显示！保持视图2的所有内容和颜色
    # cmd.hide('everything') - 注释掉
    
    # 查找ATP（使用与视图2相同的逻辑）
    all_atp_selection = 'resn ATP or resn ADP or resn AMP'
    if cmd.count_atoms(all_atp_selection) == 0:
        all_atp_selection = 'resn GTP or resn GDP or resn GMP+CTP+CDP+CMP+UTP+UDP+UMP'
    
    # 找到重合区域（与视图2完全相同）
    overlap_query = f'{query_name} and (byres ({template_name} around 8.0))'
    overlap_template = f'{template_name} and (byres ({query_name} around 8.0))'
    overlap_region = f'({overlap_query}) or ({overlap_template})'
    
    # 找到重合区域的ATP，但只选择一个
    overlap_atp = f'{all_atp_selection} and (({overlap_region}) around 6.0)'
    target_atp = ""
    
    if cmd.count_atoms(overlap_atp) > 0:
        cmd.select('temp_overlap_atp', overlap_atp)
        atp_list = []
        cmd.iterate('temp_overlap_atp and name C1\'', 'atp_list.append((chain, resi))', space={'atp_list': atp_list})
        
        if atp_list:
            # 只选择第一个ATP分子
            chain, resi = atp_list[0]
            target_atp = f'chain {chain} and resi {resi} and resn ATP+ADP+AMP+GTP+GDP+GMP+CTP+CDP+CMP+UTP+UDP+UMP'
            print(f"    选择单个ATP分子: 链{chain} 残基{resi}")
            
            # 如果有多个ATP，隐藏其他的
            if len(atp_list) > 1:
                for other_chain, other_resi in atp_list[1:]:
                    other_atp = f'chain {other_chain} and resi {other_resi} and resn ATP+ADP+AMP+GTP+GDP+GMP+CTP+CDP+CMP+UTP+UDP+UMP'
                    cmd.hide('everything', other_atp)
                print(f"    隐藏了其他 {len(atp_list)-1} 个ATP分子")
        
        cmd.delete('temp_overlap_atp')
    
    if not target_atp:
        print("    无法找到目标ATP，无法生成视图3")
        return None
    
    # 定义ATP周围4Å的氨基酸（添加侧链显示，但不改变主体颜色）
    active_site_4a = f'byres ({target_atp} around 4.0) and polymer.protein'
    
    # 在现有显示基础上，添加4Å范围内的氨基酸侧链
    if cmd.count_atoms(active_site_4a) > 0:
        # 显示4Å范围内的侧链
        cmd.show('sticks', f'{active_site_4a} and sidechain')
        
        # 侧链颜色稍微深一些，但与主链保持协调
        query_sidechains = f'{active_site_4a} and sidechain and {query_name}'
        template_sidechains = f'{active_site_4a} and sidechain and {template_name}'
        
        # 使用与主链协调的颜色
        cmd.color('forest', query_sidechains)      # 深青色（与青色主链协调）
        cmd.color('gray50', template_sidechains)   # 深橙色（与橙色主链协调）
        
        # 设置侧链stick大小
        cmd.set('stick_radius', 0.15, f'{active_site_4a} and sidechain')
        
        print(f"    添加4Å氨基酸侧链: {cmd.count_atoms(f'{active_site_4a} and name CA')} 个残基")
    
    # 添加氢键显示（保持与视图2的连续性）
    try:
        cmd.distance('atp_closeup_hbonds', f'{target_atp}', f'{active_site_4a}', mode=2, cutoff=3.2)
        cmd.hide('labels', 'atp_closeup_hbonds')
        cmd.color('yellow', 'atp_closeup_hbonds')
        cmd.set('dash_width', 2, 'atp_closeup_hbonds')
        print("    添加ATP相互作用氢键")
    except:
        print("    氢键显示失败")
    
    # 关键：超级放大，让ATP+氨基酸占据画面1/2大小
    focus_selection = f'{target_atp} or {active_site_4a}'
    
    try:
        # 超级极度聚焦 - ATP+氨基酸占据画面1/2
        cmd.zoom(target_atp, buffer=0.1)  # 极度放大ATP本身
        
        # 然后包含氨基酸，但整体仍然很大（占画面1/2）
        cmd.zoom(focus_selection, buffer=0.3)  # 超小buffer，ATP+氨基酸占画面1/2
        
        # 调整视角，确保ATP和氨基酸都清晰可见
        cmd.turn('x', 10)   # 向下倾斜，看到ATP结构
        cmd.turn('y', 20)   # 侧面角度，显示相互作用
        cmd.turn('z', -10)  # 微调旋转，优化朝向
        
        print("    超级放大：ATP+氨基酸占据画面1/2大小")
        
    except Exception as view_error:
        print(f"    视图聚焦失败: {view_error}")
        # 最极端的备用方案
        cmd.zoom(target_atp, buffer=0.05)  # 最大放大
        cmd.zoom(focus_selection, buffer=0.2)
    
    # 进一步优化显示，突出ATP和关键氨基酸
    try:
        # 让ATP更加突出
        cmd.set('stick_radius', 0.3, target_atp)  # ATP stick更粗
        cmd.set('sphere_scale', 0.6, f'{target_atp} and name P*')  # 磷酸基团更大
        
        # 氨基酸侧链也要清晰
        cmd.set('stick_radius', 0.25, f'{active_site_4a} and sidechain')  # 侧链stick更粗
        
        # 隐藏可能挡住ATP的前景蛋白质链
        # 定义ATP前方的区域（通过距离和位置判断）
        foreground_region = f'byres ({target_atp} around 8.0) and polymer.protein'
        
        # 选择性隐藏部分可能遮挡ATP的链段
        # 方法1: 隐藏ATP上方和前方的部分链段
        blocking_chains = f'{foreground_region} and not ({active_site_4a})'
        
        if cmd.count_atoms(blocking_chains) > 0:
            # 先尝试让这些链段透明
            cmd.set('cartoon_transparency', 0.8, blocking_chains)
            print(f"    设置前景链段为高透明度")
            
            # 如果仍然遮挡，可以完全隐藏一些链段
            # 选择距离ATP较远但仍在视野中的部分
            distant_blocking = f'byres ({target_atp} around 10.0) and polymer.protein and not (byres ({target_atp} around 6.0))'
            if cmd.count_atoms(distant_blocking) > 0:
                cmd.hide('cartoon', distant_blocking)
                print(f"    隐藏部分遮挡ATP的链段")
        
        print("    增强ATP和氨基酸的可视化效果，减少遮挡")
    except Exception as e:
        print(f"    优化显示时出错: {e}")
        pass
    
    # 设置渲染质量
    set_high_quality_render()
    cmd.bg_color('white')
    
    # 保存巨型放大图 - ATP为画面主导
    png_file = f'{output_prefix}_view3_atp_closeup.png'
    try:
        # 超高分辨率显示巨大的ATP细节
        cmd.png(png_file, width=4000, height=3000, dpi=300, ray=1)  # 4K分辨率
        print(f"    保存ATP巨型特写图（占画面1/2）: {png_file}")
    except Exception as save_error:
        print(f"    图像保存失败: {save_error}")
        return None
    
    # 保存PyMOL会话文件 - 可旋转查看
    session_file = f'{output_prefix}_view3_interactive.pse'
    try:
        cmd.save(session_file)
        print(f"    保存PyMOL交互式会话: {session_file}")
        print(f"    💡 提示: 在PyMOL中打开 {os.path.basename(session_file)} 可以自由旋转查看ATP结构")
    except Exception as session_error:
        print(f"    会话文件保存失败: {session_error}")
    
    # 可选：也保存其他3D格式
    try:
        # 保存VRML格式（某些3D查看器支持）
        vrml_file = f'{output_prefix}_view3_3d.wrl'
        cmd.save(vrml_file)
        print(f"    保存3D VRML文件: {vrml_file}")
    except:
        pass  # VRML保存失败不影响主要功能
    
    # 统计信息
    try:
        active_site_count = cmd.count_atoms(active_site_4a)
        atp_count = cmd.count_atoms(target_atp)
        query_active_residues = cmd.count_atoms(f'{active_site_4a} and {query_name} and name CA')
        template_active_residues = cmd.count_atoms(f'{active_site_4a} and {template_name} and name CA')
        
        print(f"    ATP放大图统计:")
        print(f"      ATP原子数: {atp_count}")
        print(f"      Query链活性位点残基 (4Å): {query_active_residues}")
        print(f"      Template链活性位点残基 (4Å): {template_active_residues}")
        
        return {
            'png_file': png_file,
            'session_file': session_file,
            'active_site_atoms': active_site_count,
            'atp_atoms': atp_count,
            'query_active_residues': query_active_residues,
            'template_active_residues': template_active_residues,
            'overlap_atp_atoms': atp_count
        }
    except Exception as stat_error:
        print(f"    统计计算失败: {stat_error}")
        return {
            'png_file': png_file,
            'session_file': session_file,
            'active_site_atoms': 0,
            'atp_atoms': 0,
            'query_active_residues': 0,
            'template_active_residues': 0,
            'overlap_atp_atoms': 0
        }

def visualize_structures_three_views(query_pdb, template_pdb, output_prefix, query_id, template_id, z_score):
    """使用PyMOL生成三个视图的结构可视化"""
    
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
    
    # 生成三个视图
    results = {}
    
    # 视图1: 蛋白质链对比
    view1_file = create_view1_protein_chains(query_name, template_name, output_prefix)
    results['view1_chains'] = view1_file
    
    # 视图2: 蛋白质链 + ATP
    view2_file, nucleotide_found = create_view2_with_atp(query_name, template_name, output_prefix)
    results['view2_with_atp'] = view2_file
    results['nucleotide_found'] = nucleotide_found
    
    # 视图3: ATP活性位点放大
    if nucleotide_found:
        view3_result = create_view3_atp_closeup(query_name, template_name, output_prefix)
        if view3_result:
            results['view3_atp_closeup'] = view3_result['png_file']
            results['active_site_atoms'] = view3_result['active_site_atoms']
            results['atp_atoms'] = view3_result['atp_atoms']
            results['query_active_residues'] = view3_result['query_active_residues']
            results['template_active_residues'] = view3_result['template_active_residues']
            results['overlap_atp_atoms'] = view3_result['overlap_atp_atoms']
        else:
            results['view3_atp_closeup'] = None
    else:
        results['view3_atp_closeup'] = None
        print("    跳过视图3: 未发现核苷酸分子")
    
    # 保存叠合后的结构
    pdb_file = f'{output_prefix}_aligned.pdb'
    cmd.save(pdb_file, 'all')
    results['pdb_file'] = pdb_file
    
    # 添加对齐信息
    results['rmsd'] = alignment_result[0]
    results['aligned_atoms'] = alignment_result[1]
    
    # 输出总结信息
    title_text = f"Query: {query_id} vs Template: {template_id} (Z-score: {z_score:.2f})"
    print(f"    {title_text}")
    print(f"    生成的文件:")
    print(f"      视图1 (蛋白质链): {os.path.basename(results['view1_chains'])}")
    print(f"      视图2 (含ATP): {os.path.basename(results['view2_with_atp'])}")
    if results['view3_atp_closeup']:
        print(f"      视图3 (ATP放大): {os.path.basename(results['view3_atp_closeup'])}")
        print(f"      活性位点原子数: {results.get('active_site_atoms', 'N/A')}")
    
    return results

def process_csv_data(csv_file_path, predicted_dir, template_dir, output_dir, max_structures=None):
    """基于CSV数据处理结构可视化 - 三视图版本"""
    
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
            # 进行三视图结构可视化
            vis_result = visualize_structures_three_views(
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
                    'view1_chains': vis_result['view1_chains'],
                    'view2_with_atp': vis_result['view2_with_atp'],
                    'view3_atp_closeup': vis_result['view3_atp_closeup'],
                    'active_site_atoms': vis_result.get('active_site_atoms', None),
                    'atp_atoms': vis_result.get('atp_atoms', None),
                    'query_active_residues': vis_result.get('query_active_residues', None),
                    'template_active_residues': vis_result.get('template_active_residues', None),
                    'overlap_atp_atoms': vis_result.get('overlap_atp_atoms', None),
                    'pdb_file': vis_result['pdb_file']
                }
                results.append(result)
                successful_count += 1
                print(f"    ✓ 成功处理 ({successful_count}/{i})")
            
        except Exception as e:
            print(f"    ✗ 处理失败: {e}")
            continue
    
    # 保存处理结果摘要
    summary_file = os.path.join(output_dir, 'visualization_summary_three_views.csv')
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Query_ID', 'Template_ID', 'Z_Score', 'Original_Filename', 'Query_PDB', 'Template_PDB', 
            'RMSD', 'Aligned_Atoms', 'Nucleotide_Found', 'View1_Chains', 'View2_With_ATP', 
            'View3_ATP_Closeup', 'Active_Site_Atoms', 'ATP_Atoms', 'Query_Active_Residues', 
            'Template_Active_Residues', 'Overlap_ATP_Atoms', 'PDB_File'
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
                os.path.basename(result['view1_chains']) if result['view1_chains'] else '',
                os.path.basename(result['view2_with_atp']) if result['view2_with_atp'] else '',
                os.path.basename(result['view3_atp_closeup']) if result['view3_atp_closeup'] else '',
                result['active_site_atoms'] if result['active_site_atoms'] else '',
                result['atp_atoms'] if result['atp_atoms'] else '',
                result['query_active_residues'] if result['query_active_residues'] else '',
                result['template_active_residues'] if result['template_active_residues'] else '',
                result['overlap_atp_atoms'] if result['overlap_atp_atoms'] else '',
                os.path.basename(result['pdb_file']) if result['pdb_file'] else ''
            ])
    
    print(f"\n🎉 三视图可视化处理完成！")
    print(f"总记录数: {len(csv_data)}")
    print(f"成功处理: {successful_count}")
    print(f"失败数量: {len(csv_data) - successful_count}")
    print(f"结果摘要保存在: {summary_file}")
    
    # 统计视图生成情况
    view1_count = sum(1 for r in results if r['view1_chains'])
    view2_count = sum(1 for r in results if r['view2_with_atp'])
    view3_count = sum(1 for r in results if r['view3_atp_closeup'])
    nucleotide_count = sum(1 for r in results if r['nucleotide_found'])
    
    print(f"\n📊 视图生成统计:")
    print(f"视图1 (蛋白质链): {view1_count}/{successful_count}")
    print(f"视图2 (含ATP): {view2_count}/{successful_count}")
    print(f"视图3 (ATP放大): {view3_count}/{successful_count}")
    print(f"包含核苷酸的结构: {nucleotide_count}/{successful_count}")
    
    return results

def main():
    """主函数"""
    # 设置路径
    base_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code"
    pipeline_dir = os.path.join(base_dir, "NO_highscore_result")
    
    csv_file_path = os.path.join(pipeline_dir, "zscore.csv")
    predicted_dir = os.path.join(pipeline_dir, "pdb_results")
    template_dir = os.path.join(pipeline_dir, "template")
    output_dir = os.path.join(pipeline_dir, "pymol_three_views_output")
    
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
    
    print("基于CSV文件进行PyMOL三视图结构可视化")
    print("生成三个视图:")
    print("  1. 蛋白质链对比")
    print("  2. 蛋白质链 + ATP")
    print("  3. ATP活性位点放大 (4Å范围)")
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
        
        print(f"\n📊 最终统计信息:")
        print(f"Z-score范围: {min(z_scores):.2f} - {max(z_scores):.2f}")
        print(f"平均Z-score: {sum(z_scores)/len(z_scores):.2f}")
        print(f"RMSD范围: {min(rmsds):.3f} - {max(rmsds):.3f} Å")
        print(f"平均RMSD: {sum(rmsds)/len(rmsds):.3f} Å")
        print(f"包含核苷酸的结构: {nucleotide_count}/{len(results)} ({nucleotide_count/len(results)*100:.1f}%)")
        
        # 活性位点统计
        active_sites = [r['active_site_atoms'] for r in results if r['active_site_atoms']]
        if active_sites:
            print(f"活性位点原子数范围: {min(active_sites)} - {max(active_sites)}")
            print(f"平均活性位点原子数: {sum(active_sites)/len(active_sites):.1f}")
    
    # 退出PyMOL
    cmd.quit()

if __name__ == "__main__":
    main()