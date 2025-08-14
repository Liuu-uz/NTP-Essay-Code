import os
import csv
import pymol
from pymol import cmd
import glob

def parse_extracted_csv(csv_file_path):
    """Parse extracted_first_rows.csv file - 修正版本：正确提取链信息"""
    results = []
    try:
        with open(csv_file_path, 'r', encoding='utf-8') as file:
            reader = csv.DictReader(file)
            for row in reader:
                # Extract query ID from Source_File (e.g.: Q67ZM7_zscores.csv -> Q67ZM7)
                query_id = row['Source_File'].replace('_zscores.csv', '')
                
                # Extract template ID from the last 5 characters of filename, take first 4
                filename = row['filename']
                # Remove .txt suffix
                filename_no_ext = filename.replace('.txt', '')
                
                if len(filename_no_ext) >= 5:
                    # Take last 5 characters
                    last_5_chars = filename_no_ext[-5:]
                    
                    # 修正：分别提取PDB ID和链ID
                    pdb_id = last_5_chars[:4]      # 前4个字符：PDB ID
                    chain_id = last_5_chars[4]     # 第5个字符：链ID
                    template_full_id = last_5_chars  # 完整的模板ID（包含链）
                    
                    print(f"  Filename: {filename} -> Last 5 chars: {last_5_chars}")
                    print(f"    PDB ID: {pdb_id}, Chain ID: {chain_id}, Full Template ID: {template_full_id}")
                    
                    results.append({
                        'query_id': query_id,
                        'template_pdb_id': pdb_id,        # PDB ID (4字符)
                        'template_chain_id': chain_id,    # 链ID (1字符)
                        'template_full_id': template_full_id,  # 完整ID (5字符)
                        'z_score': float(row['Z-score']),
                        'filename': filename
                    })
                else:
                    print(f"  Warning: Filename {filename} has less than 5 characters")
                    
        return results
    except Exception as e:
        print(f"Error parsing CSV file: {e}")
        return []

def find_query_pdb(predicted_dir, query_id):
    """Find query structure PDB file in predicted_structures directory"""
    # Try multiple possible filename formats
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

def find_template_pdb(template_dir, pdb_id):
    """Find template PDB file in template directory using PDB ID only"""
    if not os.path.exists(template_dir):
        print(f"    Template directory does not exist: {template_dir}")
        return None
    
    # Try multiple filename formats for PDB files
    possible_patterns = [
        f"{pdb_id.lower()}*.pdb",
        f"{pdb_id.upper()}*.pdb", 
        f"{pdb_id}*.pdb"
    ]
    
    for pattern in possible_patterns:
        template_pattern = os.path.join(template_dir, pattern)
        template_files = glob.glob(template_pattern)
        if template_files:
            print(f"    Found template file: {template_files[0]} (using pattern: {pattern})")
            return template_files[0]
    
    print(f"    Template file not found, tried patterns: {possible_patterns}")
    return None

def setup_pymol_session():
    """Initialize PyMOL session"""
    pymol.finish_launching(['pymol', '-qc'])
    cmd.reinitialize()

def set_high_quality_render():
    """Set high quality rendering parameters"""
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_shadows', 1)
    cmd.set('ambient', 0.2)
    cmd.set('direct', 0.8)
    cmd.set('reflect', 0.5)
    cmd.set('shininess', 50)
    cmd.set('spec_reflect', 0.8)
    cmd.set('antialias', 2)
    cmd.set('ray_opaque_background', 1)

def create_view1_with_highlighted_chains(query_name, template_name, output_prefix, 
                                       query_chain, template_chain):
    """View 1: 突出显示对齐的链"""
    print("    Generating View 1: Highlighting aligned chains")
    
    cmd.hide('everything')
    cmd.show('cartoon', 'all')  # 显示所有链
    
    if query_chain and template_chain:
        # 对齐的链用亮色
        cmd.color('cyan', f'{query_name} and chain {query_chain}')
        cmd.color('orange', f'{template_name} and chain {template_chain}')
        
        # 查询结构的其他链用灰色
        cmd.color('gray70', f'{query_name} and not chain {query_chain}')
        # 模板结构的其他链保持橙色（与对齐链相同颜色）
        cmd.color('orange', f'{template_name} and not chain {template_chain}')
        
        print(f"    Highlighted: Query chain {query_chain} (cyan), Template chain {template_chain} (orange)")
        print(f"    Query other chains shown in gray, Template other chains shown in orange")
        
        # 统计显示的链
        query_aligned_residues = cmd.count_atoms(f'{query_name} and chain {query_chain} and name CA')
        template_aligned_residues = cmd.count_atoms(f'{template_name} and chain {template_chain} and name CA')
        query_other_residues = cmd.count_atoms(f'{query_name} and not chain {query_chain} and name CA')
        template_other_residues = cmd.count_atoms(f'{template_name} and not chain {template_chain} and name CA')
        
        print(f"    Chain statistics:")
        print(f"      Query aligned chain {query_chain}: {query_aligned_residues} residues")
        print(f"      Template aligned chain {template_chain}: {template_aligned_residues} residues")
        if query_other_residues > 0:
            print(f"      Query other chains: {query_other_residues} residues (gray)")
        if template_other_residues > 0:
            print(f"      Template other chains: {template_other_residues} residues (orange)")
    else:
        # 标准着色 - 没有特定链对齐信息
        cmd.color('cyan', query_name)
        cmd.color('orange', template_name)
        print("    Standard coloring: Query (cyan), Template (orange)")
        
        # 统计所有链
        query_residues = cmd.count_atoms(f'{query_name} and name CA')
        template_residues = cmd.count_atoms(f'{template_name} and name CA')
        print(f"    Structure statistics:")
        print(f"      Query structure: {query_residues} residues")
        print(f"      Template structure: {template_residues} residues")
    
    cmd.zoom('all')
    cmd.orient()
    
    set_high_quality_render()
    cmd.bg_color('white')
    
    png_file = f'{output_prefix}_view1_chains.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    print(f"    Saved: {png_file}")
    
    return png_file

def create_view2_with_aligned_chains(query_name, template_name, output_prefix, 
                                     query_chain, template_chain):
    """View 2: 只在【模板结构的指定链】搜索并显示 ATP/类ATP 分子"""
    print("    Generating View 2: TEMPLATE-CHAIN ATP display")

    cmd.hide('everything')

    # 结构统计
    print(f"    Structure verification:")
    print(f"      Query structure '{query_name}': {cmd.count_atoms(query_name)} atoms")
    print(f"      Template structure '{template_name}': {cmd.count_atoms(template_name)} atoms")

    # 列出可用链
    def _list_chains(obj):
        chains = []
        cmd.iterate(f'{obj} and name CA', 'chains.append(chain)', space={'chains': chains})
        return sorted(list(set(chains)))

    print(f"    Chains in each structure:")
    unique_query_chains = _list_chains(query_name)
    unique_template_chains = _list_chains(template_name)
    print(f"      Query chains: {unique_query_chains}")
    print(f"      Template chains: {unique_template_chains}")

    # 基于对齐链定义 overlap 可视化区域（蛋白主链）
    if query_chain and template_chain:
        query_selection = f'{query_name} and chain {query_chain}'
        template_selection = f'{template_name} and chain {template_chain}'
        overlap_query = f'{query_selection} and (byres ({template_selection} around 8.0))'
        overlap_template = f'{template_selection} and (byres ({query_selection} around 8.0))'
        overlap_region = f'({overlap_query}) or ({overlap_template})'
    else:
        overlap_query = f'{query_name} and (byres ({template_name} around 8.0))'
        overlap_template = f'{template_name} and (byres ({query_name} around 8.0))'
        overlap_region = f'({overlap_query}) or ({overlap_template})'

    if cmd.count_atoms(overlap_region) > 0:
        cmd.show('cartoon', overlap_region)
        cmd.color('cyan', overlap_query)
        cmd.color('orange', overlap_template)
        print("    Showing overlapping protein regions")
    else:
        cmd.show('cartoon', 'all')
        cmd.color('cyan', query_name)
        cmd.color('orange', template_name)
        overlap_region = 'all'
        print("    Showing all protein chains")

    print(f"\n    ATP SEARCH - TEMPLATE STRUCTURE CHAIN {template_chain}:")
    print(f"    =====================================")

    # 允许的核苷酸三联
    nucleotide_types = ['ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', 'CTP', 'CDP', 'CMP', 'UTP', 'UDP', 'UMP']

    # 只在模板的指定链中找
    atp_selection = " or ".join([f"(resn {nt} and {template_name} and chain {template_chain})"
                                 for nt in nucleotide_types])
    verification_count = cmd.count_atoms(atp_selection)

    if verification_count == 0:
        print(f"    No nucleotides found on template chain {template_chain}.")
        # 可选回退：可在模板全结构内找（按需开启）
        # atp_selection = " or ".join([f"(resn {nt} and {template_name})" for nt in nucleotide_types])
        # verification_count = cmd.count_atoms(atp_selection)

    nucleotide_found = verification_count > 0
    print(f"      Selected atoms on template chain: {verification_count}")

    if nucleotide_found:
        # 显示核苷酸（模板侧）
        cmd.show('sticks', atp_selection)
        cmd.color('red', atp_selection)
        cmd.show('spheres', f'({atp_selection}) and name P*')
        cmd.set('sphere_scale', 0.3, f'({atp_selection}) and name P*')
        print(f"    Displayed nucleotides on TEMPLATE chain {template_chain}")

    # 聚焦视图
    focus_selection = f'{overlap_region} or ({atp_selection})' if nucleotide_found else overlap_region
    if cmd.count_atoms(focus_selection) > 0:
        cmd.zoom(focus_selection)
        cmd.orient(focus_selection)
    else:
        cmd.zoom('all'); cmd.orient()

    set_high_quality_render()
    cmd.bg_color('white')

    # 文件名同步改为 template
    png_file = f'{output_prefix}_view2_template_atp_only.png'
    cmd.png(png_file, width=1800, height=1400, dpi=300, ray=1)
    print(f"    Saved: {png_file}")

    return png_file, nucleotide_found

def create_view3_atp_closeup_aligned(query_name, template_name, output_prefix, 
                                     query_chain, template_chain):
    """View 3: 基于【模板结构指定链】的 ATP 放大视图"""
    print("    Generating View 3: ATP region magnification for TEMPLATE chain")

    # 对齐区域（用于筛选近邻）
    if query_chain and template_chain:
        query_selection = f'{query_name} and chain {query_chain}'
        template_selection = f'{template_name} and chain {template_chain}'
        overlap_query = f'{query_selection} and (byres ({template_selection} around 8.0))'
        overlap_template = f'{template_selection} and (byres ({query_selection} around 8.0))'
        overlap_region = f'({overlap_query}) or ({overlap_template})'
    else:
        overlap_query = f'{query_name} and (byres ({template_name} around 8.0))'
        overlap_template = f'{template_name} and (byres ({query_name} around 8.0))'
        overlap_region = f'({overlap_query}) or ({overlap_template})'

    # 候选核苷酸
    all_atp_selection = 'resn ATP or resn ADP or resn AMP'
    if cmd.count_atoms(all_atp_selection) == 0:
        all_atp_selection = ('resn GTP or resn GDP or resn GMP or '
                             'resn CTP or resn CDP or resn CMP or '
                             'resn UTP or resn UDP or resn UMP')

    # —— 仅在模板指定链搜索 —— #
    atp_found = False
    search_strategy = ""
    overlap_atp = ""

    if template_chain:
        templ_chain_atp = f'({all_atp_selection}) and {template_name} and chain {template_chain}'
        templ_chain_atp_n = cmd.count_atoms(templ_chain_atp)
        print(f"    ATP in template chain {template_chain}: {templ_chain_atp_n}")

        if templ_chain_atp_n > 0:
            # 优先：模板指定链上的 ATP
            overlap_atp = templ_chain_atp
            atp_found = True
            search_strategy = f"template chain {template_chain}"
        else:
            print(f"    No ATP found on template chain {template_chain}")
            # 回退：模板全结构 + 与重叠区 8Å 内
            templ_all_atp = f'({all_atp_selection}) and {template_name}'
            if cmd.count_atoms(templ_all_atp) > 0:
                overlap_atp = f'{templ_all_atp} and (({overlap_region}) around 8.0)'
                if cmd.count_atoms(overlap_atp) > 0:
                    atp_found = True
                    search_strategy = "template overlap region (8Å)"
                else:
                    # 最后回退：模板任意 ATP
                    overlap_atp = templ_all_atp
                    atp_found = True
                    search_strategy = "any template ATP"
    else:
        # 无模板链信息时：模板全结构策略
        templ_all_atp = f'({all_atp_selection}) and {template_name}'
        if cmd.count_atoms(templ_all_atp) > 0:
            overlap_atp = f'{templ_all_atp} and (({overlap_region}) around 8.0)'
            if cmd.count_atoms(overlap_atp) > 0:
                atp_found = True
                search_strategy = "template overlap region (8Å)"
            else:
                overlap_atp = templ_all_atp
                atp_found = True
                search_strategy = "any template ATP"

    if not atp_found:
        print("    Cannot find target ATP on template; aborting View 3")
        return None

    # 锁定一个 ATP 残基用于放大
    cmd.select('temp_overlap_atp', overlap_atp)
    atp_list = []
    cmd.iterate('temp_overlap_atp and name C1\'', 'atp_list.append((chain, resi))', space={'atp_list': atp_list})
    if not atp_list:
        cmd.iterate('temp_overlap_atp and name P', 'atp_list.append((chain, resi))', space={'atp_list': atp_list})
    if not atp_list:
        cmd.iterate('temp_overlap_atp', 'atp_list.append((chain, resi))', space={'atp_list': atp_list})

    if not atp_list:
        print("    Cannot identify individual ATP residue")
        cmd.delete('temp_overlap_atp')
        return None

    chain, resi = atp_list[0]
    target_atp = (f'chain {chain} and resi {resi} and '
                  f'resn ATP+ADP+AMP+GTP+GDP+GMP+CTP+CDP+CMP+UTP+UDP+UMP and {template_name}')
    print(f"    Selected TEMPLATE ATP for close-up: Chain {chain}, Residue {resi} (via {search_strategy})")

    # 隐藏其他 ATP
    if len(atp_list) > 1:
        for oc, oresi in atp_list[1:]:
            other_atp = (f'chain {oc} and resi {oresi} and '
                         f'resn ATP+ADP+AMP+GTP+GDP+GMP+CTP+CDP+CMP+UTP+UDP+UMP and {template_name}')
            cmd.hide('everything', other_atp)
        print(f"    Hidden {len(atp_list)-1} other ATP molecules for clarity")
    cmd.delete('temp_overlap_atp')

    # 4Å 活性位点（包含两侧蛋白以便对比）
    active_site_4a = f'byres ({target_atp} around 4.0) and polymer.protein'
    if cmd.count_atoms(active_site_4a) > 0:
        cmd.show('sticks', f'{active_site_4a} and sidechain')
        # 模板侧/查询侧配色
        if query_chain and template_chain:
            query_sidechains = f'{active_site_4a} and sidechain and {query_name} and chain {query_chain}'
            template_sidechains = f'{active_site_4a} and sidechain and {template_name} and chain {template_chain}'
        else:
            query_sidechains = f'{active_site_4a} and sidechain and {query_name}'
            template_sidechains = f'{active_site_4a} and sidechain and {template_name}'

        # 维持你原有风格：query=forest, template=gray50；ATP=红色
        cmd.color('forest', query_sidechains)
        cmd.color('gray50', template_sidechains)
        cmd.color('red', target_atp)
        cmd.set('stick_radius', 0.15, f'{active_site_4a} and sidechain')

        print(f"    Added 4Å amino acid sidechains: {cmd.count_atoms(f'{active_site_4a} and name CA')} residues")

    # H-bonds
    try:
        cmd.distance('atp_closeup_hbonds', f'{target_atp}', f'{active_site_4a}', mode=2, cutoff=3.2)
        cmd.hide('labels', 'atp_closeup_hbonds')
        cmd.color('yellow', 'atp_closeup_hbonds')
        cmd.set('dash_width', 2, 'atp_closeup_hbonds')
        print("    Added ATP interaction hydrogen bonds")
    except:
        print("    Hydrogen bond display failed")

    # 视图与优化
    focus_selection = f'{target_atp} or {active_site_4a}'
    try:
        cmd.zoom(target_atp, buffer=0.1)
        cmd.zoom(focus_selection, buffer=0.3)
        cmd.turn('x', 10); cmd.turn('y', 20); cmd.turn('z', -10)
        print("    Super magnification: ATP+amino acids occupy 1/2 of screen size")
    except Exception as view_error:
        print(f"    View focusing failed: {view_error}")
        cmd.zoom(target_atp, buffer=0.05)
        cmd.zoom(focus_selection, buffer=0.2)

    # 减遮挡
    try:
        cmd.set('stick_radius', 0.3, target_atp)
        cmd.set('sphere_scale', 0.6, f'{target_atp} and name P*')
        cmd.set('stick_radius', 0.25, f'{active_site_4a} and sidechain')

        foreground_region = f'byres ({target_atp} around 8.0) and polymer.protein'
        blocking_chains = f'{foreground_region} and not ({active_site_4a})'
        if cmd.count_atoms(blocking_chains) > 0:
            cmd.set('cartoon_transparency', 0.8, blocking_chains)
            distant_blocking = (f'byres ({target_atp} around 10.0) and polymer.protein '
                                f'and not (byres ({target_atp} around 6.0))')
            if cmd.count_atoms(distant_blocking) > 0:
                cmd.hide('cartoon', distant_blocking)
    except Exception as e:
        print(f"    Error during display optimization: {e}")

    set_high_quality_render()
    cmd.bg_color('white')

    png_file = f'{output_prefix}_view3_template_atp_closeup.png'
    try:
        cmd.png(png_file, width=4000, height=3000, dpi=300, ray=1)
        print(f"    Saved ATP close-up image: {png_file}")
    except Exception as save_error:
        print(f"    Image save failed: {save_error}")
        return None

    session_file = f'{output_prefix}_view3_interactive.pse'
    try:
        cmd.save(session_file)
        print(f"    Saved PyMOL interactive session: {session_file}")
    except Exception as session_error:
        print(f"    Session file save failed: {session_error}")

    # 统计
    try:
        active_site_count = cmd.count_atoms(active_site_4a)
        atp_count = cmd.count_atoms(target_atp)
        if query_chain and template_chain:
            query_active_residues = cmd.count_atoms(f'{active_site_4a} and {query_name} and chain {query_chain} and name CA')
            template_active_residues = cmd.count_atoms(f'{active_site_4a} and {template_name} and chain {template_chain} and name CA')
        else:
            query_active_residues = cmd.count_atoms(f'{active_site_4a} and {query_name} and name CA')
            template_active_residues = cmd.count_atoms(f'{active_site_4a} and {template_name} and name CA')

        print(f"    ATP magnification statistics:")
        print(f"      ATP atoms: {atp_count}")
        print(f"      Query active site residues (4Å): {query_active_residues}")
        print(f"      Template active site residues (4Å): {template_active_residues}")

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
        print(f"    Statistics calculation failed: {stat_error}")
        return {
            'png_file': png_file,
            'session_file': session_file,
            'active_site_atoms': 0,
            'atp_atoms': 0,
            'query_active_residues': 0,
            'template_active_residues': 0,
            'overlap_atp_atoms': 0
        }

def visualize_structures_with_chain_alignment(query_pdb, template_pdb, template_chain_id, 
                                            output_prefix, query_id, template_full_id, z_score):
    """使用指定链进行结构对齐和可视化"""
    
    # Clear current session
    cmd.delete('all')
    
    # Load structure files
    query_name = 'query_structure'
    template_name = 'template_structure'
    
    try:
        cmd.load(query_pdb, query_name)
        cmd.load(template_pdb, template_name)
    except Exception as e:
        print(f"    Failed to load PDB files: {e}")
        return None
    
    # 检查模板结构中是否存在指定的链
    template_chain_selection = f'{template_name} and chain {template_chain_id}'
    template_chain_atoms = cmd.count_atoms(f'{template_chain_selection} and name CA')
    
    if template_chain_atoms == 0:
        print(f"    Warning: Template chain {template_chain_id} not found in {template_pdb}")
        print(f"    Available chains in template:")
        
        # 检查模板中实际存在的链
        chains = []
        cmd.iterate(f'{template_name} and name CA', 'chains.append(chain)', space={'chains': chains})
        available_chains = list(set(chains))
        
        for chain in available_chains:
            chain_length = cmd.count_atoms(f'{template_name} and chain {chain} and name CA')
            print(f"      Chain {chain}: {chain_length} residues")
        
        if available_chains:
            # 使用第一个可用的链
            template_chain_id = available_chains[0]
            template_chain_selection = f'{template_name} and chain {template_chain_id}'
            print(f"    Using first available chain: {template_chain_id}")
        else:
            print(f"    No chains found in template, using whole structure")
            template_chain_selection = template_name
    else:
        print(f"    Template chain {template_chain_id} found: {template_chain_atoms} residues")
    
    # 找到query结构中最佳匹配的链
    query_chains = []
    cmd.iterate(f'{query_name} and name CA', 'query_chains.append(chain)', space={'query_chains': query_chains})
    available_query_chains = list(set(query_chains))
    
    best_alignment = None
    best_rmsd = float('inf')
    best_query_chain = None
    
    print(f"    Trying to align query chains to template chain {template_chain_id}:")
    
    for q_chain in available_query_chains:
        query_chain_selection = f'{query_name} and chain {q_chain}'
        query_chain_atoms = cmd.count_atoms(f'{query_chain_selection} and name CA')
        
        try:
            # 尝试对齐到指定的模板链
            alignment_result = cmd.align(query_chain_selection, template_chain_selection)
            rmsd = alignment_result[0]
            aligned_atoms = alignment_result[1]
            
            print(f"      Query chain {q_chain} ({query_chain_atoms} residues): RMSD={rmsd:.3f}Å, aligned_atoms={aligned_atoms}")
            
            if rmsd < best_rmsd and aligned_atoms > 20:  # 至少对齐20个原子
                best_rmsd = rmsd
                best_alignment = alignment_result
                best_query_chain = q_chain
                
        except Exception as e:
            print(f"      Query chain {q_chain} alignment failed: {e}")
            continue
    
    if best_alignment and best_query_chain:
        print(f"    Best alignment: Query chain {best_query_chain} vs Template chain {template_chain_id}")
        print(f"    Final alignment: RMSD = {best_rmsd:.3f}Å, Aligned atoms = {best_alignment[1]}")
        
        # 重新执行最佳对齐
        query_chain_selection = f'{query_name} and chain {best_query_chain}'
        cmd.align(query_chain_selection, template_chain_selection)
        
    else:
        print(f"    Chain-specific alignment failed, using whole structure alignment")
        try:
            best_alignment = cmd.align(query_name, template_name)
            best_query_chain = None
            print(f"    Whole structure alignment: RMSD = {best_alignment[0]:.3f}Å")
        except Exception as e:
            print(f"    All alignment attempts failed: {e}")
            return None
    
    # 生成可视化结果
    results = {}
    
    # View 1: 突出显示对齐的链
    view1_file = create_view1_with_highlighted_chains(
        query_name, template_name, output_prefix,
        best_query_chain, template_chain_id
    )
    results['view1_chains'] = view1_file
    
    # View 2: 蛋白质链 + ATP（基于对齐的链）
    view2_file, nucleotide_found = create_view2_with_aligned_chains(
        query_name, template_name, output_prefix,
        best_query_chain, template_chain_id
    )
    results['view2_with_atp'] = view2_file
    results['nucleotide_found'] = nucleotide_found
    
    # View 3: ATP放大视图
    if nucleotide_found:
        view3_result = create_view3_atp_closeup_aligned(
            query_name, template_name, output_prefix,
            best_query_chain, template_chain_id
        )
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
        print("    Skipping View 3: No nucleotide molecules found")
    
    # 保存对齐后的结构
    pdb_file = f'{output_prefix}_aligned.pdb'
    cmd.save(pdb_file, 'all')
    results['pdb_file'] = pdb_file
    
    # 添加对齐信息
    results['rmsd'] = best_alignment[0] if best_alignment else None
    results['aligned_atoms'] = best_alignment[1] if best_alignment else None
    results['query_chain'] = best_query_chain
    results['template_chain'] = template_chain_id
    
    return results

def process_csv_data_with_chain_info(csv_file_path, predicted_dir, template_dir, output_dir, max_structures=None):
    """处理CSV数据 - 修正版本：使用链信息进行对齐"""
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize PyMOL
    setup_pymol_session()
    
    # Parse CSV file with chain information
    csv_data = parse_extracted_csv(csv_file_path)
    if not csv_data:
        print("Cannot parse CSV file or file is empty")
        return []
    
    print(f"Parsed {len(csv_data)} records from CSV file")
    
    # If maximum processing number is set, only process first N
    if max_structures:
        csv_data = csv_data[:max_structures]
        print(f"Limited to processing first {max_structures} structures")
    
    results = []
    successful_count = 0
    
    for i, data in enumerate(csv_data, 1):
        query_id = data['query_id']
        template_pdb_id = data['template_pdb_id']
        template_chain_id = data['template_chain_id']
        template_full_id = data['template_full_id']
        z_score = data['z_score']
        filename = data['filename']
        
        print(f"\n[{i}/{len(csv_data)}] Processing: {query_id} vs {template_full_id} (Z-score: {z_score:.2f})")
        print(f"  Template: PDB={template_pdb_id}, Chain={template_chain_id}")
        print(f"  Original filename: {filename}")
        
        # Find query structure PDB file
        query_pdb = find_query_pdb(predicted_dir, query_id)
        if not query_pdb:
            print(f"    Query structure not found: {query_id}")
            continue
        
        # Find template structure PDB file (using PDB ID only)
        template_pdb = find_template_pdb(template_dir, template_pdb_id)
        if not template_pdb:
            print(f"    Template structure not found: {template_pdb_id}")
            continue
        
        print(f"    Query structure: {query_pdb}")
        print(f"    Template structure: {template_pdb}")
        
        # Create output file prefix
        output_prefix = os.path.join(output_dir, f"{query_id}_vs_{template_full_id}_zscore{z_score:.1f}")
        
        try:
            # 使用链信息进行结构对齐和可视化
            vis_result = visualize_structures_with_chain_alignment(
                query_pdb, template_pdb, template_chain_id,
                output_prefix, query_id, template_full_id, z_score
            )
            
            if vis_result:
                # Record results
                result = {
                    'query_id': query_id,
                    'template_pdb_id': template_pdb_id,
                    'template_chain_id': template_chain_id,
                    'template_full_id': template_full_id,
                    'z_score': z_score,
                    'filename': filename,
                    'query_pdb': query_pdb,
                    'template_pdb': template_pdb,
                    'rmsd': vis_result.get('rmsd'),
                    'aligned_atoms': vis_result.get('aligned_atoms'),
                    'query_chain': vis_result.get('query_chain'),
                    'template_chain': vis_result.get('template_chain'),
                    'nucleotide_found': vis_result.get('nucleotide_found', False),
                    'view1_chains': vis_result.get('view1_chains'),
                    'view2_with_atp': vis_result.get('view2_with_atp'),
                    'view3_atp_closeup': vis_result.get('view3_atp_closeup'),
                    'active_site_atoms': vis_result.get('active_site_atoms'),
                    'atp_atoms': vis_result.get('atp_atoms'),
                    'query_active_residues': vis_result.get('query_active_residues'),
                    'template_active_residues': vis_result.get('template_active_residues'),
                    'overlap_atp_atoms': vis_result.get('overlap_atp_atoms'),
                    'pdb_file': vis_result.get('pdb_file')
                }
                results.append(result)
                successful_count += 1
                print(f"    Success ({successful_count}/{i})")
            
        except Exception as e:
            print(f"    Processing failed: {e}")
            continue
    
    # Save processing results summary
    summary_file = os.path.join(output_dir, 'visualization_summary_with_chains.csv')
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Query_ID', 'Template_PDB_ID', 'Template_Chain_ID', 'Template_Full_ID', 'Z_Score', 
            'Original_Filename', 'Query_PDB', 'Template_PDB', 'RMSD', 'Aligned_Atoms', 
            'Query_Chain', 'Template_Chain', 'Nucleotide_Found', 'View1_Chains', 
            'View2_With_ATP', 'View3_ATP_Closeup', 'Active_Site_Atoms', 'ATP_Atoms',
            'Query_Active_Residues', 'Template_Active_Residues', 'Overlap_ATP_Atoms', 'PDB_File'
        ])
        
        for result in results:
            writer.writerow([
                result['query_id'],
                result['template_pdb_id'],
                result['template_chain_id'],
                result['template_full_id'],
                f"{result['z_score']:.2f}",
                result['filename'],
                result['query_pdb'],
                result['template_pdb'],
                f"{result['rmsd']:.3f}" if result['rmsd'] else '',
                result['aligned_atoms'] if result['aligned_atoms'] else '',
                result['query_chain'] if result['query_chain'] else '',
                result['template_chain'],
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
    
    print(f"\nChain-specific alignment processing complete!")
    print(f"Total records: {len(csv_data)}")
    print(f"Successfully processed: {successful_count}")
    print(f"Failed count: {len(csv_data) - successful_count}")
    print(f"Results summary saved to: {summary_file}")
    
    # Statistics on view generation
    view1_count = sum(1 for r in results if r['view1_chains'])
    view2_count = sum(1 for r in results if r['view2_with_atp'])
    view3_count = sum(1 for r in results if r['view3_atp_closeup'])
    nucleotide_count = sum(1 for r in results if r['nucleotide_found'])
    
    print(f"\nView generation statistics:")
    print(f"View 1 (protein chains): {view1_count}/{successful_count}")
    print(f"View 2 (with ATP): {view2_count}/{successful_count}")
    print(f"View 3 (ATP magnification): {view3_count}/{successful_count}")
    print(f"Structures containing nucleotides: {nucleotide_count}/{successful_count}")
    
    # Chain alignment statistics
    chain_specific_alignments = sum(1 for r in results if r['query_chain'] and r['template_chain'])
    whole_structure_alignments = successful_count - chain_specific_alignments
    
    print(f"\nAlignment statistics:")
    print(f"Chain-specific alignments: {chain_specific_alignments}/{successful_count}")
    print(f"Whole structure alignments: {whole_structure_alignments}/{successful_count}")
    
    if results:
        # RMSD statistics
        rmsds = [r['rmsd'] for r in results if r['rmsd']]
        if rmsds:
            print(f"RMSD range: {min(rmsds):.3f} - {max(rmsds):.3f} Å")
            print(f"Average RMSD: {sum(rmsds)/len(rmsds):.3f} Å")
        
        # Template chain usage statistics
        template_chains = [r['template_chain'] for r in results if r['template_chain']]
        if template_chains:
            from collections import Counter
            chain_counts = Counter(template_chains)
            print(f"Template chains used: {dict(chain_counts)}")
    
    return results

def main():
    """Main function"""
    # Set paths
    base_dir = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code"
    pipeline_dir = os.path.join(base_dir, "NO_highscore_result")
    
    csv_file_path = os.path.join(pipeline_dir, "zscore.csv")
    predicted_dir = os.path.join(pipeline_dir, "pdb_results")
    template_dir = os.path.join(pipeline_dir, "template")
    output_dir = os.path.join(pipeline_dir, "pymol_chain_aligned_output")
    
    # Check if files and directories exist
    if not os.path.exists(csv_file_path):
        print(f"Error: CSV file does not exist: {csv_file_path}")
        return
    
    if not os.path.exists(predicted_dir):
        print(f"Error: Predicted structures directory does not exist: {predicted_dir}")
        return
    
    if not os.path.exists(template_dir):
        print(f"Error: Template directory does not exist: {template_dir}")
        return
    
    print("PyMOL Chain-Specific Structure Alignment and Visualization")
    print("=" * 60)
    print("Features:")
    print("  ✓ Correct template chain ID extraction from filename")
    print("  ✓ Chain-specific structure alignment")
    print("  ✓ Intelligent fallback to whole structure alignment")
    print("  ✓ Three visualization views:")
    print("    1. Protein chain comparison (highlighting aligned chains)")
    print("    2. Overlapping regions + ATP/nucleotides")
    print("    3. ATP active site magnification (4Å range)")
    print("=" * 60)
    print(f"CSV file: {csv_file_path}")
    print(f"Predicted structures directory: {predicted_dir}")
    print(f"Template directory: {template_dir}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Ask whether to limit processing count (for testing)
    response = input("Limit processing count? (Enter number to limit, press Enter to process all): ").strip()
    max_structures = None
    if response.isdigit():
        max_structures = int(response)
        print(f"Will only process first {max_structures} structures")
    
    # Process files
    results = process_csv_data_with_chain_info(csv_file_path, predicted_dir, template_dir, output_dir, max_structures)
    
    # Display final statistics
    if results:
        print(f"\n" + "=" * 60)
        print("FINAL SUMMARY")
        print("=" * 60)
        
        z_scores = [r['z_score'] for r in results]
        rmsds = [r['rmsd'] for r in results if r['rmsd']]
        nucleotide_count = sum(1 for r in results if r['nucleotide_found'])
        chain_alignments = sum(1 for r in results if r['query_chain'] and r['template_chain'])
        
        print(f"Successfully processed structures: {len(results)}")
        print(f"Z-score range: {min(z_scores):.2f} - {max(z_scores):.2f}")
        print(f"Average Z-score: {sum(z_scores)/len(z_scores):.2f}")
        
        if rmsds:
            print(f"RMSD range: {min(rmsds):.3f} - {max(rmsds):.3f} Å")
            print(f"Average RMSD: {sum(rmsds)/len(rmsds):.3f} Å")
        
        print(f"Chain-specific alignments: {chain_alignments}/{len(results)} ({chain_alignments/len(results)*100:.1f}%)")
        print(f"Structures with nucleotides: {nucleotide_count}/{len(results)} ({nucleotide_count/len(results)*100:.1f}%)")
        
        # Active site statistics
        active_sites = [r['active_site_atoms'] for r in results if r['active_site_atoms']]
        if active_sites:
            print(f"Active site atom count range: {min(active_sites)} - {max(active_sites)}")
            print(f"Average active site atom count: {sum(active_sites)/len(active_sites):.1f}")
        
        # Template chain distribution
        template_chains = [r['template_chain'] for r in results if r['template_chain']]
        if template_chains:
            from collections import Counter
            chain_counts = Counter(template_chains)
            print(f"Template chain usage: {dict(chain_counts)}")
        
        print("\nOutput files generated:")
        print(f"  - Summary CSV: visualization_summary_with_chains.csv")
        print(f"  - PNG images: *_view1_chains.png, *_view2_overlap_with_atp.png, *_view3_atp_closeup.png")
        print(f"  - PyMOL sessions: *_view3_interactive.pse")
        print(f"  - Aligned structures: *_aligned.pdb")
    
    # Quit PyMOL
    cmd.quit()
    print("\nProcessing complete!")

if __name__ == "__main__":
    main()