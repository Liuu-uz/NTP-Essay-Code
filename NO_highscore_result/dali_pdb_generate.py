#!/usr/bin/env python3
"""
Extract transformation matrices from DALI HTML files and apply them to PDB structures
Usage: python dali_html_to_pdb.py <html_file> [--pdb-dir /path/to/pdb/files]
"""

import os
import sys
import re
import numpy as np
import argparse
from pathlib import Path

def parse_dali_html(html_file):
    """解析DALI HTML文件，提取比对信息和变换参数"""
    
    with open(html_file, 'r') as f:
        content = f.read()
    
    results = {}  # 使用字典避免重复
    
    print("🔍 Debug: Looking for alignment results in HTML...")
    
    # 首先查找特定格式的比对结果行
    # 例如: "1:  7tgk-D  9.9  4.4  248   593   12"
    alignment_pattern = r'(\d+):\s*([^\s]+)\s+([\d\.]+)\s+([\d\.]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+<A HREF="([^"]*)"[^>]*>PDB</A>'
    alignment_matches = re.findall(alignment_pattern, content)
    
    print(f"Found {len(alignment_matches)} alignment entries")
    
    for match in alignment_matches:
        rank = int(match[0])
        pdb_chain = match[1]
        z_score = float(match[2])
        rmsd = float(match[3])
        lali = int(match[4])
        nres = int(match[5])
        identity = int(match[6])
        pdb_url = match[7]
        
        print(f"Found alignment: {pdb_chain} Z={z_score} RMSD={rmsd}")
        
        # 提取PDB ID
        pdb_id = pdb_chain.split('-')[0] if '-' in pdb_chain else pdb_chain[:4]
        
        # 从URL中提取变换参数
        u_match = re.search(r'u=([-\d\.,]+)', pdb_url)
        t_match = re.search(r't=([-\d,]+)', pdb_url)
        
        if u_match and t_match:
            u_values = [float(x) for x in u_match.group(1).split(',')]
            t_values = [float(x) for x in t_match.group(1).split(',')]
            
            # 重构旋转矩阵 (3x3)
            rotation_matrix = np.array([
                [u_values[0], u_values[1], u_values[2]],
                [u_values[3], u_values[4], u_values[5]], 
                [u_values[6], u_values[7], u_values[8]]
            ])
            
            # 平移向量
            translation_vector = np.array(t_values)
            
            result = {
                'rank': rank,
                'pdb_chain': pdb_chain,
                'pdb_id': pdb_id,
                'z_score': z_score,
                'rmsd': rmsd,
                'lali': lali,
                'nres': nres,
                'identity': identity,
                'rotation_matrix': rotation_matrix,
                'translation_vector': translation_vector,
                'pdb_url': pdb_url,
                'source': 'alignment_line'
            }
            
            results[pdb_id] = result
            print(f"✅ Complete result for {pdb_chain}: Z={z_score}, transformations found")
    
    # 备用方法1: 查找checkbox中的变换参数
    checkbox_pattern = r'<INPUT TYPE="CHECKBOX"[^>]*VALUE="([^"]*)"[^>]*>'
    checkbox_matches = re.findall(checkbox_pattern, content, re.IGNORECASE)
    
    print(f"Found {len(checkbox_matches)} checkbox elements")
    
    for i, checkbox_value in enumerate(checkbox_matches):
        # 从checkbox_value中提取变换参数
        u_match = re.search(r'u=([-\d\.,]+)', checkbox_value)
        t_match = re.search(r't=([-\d,]+)', checkbox_value)
        
        if u_match and t_match:
            # 提取PDB信息
            cd2_match = re.search(r'cd2=([^\s]+)', checkbox_value)
            pdb_chain = cd2_match.group(1) if cd2_match else f"unknown_{i+1}"
            pdb_id = pdb_chain[:4] if len(pdb_chain) >= 4 else pdb_chain
            
            # 如果我们还没有这个结构的完整信息
            if pdb_id not in results:
                print(f"Adding from checkbox: {pdb_chain}")
                
                u_values = [float(x) for x in u_match.group(1).split(',')]
                t_values = [float(x) for x in t_match.group(1).split(',')]
                
                # 重构旋转矩阵 (3x3)
                rotation_matrix = np.array([
                    [u_values[0], u_values[1], u_values[2]],
                    [u_values[3], u_values[4], u_values[5]], 
                    [u_values[6], u_values[7], u_values[8]]
                ])
                
                # 平移向量
                translation_vector = np.array(t_values)
                
                result = {
                    'rank': len(results) + 1,
                    'pdb_chain': pdb_chain,
                    'pdb_id': pdb_id,
                    'z_score': 0.0,
                    'rmsd': 0.0,
                    'lali': 0,
                    'nres': 0,
                    'identity': 0,
                    'rotation_matrix': rotation_matrix,
                    'translation_vector': translation_vector,
                    'checkbox_value': checkbox_value,
                    'source': 'checkbox'
                }
                
                results[pdb_id] = result
    
    # 备用方法2: 查找PRE标签中的表格数据，更新已有结果
    pre_pattern = r'<PRE>(.*?)</PRE>'
    pre_matches = re.findall(pre_pattern, content, re.DOTALL | re.IGNORECASE)
    
    for pre_content in pre_matches:
        # 查找表格行: No: Chain Z rmsd lali nres %id PDB Description
        table_lines = pre_content.split('\n')
        for line in table_lines:
            # 匹配数据行格式
            match = re.match(r'\s*(\d+):\s*([^\s]+)\s+([\d\.]+)\s+([\d\.]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+.*', line)
            if match:
                rank = int(match.group(1))
                pdb_chain = match.group(2)
                z_score = float(match.group(3))
                rmsd = float(match.group(4))
                lali = int(match.group(5))
                nres = int(match.group(6))
                identity = int(match.group(7))
                
                # 提取PDB ID
                pdb_id = pdb_chain.split('-')[0] if '-' in pdb_chain else pdb_chain[:4]
                
                # 更新已有结果
                if pdb_id in results and results[pdb_id]['z_score'] == 0.0:
                    results[pdb_id].update({
                        'z_score': z_score,
                        'rmsd': rmsd,
                        'lali': lali,
                        'nres': nres,
                        'identity': identity,
                        'rank': rank
                    })
                    print(f"Updated statistics for {pdb_chain}: Z={z_score}")
    
    # 转换回列表并按rank排序
    results_list = list(results.values())
    results_list.sort(key=lambda x: x['rank'])
    
    print(f"Total unique results found: {len(results_list)}")
    
    # 显示找到的结果摘要
    for result in results_list:
        print(f"  - {result['pdb_chain']}: Z={result['z_score']:.1f}, RMSD={result['rmsd']:.1f}, Source={result['source']}")
    
    return results_list

def read_pdb_file(pdb_file, chain_id=None):
    """读取PDB文件，返回原子坐标"""
    
    atoms = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # 提取链ID
                pdb_chain = line[21].strip()
                
                # 如果指定了链ID，只读取该链
                if chain_id and pdb_chain != chain_id:
                    continue
                
                # 提取坐标
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                
                atom_info = {
                    'line': line.rstrip(),
                    'chain': pdb_chain,
                    'coords': np.array([x, y, z]),
                    'atom_name': line[12:16].strip(),
                    'residue_name': line[17:20].strip(),
                    'residue_number': int(line[22:26].strip()),
                }
                
                atoms.append(atom_info)
    
    return atoms

def apply_transformation(atoms, rotation_matrix, translation_vector):
    """对原子坐标应用旋转和平移变换"""
    
    transformed_atoms = []
    
    for atom in atoms:
        # 应用旋转矩阵和平移向量
        old_coords = atom['coords']
        new_coords = np.dot(rotation_matrix, old_coords) + translation_vector
        
        # 创建新的原子记录
        new_atom = atom.copy()
        new_atom['coords'] = new_coords
        
        transformed_atoms.append(new_atom)
    
    return transformed_atoms

def write_transformed_pdb(atoms, output_file, title="Transformed by DALI"):
    """将变换后的原子坐标写入PDB文件"""
    
    with open(output_file, 'w') as f:
        f.write(f"HEADER    {title}\n")
        f.write("REMARK   Coordinates transformed using DALI rotation matrix and translation vector\n")
        
        for atom in atoms:
            # 重构PDB行，使用新的坐标
            line = atom['line']
            new_coords = atom['coords']
            
            # 替换坐标部分
            new_line = (
                line[:30] + 
                f"{new_coords[0]:8.3f}" + 
                f"{new_coords[1]:8.3f}" + 
                f"{new_coords[2]:8.3f}" + 
                line[54:]
            )
            
            f.write(new_line + '\n')
        
        f.write("END\n")

def download_pdb_structure(pdb_id, pdb_dir):
    """下载PDB结构文件"""
    
    pdb_file = os.path.join(pdb_dir, f"{pdb_id.lower()}.pdb")
    
    if os.path.exists(pdb_file):
        print(f"PDB file already exists: {pdb_file}")
        return pdb_file
    
    print(f"Downloading PDB structure: {pdb_id}")
    
    import urllib.request
    
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    
    try:
        urllib.request.urlretrieve(url, pdb_file)
        print(f"Downloaded: {pdb_file}")
        return pdb_file
    except Exception as e:
        print(f"Failed to download {pdb_id}: {e}")
        return None

def process_dali_html(html_file, pdb_dir, output_dir, max_structures=10):
    """处理DALI HTML文件，生成变换后的PDB文件"""
    
    print(f"🔍 Parsing DALI HTML file: {html_file}")
    
    # 解析HTML文件
    results = parse_dali_html(html_file)
    
    if not results:
        print("❌ No alignment results found in HTML file")
        return
    
    print(f"📊 Found {len(results)} alignment results")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(pdb_dir, exist_ok=True)
    
    # 处理前N个最佳匹配
    for i, result in enumerate(results[:max_structures]):
        # 设置默认值，如果缺少某些信息
        pdb_chain = result.get('pdb_chain', f'unknown_{i+1}')
        z_score = result.get('z_score', 0.0)
        rmsd = result.get('rmsd', 0.0)
        lali = result.get('lali', 0)
        nres = result.get('nres', 0)
        identity = result.get('identity', 0)
        
        print(f"\n🔗 [{i+1}/{min(len(results), max_structures)}] Processing {pdb_chain}")
        if z_score > 0:
            print(f"   Z-score: {z_score:.1f}")
        if rmsd > 0:
            print(f"   RMSD: {rmsd:.1f} Å")
        
        # 提取PDB ID和链ID
        if '-' in pdb_chain:
            pdb_id, chain_id = pdb_chain.split('-')
        else:
            # 假设格式为 "1abc" + 链ID，或者只有PDB ID
            if len(pdb_chain) == 5:
                pdb_id = pdb_chain[:4]
                chain_id = pdb_chain[4]
            else:
                pdb_id = pdb_chain[:4] if len(pdb_chain) >= 4 else pdb_chain
                chain_id = None
        
        print(f"   PDB ID: {pdb_id}, Chain: {chain_id}")
        
        # 检查是否有变换矩阵
        if 'rotation_matrix' not in result or 'translation_vector' not in result:
            print(f"   ⚠️ No transformation parameters found for {pdb_chain}")
            continue
        
        # 下载PDB文件
        pdb_file = download_pdb_structure(pdb_id, pdb_dir)
        if not pdb_file:
            continue
        
        # 读取PDB文件
        print(f"   Reading PDB structure...")
        atoms = read_pdb_file(pdb_file, chain_id)
        
        if not atoms:
            print(f"   ⚠️ No atoms found for chain {chain_id} in {pdb_id}")
            # 尝试不指定链ID
            atoms = read_pdb_file(pdb_file, None)
            if atoms:
                print(f"   📋 Using all chains ({len(atoms)} atoms)")
            else:
                continue
        else:
            print(f"   Found {len(atoms)} atoms in chain {chain_id}")
        
        # 应用变换
        print(f"   Applying transformation...")
        print(f"   Rotation matrix shape: {result['rotation_matrix'].shape}")
        print(f"   Translation vector: {result['translation_vector']}")
        
        transformed_atoms = apply_transformation(
            atoms, 
            result['rotation_matrix'], 
            result['translation_vector']
        )
        
        # 生成输出文件名
        if z_score > 0:
            output_file = os.path.join(output_dir, f"transformed_{pdb_chain}_z{z_score:.1f}.pdb")
        else:
            output_file = os.path.join(output_dir, f"transformed_{pdb_chain}.pdb")
        
        # 写入变换后的PDB文件
        title_parts = [f"DALI transformation of {pdb_chain}"]
        if z_score > 0:
            title_parts.append(f"Z={z_score:.1f}")
        if rmsd > 0:
            title_parts.append(f"RMSD={rmsd:.1f}")
        
        title = ", ".join(title_parts)
        write_transformed_pdb(transformed_atoms, output_file, title)
        
        print(f"   ✅ Saved transformed structure: {output_file}")
        
        # 写入变换信息
        if z_score > 0:
            info_file = os.path.join(output_dir, f"transform_info_{pdb_chain}_z{z_score:.1f}.txt")
        else:
            info_file = os.path.join(output_dir, f"transform_info_{pdb_chain}.txt")
            
        with open(info_file, 'w') as f:
            f.write(f"DALI Transformation Information\n")
            f.write(f"==============================\n\n")
            f.write(f"Original structure: {pdb_chain}\n")
            if z_score > 0:
                f.write(f"Z-score: {z_score:.2f}\n")
            if rmsd > 0:
                f.write(f"RMSD: {rmsd:.2f} Å\n")
            if lali > 0:
                f.write(f"Aligned length: {lali} residues\n")
            if nres > 0:
                f.write(f"Total residues: {nres}\n")
            if identity > 0:
                f.write(f"Sequence identity: {identity}%\n")
            f.write(f"\nRotation matrix:\n")
            f.write(f"{result['rotation_matrix']}\n\n")
            f.write(f"Translation vector:\n")
            f.write(f"{result['translation_vector']}\n\n")
            f.write(f"Transformed PDB file: {output_file}\n")
            
            # 如果有额外的信息
            if 'checkbox_value' in result:
                f.write(f"\nOriginal checkbox value:\n")
                f.write(f"{result['checkbox_value']}\n")
        
        print(f"   📄 Saved transformation info: {info_file}")
        
    print(f"\n🎉 Processed {min(len(results), max_structures)} structures successfully!")

def main():
    parser = argparse.ArgumentParser(
        description='Extract transformations from DALI HTML and generate transformed PDB files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python dali_html_to_pdb.py a0l1A_vs_7tgkD_results.html
  python dali_html_to_pdb.py dali_results.html --pdb-dir ./pdb_files --max-structures 5
  python dali_html_to_pdb.py results.html --output-dir ./transformed_pdbs
        """
    )
    
    parser.add_argument('html_file', help='DALI HTML results file')
    parser.add_argument('--pdb-dir', default='./pdb_structures', 
                        help='Directory to store downloaded PDB files (default: ./pdb_structures)')
    parser.add_argument('--output-dir', default='./transformed_structures',
                        help='Directory to store transformed PDB files (default: ./transformed_structures)')
    parser.add_argument('--max-structures', type=int, default=10,
                        help='Maximum number of structures to transform (default: 10)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.html_file):
        print(f"❌ HTML file not found: {args.html_file}")
        sys.exit(1)
    
    print(f"🚀 Starting DALI HTML to PDB transformation...")
    print(f"📁 HTML file: {args.html_file}")
    print(f"📁 PDB directory: {args.pdb_dir}")
    print(f"📁 Output directory: {args.output_dir}")
    print(f"🔢 Max structures: {args.max_structures}")
    print("="*60)
    
    try:
        process_dali_html(args.html_file, args.pdb_dir, args.output_dir, args.max_structures)
        
        print("="*60)
        print("🎉 Transformation completed successfully!")
        print(f"📁 Transformed PDB files saved in: {args.output_dir}")
        
    except Exception as e:
        print(f"❌ Transformation failed: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()