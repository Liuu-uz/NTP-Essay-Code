import argparse
import numpy as np
import mdtraj as md
from multiprocessing import Pool
import os
import csv


def predict_atp_binding_residues(structure, target_chain_id=None):
    """
    预测ATP结合残基
    :param structure: mdtraj结构对象
    :param target_chain_id: 目标链ID（如'A'），如果为None则分析所有链
    """
    
    # ATP结合常见残基
    atp_binding_residues = {
        'K': 'high',     # 赖氨酸 - 磷酸结合
        'R': 'high',     # 精氨酸 - 磷酸结合
        'G': 'medium',   # 甘氨酸 - 柔性loop
        'S': 'medium',   # 丝氨酸 - 磷酸结合
        'T': 'medium',   # 苏氨酸 - 磷酸结合
        'N': 'medium',   # 天冬酰胺 - 腺嘌呤结合
        'D': 'medium',   # 天冬氨酸 - Mg2+配位
        'E': 'low',      # 谷氨酸 - Mg2+配位
        'H': 'low',      # 组氨酸 - 腺嘌呤结合
        'Y': 'low',      # 酪氨酸 - 芳香族相互作用
        'F': 'low',      # 苯丙氨酸 - 芳香族相互作用
        'W': 'low'       # 色氨酸 - 芳香族相互作用
    }
    
    binding_predictions = []
    
    # 获取所有链
    chains_list = list(structure.top.chains)
    
    # 如果指定了链ID，找到对应的链索引
    target_chain_idx = None
    if target_chain_id:
        # 尝试通过链ID匹配
        for chain_idx, chain in enumerate(chains_list):
            # 检查链ID是否匹配（有些文件可能有chain_id属性）
            if hasattr(chain, 'chain_id') and chain.chain_id == target_chain_id:
                target_chain_idx = chain_idx
                break
            # 如果没有chain_id属性，假设按字母顺序排列（A=0, B=1, etc.）
            elif target_chain_id.isalpha() and len(target_chain_id) == 1:
                expected_idx = ord(target_chain_id.upper()) - ord('A')
                if chain_idx == expected_idx:
                    target_chain_idx = chain_idx
                    break
        
        if target_chain_idx is None:
            print(f"警告: 未找到链 {target_chain_id}，将分析第一条蛋白质链")
            # 如果找不到指定链，使用第一条含有蛋白质的链
            for chain_idx, chain in enumerate(chains_list):
                if any(r.is_protein for r in chain.residues):
                    target_chain_idx = chain_idx
                    break
    else:
        # 如果没有指定链，使用第一条含有蛋白质的链
        for chain_idx, chain in enumerate(chains_list):
            if any(r.is_protein for r in chain.residues):
                target_chain_idx = chain_idx
                break
    
    if target_chain_idx is None:
        print("错误: 未找到蛋白质链")
        return binding_predictions
    
    # 获取目标链对象
    target_chain = chains_list[target_chain_idx]
    
    # 获取目标链的蛋白质残基 - 这里是关键修复！
    protein_residues = [r for r in target_chain.residues if r.is_protein]
    
    print(f"正在分析链 {target_chain_id or target_chain_idx}，包含 {len(protein_residues)} 个蛋白质残基")
    
    # 重要：现在我们只对目标链的残基进行编号！
    for i, residue in enumerate(protein_residues):
        residue_code = residue.name
        one_letter = get_one_letter_code(residue_code)
        
        # 基本预测
        binding_probability = 'N'  # 默认不结合
        
        if one_letter in atp_binding_residues:
            priority = atp_binding_residues[one_letter]
            
            # 检查序列环境（在当前链内）
            context_score = analyze_sequence_context(protein_residues, i, atp_binding_residues)
            
            # 结合残基类型和环境评分
            if priority == 'high' and context_score > 0.3:
                binding_probability = 'B'  # 高可能性结合
            elif priority == 'high' and context_score > 0.1:
                binding_probability = 'B'
            elif priority == 'medium' and context_score > 0.4:
                binding_probability = 'B'
            elif priority == 'low' and context_score > 0.6:
                binding_probability = 'B'
        
        binding_predictions.append({
            'residue_index': i + 1,  # 在目标链内的编号：1, 2, 3...
            'residue_code': one_letter,
            'residue_name': residue_code,
            'binding': binding_probability,
            'resSeq': residue.resSeq
        })
    
    return binding_predictions


def analyze_sequence_context(residues, position, atp_binding_residues):
    """
    分析残基的序列环境，寻找ATP结合motif
    """
    window_size = 5  # 前后各5个残基
    context_score = 0.0
    
    start = max(0, position - window_size)
    end = min(len(residues), position + window_size + 1)
    
    # 计算窗口内ATP结合残基的密度
    atp_residue_count = 0
    total_count = 0
    
    for i in range(start, end):
        if i != position:  # 不包括当前残基
            residue_code = residues[i].name
            one_letter = get_one_letter_code(residue_code)
            total_count += 1
            
            if one_letter in atp_binding_residues:
                weight = {'high': 1.0, 'medium': 0.6, 'low': 0.3}[atp_binding_residues[one_letter]]
                atp_residue_count += weight
    
    if total_count > 0:
        context_score = atp_residue_count / total_count
    
    # 检查特定motif
    # Walker A motif: GXXXXGK[S/T]
    if check_walker_a_motif(residues, position):
        context_score += 0.5
    
    # Walker B motif: hhhhD (h=疏水性残基)
    if check_walker_b_motif(residues, position):
        context_score += 0.3
    
    return context_score


def check_walker_a_motif(residues, position):
    """
    检查Walker A motif: GXXXXGK[S/T]
    """
    if position >= 7:
        sequence = []
        for i in range(position - 7, min(len(residues), position + 1)):
            sequence.append(get_one_letter_code(residues[i].name))
        
        seq_str = ''.join(sequence)
        # 简化的Walker A检查
        if 'G' in seq_str and 'K' in seq_str and ('S' in seq_str or 'T' in seq_str):
            return True
    
    return False


def check_walker_b_motif(residues, position):
    """
    检查Walker B motif: hhhhD
    """
    hydrophobic = ['A', 'V', 'L', 'I', 'M', 'F', 'W', 'Y', 'P']
    
    if position >= 4:
        sequence = []
        for i in range(position - 4, min(len(residues), position + 1)):
            sequence.append(get_one_letter_code(residues[i].name))
        
        # 检查是否有足够的疏水性残基后跟天冬氨酸
        if len(sequence) >= 5:
            hydrophobic_count = sum(1 for aa in sequence[:-1] if aa in hydrophobic)
            if hydrophobic_count >= 2 and sequence[-1] == 'D':
                return True
    
    return False


def get_one_letter_code(three_letter):
    """
    将三字母氨基酸代码转换为单字母代码
    """
    code_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return code_map.get(three_letter, 'X')


def parse_filename(filename):
    """
    解析文件名获取蛋白质ID和链ID
    例如: 1d4x_A_af.cif -> protein_id='1d4x', chain_id='A'
    """
    basename = os.path.splitext(filename)[0]  # 去掉.cif
    parts = basename.split('_')
    
    if len(parts) >= 2:
        protein_id = parts[0]
        chain_id = parts[1]
        return protein_id, chain_id
    else:
        # 如果文件名格式不符合预期，返回整个文件名作为protein_id
        return basename, None


def process_cif(file):
    """
    处理单个CIF文件，预测ATP结合位点
    """
    try:
        # 加载结构
        structure = md.load(file)
        
        # 解析文件名获取蛋白质ID和链ID
        filename = os.path.basename(file)
        protein_id, target_chain_id = parse_filename(filename)
        
        print(f"处理文件: {filename}")
        print(f"解析得到 - 蛋白质ID: {protein_id}, 链ID: {target_chain_id}")
        
        # 调试：打印所有链的信息
        chains_list = list(structure.top.chains)
        print(f"文件中共有 {len(chains_list)} 条链")
        
        for i, chain in enumerate(chains_list):
            protein_residues_in_chain = [r for r in chain.residues if r.is_protein]
            if protein_residues_in_chain:
                # 获取前10个残基的序列作为预览
                preview_seq = ''.join([get_one_letter_code(r.name) for r in protein_residues_in_chain[:10]])
                print(f"  链 {i}: {len(protein_residues_in_chain)} 个蛋白质残基, 开始序列: {preview_seq}...")
        
        # 预测ATP结合残基（只分析指定链）
        predictions = predict_atp_binding_residues(structure, target_chain_id)
        
        if not predictions:
            print(f"警告: 文件 {filename} 中没有找到可分析的残基")
            return []
        
        # 返回结果列表，每个残基一行
        results = []
        for pred in predictions:
            results.append({
                'Prot.ID': f"{protein_id}_{target_chain_id}" if target_chain_id else protein_id,
                'NO': pred['residue_index'],
                'Residue': pred['residue_code'],
                'ATP_Binding_Site': pred['binding']
            })
        
        return results
        
    except Exception as e:
        print(f"处理文件 {file} 时出错: {str(e)}")
        return []


def get_cif_files(folder):
    """
    获取文件夹中所有CIF文件的路径列表
    """
    if not os.path.exists(folder):
        raise ValueError(f"文件夹不存在: {folder}")
    
    all_files = sorted(os.listdir(folder), key=lambda x: x.lower())
    cif_files = [
        os.path.join(folder, f) 
        for f in all_files 
        if f.lower().endswith('.cif')
    ]
    return cif_files


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='预测CIF文件中的ATP结合位点')
    parser.add_argument('--folder', required=True, help='包含CIF文件的目录')
    parser.add_argument('--processes', type=int, default=4, help='并行进程数 (默认: 4)')
    parser.add_argument('--output', default='atp_binding_sites_alphafold.csv', help='输出CSV文件名')
    parser.add_argument('--format', choices=['csv', 'txt'], default='csv', help='输出格式 (csv 或 txt)')
    
    args = parser.parse_args()

    # 获取CIF文件列表
    try:
        cif_files = get_cif_files(args.folder)
        print(f"找到 {len(cif_files)} 个CIF文件")
        
        if not cif_files:
            print("警告: 在指定文件夹中未找到CIF文件")
            exit(1)
            
    except ValueError as e:
        print(f"错误: {e}")
        exit(1)

    # 并行处理
    print(f"开始分析，使用 {args.processes} 个进程...")
    
    all_results = []
    
    with Pool(processes=args.processes) as pool:
        results = pool.imap(process_cif, cif_files)
        
        for i, result in enumerate(results, 1):
            print(f"完成 {i}/{len(cif_files)}")
            all_results.extend(result)
    
    # 输出结果
    if args.format == 'csv':
        with open(args.output, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Prot.ID', 'NO', 'Residue', 'ATP_Binding_Site']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for row in all_results:
                writer.writerow(row)
    
    else:  # txt格式
        output_file = args.output.replace('.csv', '.txt')
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(f"{'Prot.ID':<15} {'NO':<5} {'Residue':<8} {'ATP Binding Site':<15}\n")
            f.write("-" * 50 + "\n")
            
            for row in all_results:
                f.write(f"{row['Prot.ID']:<15} {row['NO']:<5} {row['Residue']:<8} {row['ATP_Binding_Site']:<15}\n")
    
    print(f"分析完成！结果保存在 {args.output}")
    
    # 统计信息
    total_residues = len(all_results)
    binding_residues = len([r for r in all_results if r['ATP_Binding_Site'] == 'B'])
    
    print(f"统计信息:")
    print(f"  总残基数: {total_residues}")
    print(f"  预测ATP结合残基数: {binding_residues}")
    print(f"  ATP结合残基比例: {binding_residues/total_residues*100:.1f}%")