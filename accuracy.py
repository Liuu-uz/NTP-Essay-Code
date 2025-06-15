from bs4 import BeautifulSoup
import os
import re

def extract_true_binding_positions(output_file_path, htm_file_path):
    """
    从output.txt文件中提取真实的结合位点
    返回: [(file_name, chain, [seq_pos_list]), ...]
    """
    # 从htm文件路径中提取目标PDB文件名和链信息
    htm_filename = os.path.basename(htm_file_path)
    pdb_name = htm_filename.split('_')[0].upper()
    target_chain = htm_filename.split('_')[1].split('.')[0]
    target_pdb_file = f"{pdb_name}.pdb"
    
    true_positions = []
    
    try:
        with open(output_file_path, 'r', encoding='utf-8') as file:
            output_data = file.read()
    except FileNotFoundError:
        print(f"错误: 文件 {output_file_path} 不存在")
        return []
    
    # 使用正则表达式提取目标PDB文件的信息
    pattern = re.compile(
        r'======= processing file：\s*' + re.escape(target_pdb_file) + r'\s*=======(.*?)=======',
        re.DOTALL
    )
    
    matches = pattern.findall(output_data)
    if not matches:
        print(f"警告: 在output.txt中未找到{target_pdb_file}的信息")
        return []
    
    pdb_section = matches[0]
    
    # 提取所有磷酸基团附近的残基
    phospho_pattern = re.compile(
        r'Neighbor residues near phospho group.*?chain\s+(\d+).*?:\s*\n(.*?)\n\n',
        re.DOTALL
    )
    
    phospho_matches = phospho_pattern.findall(pdb_section)
    
    # 链映射关系
    chain_mapping = {"0": "A", "1": "B", "2": "C", "3": "D", "4": "E", "5": "F"}
    
    for chain_id, residues_block in phospho_matches:
        mapped_chain = chain_mapping.get(chain_id, chain_id)
        if mapped_chain != target_chain:
            continue
            
        # 提取所有残基位置
        residue_positions = []
        for line in residues_block.splitlines():
            if line.startswith("SeqPos"):
                continue
            parts = line.split()
            if parts and parts[0].isdigit():
                residue_positions.append(int(parts[0]))
        
        if residue_positions:
            true_positions.append((target_pdb_file, mapped_chain, residue_positions))
    
    return true_positions

def extract_predicted_binding_positions(htm_file_path):
    """
    解析HTM文件，提取标记为B的预测结合位点
    返回: [(chain, res_num, molecule_name), ...]
    """
    with open(htm_file_path, 'r', encoding='utf-8') as file:
        soup = BeautifulSoup(file, 'html.parser')
    
    predicted_positions = []
    rows = soup.find_all('tr')
    
    # 从文件名提取链信息
    htm_filename = os.path.basename(htm_file_path)
    chain = htm_filename.split('_')[1].split('.')[0]

    # 遍历所有行，筛选出标记为B的结合点
    for row in rows[1:]:  # 跳过表头
        cols = row.find_all('td')
        if len(cols) < 11:  # 确保有足够的列
            continue
            
        try:
            res_num = int(cols[0].text.strip())  # 提取Res #位置
        except ValueError:
            continue
            
        # 检查每个小分子的结合状态
        molecules = [
            ('ATP', 2),   # ATP binding res.在第2列
            ('ADP', 4),   # ADP binding res.在第4列  
            ('AMP', 6),   # AMP binding res.在第6列
            ('GTP', 8),   # GTP binding res.在第8列
            ('GDP', 10)   # GDP binding res.在第10列
        ]
        
        for molecule_name, binding_col_idx in molecules:
            if binding_col_idx < len(cols):
                binding_status = cols[binding_col_idx].text.strip()
                if binding_status == "B":  # 只有标记为B的才是预测的结合点
                    predicted_positions.append((chain, res_num, molecule_name))
    
    return predicted_positions

def compare_predictions_with_true(predicted_positions, true_positions):
    """
    比较预测结果和真实结果
    predicted_positions: [(chain, res_num, molecule_name)]
    true_positions: [(file_name, chain, [seq_pos_list])]
    """
    # 提取所有真实的残基位置 (去重)
    true_positions_set = set()
    for _, chain, positions in true_positions:
        for pos in positions:
            true_positions_set.add((chain, pos))
    
    # 提取所有预测的残基位置 (按位置去重)
    predicted_positions_set = set()
    for chain, predicted_pos, _ in predicted_positions:
        predicted_positions_set.add((chain, predicted_pos))
    
    # 计算各种指标
    TP = 0  # 真正例：预测为结合且实际结合
    FP = 0  # 假正例：预测为结合但实际不结合
    FN = 0  # 假反例：实际结合但未预测
    
    # 计算TP：正确预测的真实位置数量
    TP = len(true_positions_set & predicted_positions_set)
    
    # 计算FP：预测但不在真实位置的数量
    FP = len(predicted_positions_set - true_positions_set)
    
    # 计算FN：真实但未被预测的数量
    FN = len(true_positions_set - predicted_positions_set)
    
    # 计算指标
    total_true_positions = len(true_positions_set)
    total_predictions = len(predicted_positions_set)  # 按位置去重后的预测数量
    
    precision = TP / total_predictions if total_predictions > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    # 找出所有未预测到的真实位点
    missed_positions = list(true_positions_set - predicted_positions_set)
    
    # 找出所有错误的预测位点
    false_positions = list(predicted_positions_set - true_positions_set)
    
    # 获取每个错误预测位点的分子列表
    false_positions_with_molecules = []
    for chain, pos in false_positions:
        molecules = [mol for c, p, mol in predicted_positions if c == chain and p == pos]
        false_positions_with_molecules.append((chain, pos, ", ".join(molecules)))
    
    print("\n" + "=" * 60)
    print("预测准确性分析结果")
    print("=" * 60)
    print(f"真实结合位点数量 (TP + FN): {total_true_positions} = {TP} (TP) + {FN} (FN)")
    print(f"预测结合位点数量 (TP + FP): {total_predictions} = {TP} (TP) + {FP} (FP)")
    print(f"真正例 (TP): {TP}")
    print(f"假正例 (FP): {FP}")
    print(f"假反例 (FN): {FN}")
    print(f"精确率 (Precision): {precision:.4f} (TP/(TP+FP))")
    print(f"召回率 (Recall): {recall:.4f} (TP/(TP+FN))")
    print(f"F1-score: {f1_score:.4f}")
    
    # 输出详细统计信息
    if missed_positions:
        print("\n未预测到的真实位点 (FN):")
        for chain, pos in missed_positions:
            print(f"  链 {chain}, 位置 {pos}")
    
    if false_positions_with_molecules:
        print("\n错误的预测位点 (FP):")
        for chain, pos, molecules in false_positions_with_molecules:
            print(f"  链 {chain}, 位置 {pos}, 分子: {molecules}")
    
    print("=" * 60)
    
    return precision, recall, f1_score

def main():
    # 文件路径
    htm_file_path = 'web_predict/1i58_A.htm'
    output_file_path = 'output.txt'
    
    print("开始蛋白质结合位点预测准确性分析")
    
    # 1. 提取真实结合位点
    true_positions = extract_true_binding_positions(output_file_path, htm_file_path)
    if not true_positions:
        print("错误: 未找到真实结合位点")
        return
    
    # 合并所有真实位置
    all_true_positions = []
    for _, chain, positions in true_positions:
        all_true_positions.extend([(chain, pos) for pos in positions])
    
    print(f"找到 {len(set(all_true_positions))} 个真实结合位点:")
    for chain, pos in set(all_true_positions):
        print(f"  链 {chain}, 位置 {pos}")
    
    # 2. 提取预测结合位点
    predicted_positions = extract_predicted_binding_positions(htm_file_path)
    if not predicted_positions:
        print("错误: 未找到预测结合位点")
        return
    
    print(f"\n找到 {len(predicted_positions)} 个预测结合位点:")
    for chain, pos, mol in predicted_positions:
        print(f"  链 {chain}, 位置 {pos}, 分子 {mol}")
    
    # 3. 比较结果
    compare_predictions_with_true(predicted_positions, true_positions)

if __name__ == "__main__":
    main()