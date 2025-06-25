from bs4 import BeautifulSoup
import os
import re

def extract_true_binding_positions(output_file_path, htm_file_path):
    """
    Extract true binding sites from output.txt file
    Returns: [(file_name, chain, [seq_pos_list]), ...]
    """
    # Extract target PDB file name and chain information from htm file path
    htm_filename = os.path.basename(htm_file_path)
    pdb_name = htm_filename.split('_')[0].upper()
    target_chain = htm_filename.split('_')[1].split('.')[0]
    target_pdb_file = f"{pdb_name}.pdb"
    
    true_positions = []
    
    try:
        with open(output_file_path, 'r', encoding='utf-8') as file:
            output_data = file.read()
    except FileNotFoundError:
        print(f"Error: File {output_file_path} does not exist")
        return []
    
    # Use regular expression to extract target PDB file information
    pattern = re.compile(
        r'======= processing file：\s*' + re.escape(target_pdb_file) + r'\s*=======(.*?)=======',
        re.DOTALL
    )
    
    matches = pattern.findall(output_data)
    if not matches:
        print(f"Warning: Information for {target_pdb_file} not found in output.txt")
        return []
    
    pdb_section = matches[0]
    
    # Extract all residues near phosphate groups
    phospho_pattern = re.compile(
        r'Neighbor residues near phospho group.*?chain\s+(\d+).*?:\s*\n(.*?)\n\n',
        re.DOTALL
    )
    
    phospho_matches = phospho_pattern.findall(pdb_section)
    
    # Chain mapping relationship
    chain_mapping = {"0": "A", "1": "B", "2": "C", "3": "D", "4": "E", "5": "F"}
    
    for chain_id, residues_block in phospho_matches:
        mapped_chain = chain_mapping.get(chain_id, chain_id)
        if mapped_chain != target_chain:
            continue
            
        # Extract all residue positions
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
    Parse HTM file and extract predicted binding sites marked as 'B'
    Returns: [(chain, res_num, molecule_name), ...]
    """
    with open(htm_file_path, 'r', encoding='utf-8') as file:
        soup = BeautifulSoup(file, 'html.parser')
    
    predicted_positions = []
    rows = soup.find_all('tr')
    
    # Extract chain information from filename
    htm_filename = os.path.basename(htm_file_path)
    chain = htm_filename.split('_')[1].split('.')[0]

    # Iterate through all rows and filter binding sites marked as 'B'
    for row in rows[1:]:  # Skip header row
        cols = row.find_all('td')
        if len(cols) < 11:  # Ensure sufficient columns
            continue
            
        try:
            res_num = int(cols[0].text.strip())  # Extract Res # position
        except ValueError:
            continue
            
        # Check binding status for each small molecule
        molecules = [
            ('ATP', 2),   # ATP binding res. in column 2
            ('ADP', 4),   # ADP binding res. in column 4  
            ('AMP', 6),   # AMP binding res. in column 6
            ('GTP', 8),   # GTP binding res. in column 8
            ('GDP', 10)   # GDP binding res. in column 10
        ]
        
        for molecule_name, binding_col_idx in molecules:
            if binding_col_idx < len(cols):
                binding_status = cols[binding_col_idx].text.strip()
                if binding_status == "B":  # Only those marked as 'B' are predicted binding sites
                    predicted_positions.append((chain, res_num, molecule_name))
    
    return predicted_positions

def compare_predictions_with_true(predicted_positions, true_positions):
    """
    Compare prediction results with true results
    predicted_positions: [(chain, res_num, molecule_name)]
    true_positions: [(file_name, chain, [seq_pos_list])]
    """
    # Extract all true residue positions (deduplicated)
    true_positions_set = set()
    for _, chain, positions in true_positions:
        for pos in positions:
            true_positions_set.add((chain, pos))
    
    # Extract all predicted residue positions (deduplicated by position)
    predicted_positions_set = set()
    for chain, predicted_pos, _ in predicted_positions:
        predicted_positions_set.add((chain, predicted_pos))
    
    # Calculate various metrics
    TP = 0  # True Positive: predicted as binding and actually binding
    FP = 0  # False Positive: predicted as binding but actually not binding
    FN = 0  # False Negative: actually binding but not predicted
    
    # Calculate TP: number of correctly predicted true positions
    TP = len(true_positions_set & predicted_positions_set)
    
    # Calculate FP: number of predictions not in true positions
    FP = len(predicted_positions_set - true_positions_set)
    
    # Calculate FN: number of true positions not predicted
    FN = len(true_positions_set - predicted_positions_set)
    
    # Calculate metrics
    total_true_positions = len(true_positions_set)
    total_predictions = len(predicted_positions_set)  # Number of predictions after deduplication by position
    
    precision = TP / total_predictions if total_predictions > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    # Find all unpredicted true binding sites
    missed_positions = list(true_positions_set - predicted_positions_set)
    
    # Find all incorrectly predicted binding sites
    false_positions = list(predicted_positions_set - true_positions_set)
    
    # Get molecule list for each incorrectly predicted binding site
    false_positions_with_molecules = []
    for chain, pos in false_positions:
        molecules = [mol for c, p, mol in predicted_positions if c == chain and p == pos]
        false_positions_with_molecules.append((chain, pos, ", ".join(molecules)))
    
    print("\n" + "=" * 60)
    print("Prediction Accuracy Analysis Results")
    print("=" * 60)
    print(f"Number of true binding sites (TP + FN): {total_true_positions} = {TP} (TP) + {FN} (FN)")
    print(f"Number of predicted binding sites (TP + FP): {total_predictions} = {TP} (TP) + {FP} (FP)")
    print(f"True Positive (TP): {TP}")
    print(f"False Positive (FP): {FP}")
    print(f"False Negative (FN): {FN}")
    print(f"Precision: {precision:.4f} (TP/(TP+FP))")
    print(f"Recall: {recall:.4f} (TP/(TP+FN))")
    print(f"F1-score: {f1_score:.4f}")
    
    # Output detailed statistics
    if missed_positions:
        print("\nUnpredicted true binding sites (FN):")
        for chain, pos in missed_positions:
            print(f"  Chain {chain}, Position {pos}")
    
    if false_positions_with_molecules:
        print("\nIncorrectly predicted binding sites (FP):")
        for chain, pos, molecules in false_positions_with_molecules:
            print(f"  Chain {chain}, Position {pos}, Molecules: {molecules}")
    
    print("=" * 60)
    
    return precision, recall, f1_score

def main():
    # File paths
    htm_file_path = 'web_predict/1i58_A.htm'
    output_file_path = 'output.txt'
    
    print("Starting protein binding site prediction accuracy analysis")
    
    # 1. Extract true binding sites
    true_positions = extract_true_binding_positions(output_file_path, htm_file_path)
    if not true_positions:
        print("Error: No true binding sites found")
        return
    
    # Merge all true positions
    all_true_positions = []
    for _, chain, positions in true_positions:
        all_true_positions.extend([(chain, pos) for pos in positions])
    
    print(f"Found {len(set(all_true_positions))} true binding sites:")
    for chain, pos in set(all_true_positions):
        print(f"  Chain {chain}, Position {pos}")
    
    # 2. Extract predicted binding sites
    predicted_positions = extract_predicted_binding_positions(htm_file_path)
    if not predicted_positions:
        print("Error: No predicted binding sites found")
        return
    
    print(f"\nFound {len(predicted_positions)} predicted binding sites:")
    for chain, pos, mol in predicted_positions:
        print(f"  Chain {chain}, Position {pos}, Molecule {mol}")
    
    # 3. Compare results
    compare_predictions_with_true(predicted_positions, true_positions)

if __name__ == "__main__":
    main()