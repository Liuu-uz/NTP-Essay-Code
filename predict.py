import argparse
import subprocess
import re
import os
import mdtraj as md
import numpy as np
import pandas as pd

def run_process_pdb(pdb_file):
    """调用 proc.py 运行 process_pdb 并获取输出"""
    try:
        result = subprocess.run(
            ["python", "proc.py", pdb_file], 
            capture_output=True, text=True
        )
        return result.stdout.splitlines()
    except Exception as e:
        print(f"❌ 运行 proc.py 失败: {e}")
        return []

def map_residues_to_sequence(proc_output):
    """解析 `proc.py` 输出，提取蛋白质氨基酸的 Residue Index 映射到 `sequence` 里的编号"""
    sequence_residues = []
    sequence_start, sequence_end = None, None
    sequence_string = ""

    for line in proc_output:
        # ✅ 解析 `range: X Y`，确定蛋白质氨基酸的 Residue Index 范围
        match = re.search(r"range:\s*(\d+)\s+(\d+)", line)
        if match:
            sequence_start, sequence_end = int(match.group(1)), int(match.group(2))

        # ✅ 解析 `sequence (chain X ...)`，提取氨基酸序列
        match = re.search(r"sequence \(chain \w+ \d+\.\d+\):\s*([\w]+)", line)
        if match:
            sequence_string = match.group(1)

    # ✅ 重新编号氨基酸 Residue Index，从 `sequence_start` 开始
    if sequence_string and sequence_start is not None:
        sequence_residues = [(i + sequence_start, aa) for i, aa in enumerate(sequence_string)]
    
    return sequence_residues

def parse_proc_output(proc_output, target_chain):
    """解析 `proc.py` 输出，提取氨基酸、金属离子、小分子 (ATP, GTP, UDP, ADP)"""
    binding_residues = []
    metal_ions = []
    nucleotide_molecules = []  # ATP, GTP, UDP, ADP, AMP 等小分子

    AMINO_ACIDS = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", 
                   "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", 
                   "TYR", "VAL"}

    NUCLEOTIDES = {"ATP", "GTP", "ADP", "AMP", "UDP", "NAD", "GDP", "FAD"}

    chain_id_mapping = {"0": "A", "1": "B", "2": "C", "3": "D"}

    for line in proc_output:
        # ✅ 解析金属离子
        if "ions:" in line:
            matches = re.findall(r"(\w+\d+)\s*\((\d+), \d+\)", line)
            for match in matches:
                ion_name, ion_index = match
                metal_ions.append((ion_name, int(ion_index)))

        # ✅ 解析小分子 (ATP, GTP, UDP, ADP)
        if "residue:" in line:
            match = re.search(r"Residue index:\s*(\d+),\s*residue:\s*([\w\d]+)", line)
            if match:
                residue_index = int(match.group(1))
                residue_name = match.group(2).strip()
                for nucleotide in NUCLEOTIDES:
                    if residue_name.startswith(nucleotide):
                        nucleotide_molecules.append((residue_name, residue_index))
                        break  # 避免重复存入

        # ✅ 解析蛋白质氨基酸
        if "Residue index:" in line:
            match = re.search(r"Residue index:\s*(\d+),\s*residue:\s*([\w\d]+).*sequence \(chain (\w+)", line)
            if match:
                residue_index = int(match.group(1))
                residue_name = match.group(2).strip()
                chain_id = match.group(3).strip()
                chain_id = chain_id_mapping.get(chain_id, chain_id)

                if residue_name in AMINO_ACIDS and chain_id == target_chain:
                    binding_residues.append((residue_index, residue_name, chain_id))

    return binding_residues, metal_ions, nucleotide_molecules

def renumber_binding_sites(binding_residues, sequence_residues):
    """将 binding_residues (PDB 里的 Residue Index) 映射到 sequence_residues 里的编号"""
    renumbered_binding_sites = []

    residue_index_map = {pdb_index: i + 1 for i, (pdb_index, _) in enumerate(sequence_residues)}

    for pdb_index, residue_name, chain_id in binding_residues:
        if pdb_index in residue_index_map:
            new_index = residue_index_map[pdb_index]
            renumbered_binding_sites.append((new_index, residue_name, chain_id, pdb_index))

    return renumbered_binding_sites

def compute_distances(structure, residue_atoms, ligand_atoms):
    """计算氨基酸磷原子 vs. 金属离子/ATP 分子的距离"""
    pairs = [[r, l] for r in residue_atoms for l in ligand_atoms]
    distances = md.compute_distances(structure, pairs)
    return distances

def predict_atp_binding(pdb_file, target_chain):
    """调用 process_pdb 解析 PDB 文件并预测 ATP 结合位点"""
    proc_output = run_process_pdb(pdb_file)
    sequence_residues = map_residues_to_sequence(proc_output)
    binding_residues, metal_ions, nucleotide_molecules = parse_proc_output(proc_output, target_chain)

    if not binding_residues:
        print(f"⚠️ 未找到任何氨基酸（链 {target_chain}），可能解析错误！")
        return
    if not metal_ions and not nucleotide_molecules:
        print(f"⚠️ 未找到金属离子或 ATP 相关分子，可能解析错误！")

    structure = md.load(pdb_file)
    renumbered_binding_sites = renumber_binding_sites(binding_residues, sequence_residues)

    binding_sites = []
    for new_index, residue_name, chain_id, original_index in renumbered_binding_sites:
        residue_atoms = [atom.index for atom in structure.top.residue(original_index).atoms if atom.element.symbol == "P"]
        ligand_atoms = [atom.index for _, idx in metal_ions + nucleotide_molecules for atom in structure.top.residue(idx).atoms]

        distances = compute_distances(structure, residue_atoms, ligand_atoms) if ligand_atoms else np.array([])
        is_binding = "Yes" if distances.size > 0 and np.any(distances <= 4.0) else "No"

        binding_sites.append({"Residue Index": new_index, "Residue Name": residue_name, "Chain": chain_id, "Is Binding": is_binding})

    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{os.path.basename(pdb_file).split('.')[0]}_{target_chain}_atp_sites.csv")
    pd.DataFrame(binding_sites).to_csv(output_file, index=False)
    print(f"✅ 发现 {len(binding_sites)} 个 ATP 结合位点（链 {target_chain}），已保存至: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="预测 ATP 结合位点")
    parser.add_argument("pdb_file", help="PDB 文件路径")
    parser.add_argument("chain_id", help="目标链 ID（如 A, B, C, D）")
    args = parser.parse_args()
    predict_atp_binding(args.pdb_file, args.chain_id)

if __name__ == "__main__":
    main()