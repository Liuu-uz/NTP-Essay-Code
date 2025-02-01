import argparse
import mdtraj as md
import pandas as pd
import os
import subprocess

# ATP 结合残基集合
ATP_BINDING_RESIDUES = {"ARG", "LYS", "HIS", "ASP", "GLU", "SER", "THR", "TYR"}

def validate_pdb(pdb_file):
    """
    调用 validate_pdb.py 确保 PDB 结构完整
    """
    result = subprocess.run(["python", "validate_pdb.py", pdb_file], capture_output=True, text=True)
    if "❌ 解析失败" in result.stdout:
        print(f"❌ {pdb_file} 结构无效，终止 ATP 结合点预测")
        return False
    print(f"✅ PDB 文件 {pdb_file} 通过验证")
    return True

def get_proc_results(pdb_file):
    """
    调用 proc.py 解析 PDB 并获取磷酸基团位点 & 金属离子
    """
    result = subprocess.run(["python", "proc.py", pdb_file], capture_output=True, text=True)
    phosphate_sites = []
    metal_ions = []
    
    for line in result.stdout.split("\n"):
        if "Residue index" in line:
            parts = line.split(", ")
            residue_info = {p.split(": ")[0]: p.split(": ")[1] for p in parts if ": " in p}
            residue_index = int(residue_info.get("Residue index", "-1")) + 1  # ✅ 调整 Residue_Index +1
            residue_name = residue_info.get("residue", "")
            
            if residue_name.startswith("GTP") or residue_name.startswith("DTP"):
                phosphate_sites.append(residue_index)
            
            if "ions:" in line:
                metal_ions.append(residue_info["ions"].split(" ")[0])  # 获取离子名

    print(f"🔍 磷酸基团位点: {phosphate_sites}")
    print(f"🔍 金属离子: {metal_ions}")
    
    return phosphate_sites, metal_ions

def predict_atp_binding(pdb_file, chain_id, output_dir="results"):
    """
    预测 ATP 结合点：
    - 先调用 validate_pdb.py 进行 PDB 结构验证
    - 运行 proc.py 获取关键残基信息
    - 结合 `mdtraj` 预测 ATP 结合位点，并让 Residue_Index 从 1 开始
    """
    if not validate_pdb(pdb_file):
        return

    phosphate_sites, metal_ions = get_proc_results(pdb_file)

    traj = md.load(pdb_file)
    top = traj.topology
    pdb_id = os.path.basename(pdb_file).split(".")[0]

    chain = [c for c in top.chains if c.index == chain_id]
    if not chain:
        print(f"❌ Chain {chain_id} 不存在于 {pdb_file}")
        return
    
    # **核心逻辑：筛选与磷酸结合的氨基酸**
    binding_sites = []
    for res in chain[0].residues:
        real_residue_index = res.index + 1  # ✅ Residue_Index +1 让其从 1 开始
        if real_residue_index in phosphate_sites or res.name in ATP_BINDING_RESIDUES:
            binding_sites.append([pdb_id, chain_id, real_residue_index, res.name])
    
    # 结果存储
    df = pd.DataFrame(binding_sites, columns=["PDB_ID", "Chain", "Residue_Index", "Residue_Name"])
    os.makedirs(output_dir, exist_ok=True)
    csv_file = os.path.join(output_dir, f"{pdb_id}_binding_sites.csv")
    df.to_csv(csv_file, index=False)

    print(f"✅ ATP 结合点预测完成: {csv_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="预测 ATP 结合位点")
    parser.add_argument("pdb_file", type=str, help="PDB 文件路径")
    parser.add_argument("chain_id", type=int, help="目标链索引")
    args = parser.parse_args()

    predict_atp_binding(args.pdb_file, args.chain_id)