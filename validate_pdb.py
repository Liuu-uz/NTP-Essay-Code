import argparse
import mdtraj as md

def validate_pdb(pdb_file):
    """
    检查 PDB 文件的完整性：
    - 确保存在目标链
    - 检查是否包含蛋白质
    """
    try:
        traj = md.load(pdb_file)
        top = traj.topology
        num_residues = sum(1 for _ in top.residues)  # ✅ 解决 generator 问题
        num_atoms = sum(1 for _ in top.atoms)  # ✅ 解决 generator 问题
        protein_atoms = list(top.select("protein"))  # ✅ 转换为列表，避免 generator 错误


        print(f"✅ {pdb_file} 加载成功")
        print(f"🔹 总残基数: {num_residues}, 总原子数: {num_atoms}")

        if len(protein_atoms) == 0:
            print(f"⚠️ {pdb_file} 可能没有蛋白质结构！")

    except Exception as e:
        print(f"❌ 解析 {pdb_file} 失败: {e}")

validate_pdb("pdb_files/7V0F.pdb")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="检查 PDB 结构完整性")
    parser.add_argument("pdb_file", type=str, help="PDB 文件路径")
    args = parser.parse_args()

    validate_pdb(args.pdb_file)