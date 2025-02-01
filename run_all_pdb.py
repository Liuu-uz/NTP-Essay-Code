import os
import subprocess

# ✅ 读取 PDB ID 列表（从 txt 文件）
def load_pdb_list(file_path="pdb_list.txt"):
    """从 txt 文件读取 PDB ID 列表"""
    if not os.path.exists(file_path):
        print(f"❌ 文件 {file_path} 不存在，请提供 PDB ID 列表！")
        return []
    
    with open(file_path, "r") as f:
        pdb_ids = [line.strip().upper() for line in f.readlines() if line.strip()]
    
    return pdb_ids

# ✅ 1️⃣ 下载 PDB 结构
def fetch_pdb(pdb_id):
    """调用 fetch_pdb.py 下载 PDB 文件"""
    pdb_code = pdb_id[:-1]  # 提取 PDB ID（去掉最后一个字符）
    pdb_file = f"pdb_files/{pdb_code}.pdb"
    if not os.path.exists(pdb_file):
        print(f"🔍 下载 PDB: {pdb_code} ...")
        subprocess.run(["python", "fetch_pdb.py", pdb_code])  # 运行 fetch_pdb.py
    else:
        print(f"✅ PDB {pdb_code} 已存在")

# ✅ 2️⃣ 验证 PDB 结构
def validate_pdb(pdb_id):
    """调用 validate_pdb.py 检查 PDB 结构"""
    pdb_code = pdb_id[:-1]
    pdb_file = f"pdb_files/{pdb_code}.pdb"
    
    if not os.path.exists(pdb_file):
        print(f"❌ PDB 文件 {pdb_file} 不存在，跳过验证")
        return False
    
    print(f"🔍 验证 PDB 结构: {pdb_code} ...")
    result = subprocess.run(["python", "validate_pdb.py", pdb_file], capture_output=True, text=True)
    if "❌" in result.stdout:
        print(f"❌ {pdb_code} 验证失败，跳过处理")
        return False
    
    print(f"✅ PDB {pdb_code} 通过验证")
    return True

# ✅ 3️⃣ 解析 PDB 结构并处理
def process_pdb(pdb_id):
    """调用 proc.py 解析 PDB 结构"""
    pdb_code = pdb_id[:-1]
    pdb_file = f"pdb_files/{pdb_code}.pdb"
    
    if not os.path.exists(pdb_file):
        print(f"❌ PDB 文件 {pdb_file} 不存在，跳过处理")
        return
    
    print(f"🔍 解析 {pdb_code} ...")
    subprocess.run(["python", "proc.py", pdb_file])  # 运行 proc.py

# ✅ 4️⃣ 预测 ATP 结合点
def predict_atp_binding(pdb_id):
    """调用 predict.py 预测 ATP 结合点"""
    pdb_code = pdb_id[:-1]  # 提取 PDB ID
    chain_id = pdb_id[-1]  # 提取链 ID
    pdb_file = f"pdb_files/{pdb_code}.pdb"
    
    if not os.path.exists(pdb_file):
        print(f"❌ PDB 文件 {pdb_file} 不存在，跳过预测")
        return
    
    print(f"🔍 预测 ATP 结合点: {pdb_code} (Chain {chain_id}) ...")
    subprocess.run(["python", "predict.py", pdb_file, str(ord(chain_id) - ord('A'))])  # 运行 predict.py

# ✅ 运行所有 PDB ID
if __name__ == "__main__":
    os.makedirs("pdb_files", exist_ok=True)  # 确保 PDB 目录存在
    
    pdb_ids = load_pdb_list("pdb_list.txt")  # 读取 PDB ID

    if not pdb_ids:
        print("❌ 没有找到 PDB ID，请检查 `pdb_list.txt` 文件！")
    else:
        for pdb in pdb_ids:
            fetch_pdb(pdb)  # 下载 PDB
            if validate_pdb(pdb):  # 先验证 PDB
                process_pdb(pdb)  # 处理 PDB
                predict_atp_binding(pdb)  # 预测 ATP 结合点

        print("✅ 批量处理完成！")