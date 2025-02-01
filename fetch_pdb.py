import os
import requests

def fetch_pdb(pdb_id, save_dir="pdb_files"):
    """
    从 RCSB PDB 下载 PDB 结构，并保存到本地
    """
    file_path = os.path.join(save_dir, f"{pdb_id}.pdb")
    if os.path.exists(file_path):
        print(f"✅ {pdb_id}.pdb 已存在，跳过下载")
        return file_path

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    os.makedirs(save_dir, exist_ok=True)

    response = requests.get(url)
    if response.status_code == 200:
        with open(file_path, "w") as f:
            f.write(response.text)
        print(f"✅ 下载完成: {pdb_id}.pdb -> {file_path}")
    else:
        print(f"❌ 下载失败: {pdb_id}.pdb (可能 PDB ID 不存在)")

    return file_path

fetch_pdb("7V0F")