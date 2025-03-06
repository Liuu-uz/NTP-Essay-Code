import os
import requests
from Bio.PDB import MMCIFParser
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.Data.IUPACData import protein_letters_3to1

# 定义文件夹
cif_dir = os.path.abspath("cif_files")
results_dir = os.path.abspath("results")
os.makedirs(cif_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

# 读取 PDB ID 列表
pdb_list_file = "pdb_list.txt"
if not os.path.exists(pdb_list_file):
    print(f"❌ 文件 {pdb_list_file} 不存在！")
    exit(1)

with open(pdb_list_file, "r") as f:
    pdb_data = [line.strip() for line in f.readlines() if line.strip()]

# 解析 PDB ID 和链
def parse_pdb_entry(entry):
    pdb_id = entry[:-1].lower()  # 提取 PDB ID（去掉最后一个字符）
    chain = entry[-1].upper()  # 提取链 ID（最后一个字符）
    return pdb_id, chain

pdb_entries = [parse_pdb_entry(entry) for entry in pdb_data]

# 下载 CIF 文件
def fetch_cif(cif_id, save_dir=cif_dir):
    file_path = os.path.join(save_dir, f"{cif_id}.cif")
    if os.path.exists(file_path):
        print(f"✅ {cif_id}.cif 已存在，跳过下载")
        return file_path

    url = f"https://files.rcsb.org/download/{cif_id}.cif"
    response = requests.get(url)

    if response.status_code == 200:
        with open(file_path, "w") as f:
            f.write(response.text)
        print(f"✅ 下载完成: {cif_id}.cif -> {file_path}")
        return file_path
    else:
        print(f"❌ 下载失败: {cif_id}.cif (可能 CIF ID 不存在)")
        return None

# 解析 CIF 并转换为 FASTA
def cif_to_fasta(cif_file, pdb_id, selected_chain):
    parser = MMCIFParser(QUIET=True)

    try:
        structure = parser.get_structure(pdb_id, cif_file)
    except Exception as e:
        print(f"❌ 解析 CIF 文件 {cif_file} 失败: {e}")
        return None

    sequence = []
    found_chain = False

    for model in structure:
        for chain in model:
            if chain.id == selected_chain:
                found_chain = True
                print(f"🔍 处理链 {chain.id} (PDB: {pdb_id})")
                for residue in chain:
                    if residue.id[0] == " ":  # 只处理标准氨基酸
                        residue_name = residue.resname.strip().capitalize()
                        if residue_name in protein_letters_3to1:
                            sequence.append(protein_letters_3to1[residue_name])
                        else:
                            print(f"⚠️ 跳过未知氨基酸: {residue_name}")
                break  # 找到目标链后停止

    if not found_chain:
        print(f"❌ 未找到链 {selected_chain} (PDB: {pdb_id})")
        return None

    if sequence:
        fasta_sequence = "".join(sequence)
        fasta_file = os.path.join(results_dir, f"{pdb_id}_{selected_chain}.fasta")
        record = SeqRecord(Seq(fasta_sequence), id=f"{pdb_id}_{selected_chain}", description="")
        with open(fasta_file, "w") as f:
            FastaIO.FastaWriter(f).write_record(record)
        print(f"✅ FASTA 文件已保存: {fasta_file}")
        return fasta_file
    else:
        print(f"❌ {pdb_id}_{selected_chain} 序列为空，未生成 FASTA")
        return None

# 处理所有 PDB ID
for pdb_id, chain in pdb_entries:
    cif_file = fetch_cif(pdb_id)  # 下载 CIF 文件
    if cif_file:
        cif_to_fasta(cif_file, pdb_id, chain)  # 解析并转换为 FASTA

# 输出转换后的 FASTA 结果
fasta_files = [os.path.join(results_dir, f) for f in os.listdir(results_dir) if f.endswith(".fasta")]
for fasta_file in fasta_files:
    print(f"\n📄 FASTA 文件内容 ({fasta_file}):")
    with open(fasta_file, "r") as f:
        print(f.read())