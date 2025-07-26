from Bio.PDB import MMCIFParser, PDBIO
import sys
import os

if len(sys.argv) != 3:
    print("Usage: python convert_cif_to_pdb.py input.cif output.pdb")
    sys.exit(1)

cif_file = sys.argv[1]
pdb_file = sys.argv[2]

if not os.path.isfile(cif_file):
    print(f"❌ File not found: {cif_file}")
    sys.exit(1)

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("model", cif_file)

# 新增：强制设置单链为A，并确保残基编号连续
for model in structure:
    chains = list(model.get_chains())
    if chains:
        # 重命名第一条链为A
        first_chain = chains[0]
        first_chain.id = 'A'
        
        # 确保残基编号连续（从1开始）
        residues = list(first_chain.get_residues())
        for i, residue in enumerate(residues, 1):
            residue.id = (' ', i, ' ')  # (hetflag, resseq, icode)
        
        # 如果有多条链，只保留第一条
        for chain in chains[1:]:
            model.detach_child(chain.get_id())

io = PDBIO()
io.set_structure(structure)
io.save(pdb_file)

print(f"✅ Converted: {cif_file} → {pdb_file}")

# 新增：验证输出
with open(pdb_file, 'r') as f:
    atom_lines = [line for line in f if line.startswith('ATOM')]
    if atom_lines:
        print(f"📊 Generated {len(atom_lines)} ATOM records")
        print(f"📊 Chain ID: {atom_lines[0][21]}")
        print(f"📊 First residue: {atom_lines[0][22:26].strip()}")
        print("First line:", atom_lines[0].strip())

