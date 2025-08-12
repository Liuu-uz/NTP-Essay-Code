from Bio.PDB import MMCIFParser, PDBIO
import sys
import os

if len(sys.argv) != 3:
    print("Usage: python convert_cif_to_pdb.py input.cif output.pdb")
    sys.exit(1)

cif_file = sys.argv[1]
pdb_file = sys.argv[2]

if not os.path.isfile(cif_file):
    print(f"Error: File not found: {cif_file}")
    sys.exit(1)

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("model", cif_file)

# Added: Force set single chain to A and ensure consecutive residue numbering
for model in structure:
    chains = list(model.get_chains())
    if chains:
        # Rename first chain to A
        first_chain = chains[0]
        first_chain.id = 'A'
        
        # Ensure consecutive residue numbering (starting from 1)
        residues = list(first_chain.get_residues())
        for i, residue in enumerate(residues, 1):
            residue.id = (' ', i, ' ')  # (hetflag, resseq, icode)
        
        # If multiple chains exist, keep only the first one
        for chain in chains[1:]:
            model.detach_child(chain.get_id())

io = PDBIO()
io.set_structure(structure)
io.save(pdb_file)

print(f"Converted: {cif_file} -> {pdb_file}")

# Added: Validate output
with open(pdb_file, 'r') as f:
    atom_lines = [line for line in f if line.startswith('ATOM')]
    if atom_lines:
        print(f"Generated {len(atom_lines)} ATOM records")
        print(f"Chain ID: {atom_lines[0][21]}")
        print(f"First residue: {atom_lines[0][22:26].strip()}")
        print("First line:", atom_lines[0].strip())