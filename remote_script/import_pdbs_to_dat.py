import os
import subprocess

PDB_DIR = "predicted_structures"
DAT_DIR = "query_structures_DAT"
IMPORT_PL = "/home/wenhao/6tx0/software/dali/DaliLite.v5/bin/import.pl"

os.makedirs(DAT_DIR, exist_ok=True)

def import_one(pdb_path, pdb_id, dat_dir):
    print(f"📥 Importing {pdb_path} as ID: {pdb_id}")
    cmd = [
        IMPORT_PL,
        "--pdbfile", pdb_path,
        "--pdbid", pdb_id,
        "--dat", dat_dir,
        "--clean"
    ]
    subprocess.run(cmd, check=True)

# Iterate through all .pdb files
for fname in os.listdir(PDB_DIR):
    if fname.endswith(".pdb"):
        basename = os.path.splitext(fname)[0]
         # Take first 4 characters as pdbid, convert to lowercase, pad to valid ID
        pdb_id = (basename[:4].lower() + "xxx")[:4]  # Ensure length 4
        full_path = os.path.join(PDB_DIR, fname)
        import_one(full_path, pdb_id, DAT_DIR)
