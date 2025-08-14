import os
import requests

def fetch_pdb(pdb_id, save_dir="pdb_files"):
    """
    Download PDB structure from RCSB PDB and save to local directory
    """
    file_path = os.path.join(save_dir, f"{pdb_id}.pdb")
    if os.path.exists(file_path):
        print(f"✅ {pdb_id}.pdb already exists, skipping download")
        return file_path
    
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    os.makedirs(save_dir, exist_ok=True)
    
    response = requests.get(url)
    if response.status_code == 200:
        with open(file_path, "w") as f:
            f.write(response.text)
        print(f"✅ Download completed: {pdb_id}.pdb -> {file_path}")
    else:
        print(f"❌ Download failed: {pdb_id}.pdb (PDB ID may not exist)")
    
    return file_path

# Usage example
fetch_pdb("3I2B")