# remote_scripts

A collection of remote bioinformatics pipeline scripts for protein structure prediction, structural alignment, and database integration, designed to run on the remote server:

```
ssh webserver@coulomb.phys.ucl.ac.uk
```

---

## Folder Contents

```
remote_scripts/
├── alphafold_prediction.py     # Submit and manage AlphaFold structure predictions
├── batch_remote_pipeline.py    # Batch execution of the full remote analysis pipeline
├── convert_cif_to_pdb.py       # Convert CIF structure files to PDB format
├── DALI_align.py               # Run DALI structural alignments and parse results
├── import_pdbs_to_dat.py       # Import PDB files into Protenix DAT format
├── remote_pipeline.py          # Single-sequence remote processing workflow
```

---

## Remote Environment Information

- **Server:** `coulomb.phys.ucl.ac.uk`
- **Username:** `webserver`
- **Password:** `123456` *(do not hardcode in scripts — enter manually when prompted)*
- **Default:** Login starts in a conda environment
- **Installed tools:**
  - **Protenix** (requires activating `protoneix` environment)
  - **DALI**
  - **SUPFAM**
- **Protenix environment activation:**
  ```bash
  conda activate protoneix
  ```

---

## Dependencies (remote server)

Already installed on the remote server:
- Protenix
- DALI
- SUPFAM
- Python 3.x with standard libraries
- Conda for environment management

---

## How to Use

### 1. Connect to the server
```bash
ssh webserver@coulomb.phys.ucl.ac.uk
```
Enter password: `123456`.

### 2. Activate Protenix environment (if needed)
```bash
conda activate protoneix
```

### 3. Run individual scripts
```bash
python alphafold_prediction.py
python DALI_align.py
python convert_cif_to_pdb.py
python import_pdbs_to_dat.py
```

### 4. Run the full pipeline for one sequence
```bash
python remote_pipeline.py
```

### 5. Run batch pipeline for multiple sequences
```bash
python batch_remote_pipeline.py
```

---

## Script Descriptions

### `alphafold_prediction.py`
Runs AlphaFold structure prediction jobs on the remote server.

### `batch_remote_pipeline.py`
Automates the processing of multiple sequences in a batch, including:
- AlphaFold prediction
- DALI structural alignment
- Protenix database import

### `convert_cif_to_pdb.py`
Converts **CIF** files (common AlphaFold output) to **PDB** format for compatibility with downstream tools.

### `DALI_align.py`
Runs DALI to compare protein structures and calculates Z-scores, RMSD, and alignments.

### `import_pdbs_to_dat.py`
Imports predicted or existing PDB files into **Protenix DAT** format for use in Protenix.

### `remote_pipeline.py`
A single-sequence pipeline that performs:
1. AlphaFold prediction
2. CIF→PDB conversion
3. DALI alignment
4. Protenix DAT import

---