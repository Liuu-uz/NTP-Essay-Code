# A Stepwise Protein Identification Framework with Structural Visualization

This repository contains a stepwise pipeline for protein identification, combining sequence analysis, superfamily classification, structure prediction, and structural visualization using tools such as SUPFAM, DALI, and Protenix.

---

## Repository Structure

```
repo/
├─ README.md                      # This file              
│
├─ batch_results/                 # DALI results: all zscore CSV files per sequence
├─ prediction_sequences/          # JSON input files for structure prediction (per sequence)
├─ DAT_template/                   # Remote folder containing DALI template files
├─ fasta_files/                    # Initial raw FASTA data (needs preprocessing)
├─ fasta_sequence/                 # Clean FASTA sequences generated from fasta_files
├─ html_results/                    # SUPFAM analysis results (HTML format)
├─ ntp_processing_sequences_with_ligand/ # Input files for Protenix (YES branch, visualization only)
├─ pdb_files/                       # Template PDB files for alignment
├─ supfam_pdb_files/                # Protenix results (YES branch) without ATP
├─ yes_pdb_results/                  # Protenix results (YES branch) with ATP bound
├─ three_view_analysis_output/      # Structural visualization outputs for YES branch
├─ NO_highscore_results/            # Visualization for sequences that found superfamily via prediction + DALI
│  ├─ README.md                     # Instructions for NO_highscore_results workflow
│  └─ <source codes>
├─ remote_script/                   # Scripts to be sent to remote server
│  ├─ README.md
│  └─ <source codes>
├─ seq-supfam/                      # SUPFAM sequence processing and visualization
│  ├─ README.md
│  └─ <source codes>
│
├─ download_fasta.py                # Script to download FASTA files from remote server
├─ fetch_pdb.py                     # Script to download PDB structures from online databases
├─ process_fasta.py                  # Script to preprocess raw FASTA files into fasta_sequence
```

---

## 🔄 Workflow Overview

The framework follows these main steps:

1. **FASTA file acquisition**
   - Use `download_fasta.py` to fetch sequence files from the remote server.
   - Alternatively, retrieve raw data and place them in `fasta_files/`.

2. **FASTA preprocessing**
   - Run `process_fasta.py` to clean and standardize sequences into `fasta_sequence/`.

3. **Superfamily classification (SUPFAM)**
   - Process sequences with the **seq-supfam** module.
   - Store SUPFAM results in `html_results/`.

4. **Structure prediction**
   - Use `prediction_sequences/` as input for AlphaFold (via `remote_script/`).
   - Predicted structures are aligned against template PDB files in `DAT_template/`.

5. **Structural alignment (DALI)**
   - DALI compares predicted/query structures with templates.
   - Output z-score CSVs in `batch_results/`.

6. **Protenix processing**
   - For sequences with ATP-binding (YES branch):
     - Input: `ntp_processing_sequences_with_ligand/`
     - Output without ATP: `supfam_pdb_files/`
     - Output with ATP: `yes_pdb_results/`
     - Visualization: `three_view_analysis_output/`

7. **NO highscore visualization**
   - For sequences that only find superfamily through prediction + DALI:
     - Output in `NO_highscore_results/` (see its README for details)

---

## External Dependencies

- **SUPFAM** (sequence classification)
- **DALI** (structural alignment)
- **Protenix** (ligand visualization, DAT handling)
- **AlphaFold** (structure prediction, run remotely)
- Python packages:
  - `pandas`
  - `requests`
  - `biopython`
  - `matplotlib`
  - `openpyxl`

---

## Example Usage

### 1. Download sequences from remote server
```bash
python download_fasta.py
```

### 2. Process raw FASTA files
```bash
python process_fasta.py
```

### 3. Fetch PDB templates
```bash
python fetch_pdb.py
```

### 4. Run remote pipeline for prediction + alignment
```bash
# In remote_script/ folder, see README for environment activation
ssh webserver@coulomb.phys.ucl.ac.uk
python remote_pipeline.py
```

### 5. Visualize results
- YES branch: see `three_view_analysis_output/`
- NO highscore branch: see `NO_highscore_results/`

---