# NO_highscore_result

This is a PyMOL-based batch visualization pipeline for protein structure alignment.  
It reads query–template pairs from `zscore.csv`, aligns structures at the chain level, highlights the overlapping region, focuses on ATP/nucleotide binding sites**, and exports high-resolution images**, aligned PDB files, and a summary table.

---

## Folder Structure

NO_highscore_result/
├── pdb_results/                 # Predicted/query PDB files
├── template/                    # Template/reference PDB files
├── pymol_chain_aligned_output/  # Generated output (images/PSE/aligned PDB/summary)
├── zscore.csv                   # Input CSV (see format below)
└── pymol_result.py               # Main PyMOL script
```

---

## Dependencies
- **PyMOL** with Python integration (`pymol` command must work in terminal)
- Python standard library modules: `csv`, `glob`, `os` (already imported in the script)

> The script runs **inside PyMOL** — no extra third-party Python packages are needed.

---

## Input Format (`zscore.csv`)

The script **requires** the following columns (case-sensitive):

- `Source_File` – e.g. `Q67ZM7_zscores.csv` (the script extracts `Q67ZM7` as `query_id`)
- `filename` – alignment record file; the **last 5 characters** must contain PDB ID + chain ID  
  - Example: `...7tgkD.txt` → `pdb_id = 7tgk`, `chain_id = D`
- `Z-score` – floating-point number (stored in the summary, not used in alignment)

Other columns are ignored.

--

## Pre-run Setup

In `main()` inside `pymol_result.py`, adjust paths to match your environment:

```python
base_dir = "/path/to/project/root"
pipeline_dir = os.path.join(base_dir, "NO_highscore_result")

csv_file_path = os.path.join(pipeline_dir, "zscore.csv")
predicted_dir = os.path.join(pipeline_dir, "pdb_results")
template_dir  = os.path.join(pipeline_dir, "template")
output_dir    = os.path.join(pipeline_dir, "pymol_chain_aligned_output")
```

---

From the `NO_highscore_result/` directory:

```bash
pymol -cq pymol_result.py
```

- **`-c`** – no GUI (command-line mode)  
- **`-q`** – quiet mode  
- **`-q`** + ray tracing in the script ensures high-quality rendering

The script will prompt whether to limit the number of processed entries:  
- Press **Enter** → process all  
- Type a number → process the first N entries

---

## Output Files

Generated in `pymol_chain_aligned_output/` for each query–template pair:

1. View 1 – Chain alignment
   - `*_view1_chains.png`  
   - Coloring: Query chain = `cyan`, Template chain = `orange`; other chains = `gray` / `orange` accordingly
2. View 2 – Overlap + nucleotide
   - `*_view2_overlap_with_atp.png`  
   - Shows aligned region and ATP/ADP/AMP/GTP/GDP/GMP/CTP/CDP/CMP/UTP/UDP/UMP if present
3. View 3 – ATP/nucleotide close-up
   - `*_view3_atp_closeup.png` (4000×3000, ray-traced, white background)  
   - `*_view3_interactive.pse` (PyMOL session file)  
   - Coloring: Query side chain = `forest`, Template side chain = `gray50`, nucleotide = `red`
4. Aligned PDB
   - `*_aligned.pdb`  

---