# seq-supfam

A pipeline for NTP-related SUPFAM sequence analysis, structural prediction, and visualization.  
This toolkit integrates sequence processing, SUPFAM superfamily annotation, remote AlphaFold prediction, and visualization of NTP-related structural features.

---

## Folder Structure

seq-supfam/
├── aotu_supfam.py                  # Automated SUPFAM superfamily annotation
├── remote_alphafold.py             # Submit jobs to AlphaFold remotely for structure prediction
├── sequence_to_fasta.py            # Convert sequences into FASTA format
├── supfam_checker.py                # Check SUPFAM annotation results for NTP-related families
├── supfam_visualisation.py          # Visualize SUPFAM-based NTP structure matches
├── NTP_Analysis_Report.xlsx         # Report with NTP-related sequence statistics and results
├── SUPFAM based NTP processing.xlsx # SUPFAM annotation data and processing results
```

---

## Dependencies

- Required Python packages:
  - `pandas`
  - `openpyxl`
  - `requests`
  - `biopython`
  - `matplotlib` (for visualization)
- SUPFAM annotation tool / database access
- (Optional) AlphaFold API or server access

---

## Main Scripts Overview

### `sequence_to_fasta.py`
- **Purpose:** Converts raw sequence data into FASTA format for downstream processing.
- **Input:** Sequence list (Excel/CSV)
- **Output:** FASTA file

### `aotu_supfam.py`
- **Purpose:** Automates submission of sequences to the SUPFAM database for superfamily annotation.
- **Output:** Annotation results in CSV/Excel

### `supfam_checker.py`
- **Purpose:** Filters SUPFAM annotation results to identify **NTP-related** superfamilies.
- **Output:** Filtered list of NTP-related entries

### `remote_alphafold.py`
- **Purpose:** Submits selected sequences to a remote AlphaFold service for structure prediction.
- **Output:** PDB files

### `supfam_visualisation.py`
- **Purpose:** Generates visual summaries (plots/charts) of SUPFAM annotation and NTP-related sequence distribution.

---

## 📊 Data Files

- **`NTP_Analysis_Report.xlsx`**  
  Contains processed results, statistics, and summary of NTP-related proteins.
  
- **`SUPFAM based NTP processing.xlsx`**  
  Detailed processing data linking SUPFAM results to NTP-related annotations.

---
