A Stepwise Protein Identification Framework with Structural Visualization

Here is a brief introduction of my codes and results content
repo/
├─ README.md                      # this file              
├─ batch_results/                 # results of DALI, collecting all zscore csv files for each sequence
├─ prediction_sequences/          # json files for each sequence. Input files to prediction structure
├─ DAT_template/                  # remote folder contains template files used for DALI
├─ fasta_files/                   # initial data，need to pre-process
├─ fasta_sequence/                # generate from fasta_files
├─ html_results/                  # results of supfam
├─ fasta_sequence/                # generate from fasta_files
├─ ntp_processing_sequences_with_ligand/                #Input files for protenix(YES branch, only for visualization)
├─ pdb_files/                     # all template pdb files
├─ supfam_pdb_files/              # results of protenix(YES branch) without ATP
├─ yes_pdb_results/              # results of protenix(YES branch) with ATP
├─ three_view_analysis_output/    # visualization for YES branch
├─ NO_highscore_results/          # visualization for sequences that find superfamily by prediction and DALI process
│  ├─ README.md
│  └─ <source codes>
├─ remote_script/                 # code files need to be send to remote server
│  ├─ README.md
│  └─ <source codes>
├─seq-supfam 
│  ├─ README.md     
│  └─ <source codes> 
├─ <download_fasta.py>/           # 
