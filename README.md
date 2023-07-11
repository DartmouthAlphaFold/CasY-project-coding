# Private storage of all codes used for CasY structure prediction project

### Folder ./alphafold
Scripts for running AlphaFold on Dartmouth's Discovery cluster:
  - `run_alphafold_singularity_gpu.py`: run AlphaFold Singularity instance
  - `visualize_alphafold_results.py`: visualize AlphaFold results
  - `run_alphafold_singularity_gpu.sh`: submit the job via Slurm
  - `alphafold_readme.md`: instruction 
  - and other scripts in testing

### Folder ./chimerax
Scripts for comparing protein structures using ChimeraX
  - `matchmaker.py`: pairwise matchmake a given set of proteins and generate distance matrices
  - `heatmap.py`: build a phylogeny using the distance matrices
  - `chimerax_readme.md`: instruction

### Folder ./phmmer
Scripts for performing pHMMER protein search and downloading sequences, using APIs hosted by European Bioinformatics Institute
  - `phmmer.py`: search against RefProt, UniProtKB, SwissProt, PDB, AlphaFoldDB, etc [databases](https://www.ebi.ac.uk/Tools/hmmer/search/phmmer)
  - `phmmer_mgnify.py`: search against [MGnify database](https://www.ebi.ac.uk/metagenomics/sequence-search/search/phmmer) 
  - `sync_*.py`: scripts in testing for speed optimization
  - `phmmer_readme.md`: instruction
  - and other scripts in testing
