# vcs_2019
# Nick Hasle, 11.22.19
# Visual Cell Sorting: A High-throughput, Microscope-based Method to Dissect Cellular Heterogeneity
## GitHub Repository readme

## A note on data storage
All files > 20M in size were bzipped OR are stored in the GEO repository (accession number XXXXX). Data directories with missing data include a README.txt file that denotes which GEO repository files ought to be downloaded.

## Directories
### bin
contains monocle3 and batchelor github downloads for use in vcs_paclitaxel

### metamorph
sitemaps: htacquir.cfg files that were copied into the MetaMorph directory to load circular/elliptical site maps for various experiments
journals: .pdf and .jnl files with the journals loaded into the Plate Acquisition dialog box for the experiments

### vcs_proof_of_concept
activation_accuracy: code and .fcs files used to assess accuracy of activation
activation_fourbin: code and .fcs files used to assess four bin activation
activation_toxicity: code and .fcs files used to assess activation-induced toxicity
activation_bulk_RNAseq: code and data used to analyze bulk RNA sequencing data

### vcs_nls
enrich2: configuration file used for enrich2
nls_imaging: Metamorph output files from analysis of NLS library pre-sort imaging data from Replicate 2, Technical Replicate 1
nls_post_sort_validation: Metamorph output files from analysis of the sorted NLS library
nls_scores: code and enrich2 output used to calculate variant scores
nls_validation: code and Metamorph output used to calculate variant NC ratios
nls_prediction: code and data used to train and validate the NLS prediction model

### vcs_paclitaxel: 
imaging: imaging data from Experiment 1 (Visual Cell Sorting experiment) taken directly before Visual Cell Sorting
scrnaseq: code and data used to analyze the single cell RNA sequencing data of Experiment 1 (cells treated with paclitaxel and sorted according to their nuclear shape) and Experiment 2 (cells treated with paclitaxel, activated as in experiment 1, and sorted into the same tube; i.e. left unseparated).