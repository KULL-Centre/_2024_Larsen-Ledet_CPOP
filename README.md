# Systematic characterization of indel variants using a yeast-based protein folding sensor
## Introduction
This respository contains all data (except from the raw FASTQ files, which are available at the NCBI Gene Expression Omnibus (GEO) repository (accession number: GSE270811)) and code to repeat the processing and analysis of the CPOP data in Larsen-Ledet et al.: "Systematic characterization of indel variants using a yeast-based protein folding sensor".

## Overview of files
*Output files*
* **cpop_data.csv** - CPOP scores and standard deviations for DHFR indel, synonymous and nonsense variants.
* **cpop_data_pre_rescale** - Raw CPOP scores and standard deviations for DHFR indel, synonymous and nonsense variants prior to rescaling.
* **cpop_data_ROC_[ins|del].csv** - CPOP scores for ROC curves, where duplicated indel variants on protein level have been removed.
* **tile[1-5].csv** - Counts per tile for DHFR indel, synonymous and nonsense variants for each replicate and condition.
  
*Input files*
* **[ins|del]_dplddt_ddg.csv** - dpLDDT and ddG predictions for DHFR indel variants.
* **rSASA.csv** - Relative solvent accessible surface area (rSASA) for each residue in DHFR.
* **mtx_dist.csv** - Distance (Å) of each residue in DHFR to the MTX binding site.
* **[ins|del]_esm1b** - ESM1b predictions for DHFR indel variants.
* **del_sequence_alignment** - MSA generated with HHblits of DHFR homologs with deletions per position.

*Excel files*
* **CPOP_primers_annealing.temp..xlsx** - Primers and annealing temperatures for the first PCR in amplicon preparation.
* **SupplementalFile1.xlsx** - All data files combined in a single Excel file.

## Processing of raw sequencing data
The function.py file is used to call DHFR variants and calculate CPOP scores. The script takes raw FASTQ files as input. The output is a dataset with CPOP scores and standard deviations for DHFR indel, synonymous and nonsense variants.

## Data analysis and plotting
The CPOP_data_analysis.R file is used to produce all plots in the main figures, and the CPOP_data_analysis_supplementary.R file is used to produce all plots in the supplementary figures. Both files take the dataset with CPOP scores and standard deviations as input.

## Preprint
https://doi.org/10.1101/2024.07.11.603017
