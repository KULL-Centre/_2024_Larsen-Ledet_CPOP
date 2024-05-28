# Systematic characterization of indel variants using a yeast-based protein folding sensor
## Introduction
This respository contains all data and code to repeat the processing and analysis of the CPOP data in Larsen-Ledet et al.: "Systematic characterization of indel variants using a yeast-based protein folding sensor"

## Overview of files
* cpop_data.csv - CPOP scores and standard deviations for DHFR indel, synonymous and nonsense variants.
* del/ins_dplddt_ddg.csv - dpLDDT and ddG predictions for DHFR indel variants.
* rSASA.csv - Relative solvent accessible surface area (rSASA) for each residue in DHFR.
* mtx_dist.csv - Distance (Å) of each residue in DHFR to the MTX binding site.
* cpop_data_ROC_del/ins.csv - CPOP scores for ROC curves, where duplicated indel variants on protein level have been removed.
* tile1-5.csv - Counts for DHFR indel, synonymous and nonsense variants for each replicate and condition.
* CPOP_primers_annealing.temp..xlsx - Primers and annealing temperatures for the first PCR in amplicon preparation.
* CPOP_data_combined.xlsx - All data files combined in a single Excel sheet

## Processing of raw sequencing data
The function.py file is used to call DHFR variants and calculate CPOP scores. The script takes raw FASTQ files as input, which are available at the NCBI Gene Expression Omnibus (GEO) repository (accession number: XXXX).

## Data analysis and plotting
The CPOP_data_analysis.R file is used to analyze and produce all plots of the main figures. The CPOP_data_analysis_supplementary.R file is used to analyze and produce all plots of the supplementary figures.
