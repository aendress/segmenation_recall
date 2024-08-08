# Recall experiment for word segmentation

## experiments 
This folder contains experimental file, but is not uploaded.

The Psyscope files are not uploaded since they cannot be run any longer anyhow on current versions of Mac OS X. As a result, we cannot make sure that they actually work. 

The online version of the recall experiment is available at https://www.testable.org/library#496496

## results

This folder contains data and analysis scripts.


### Analysis scripts

####  `segmentation_recall_combined.Rmd`
`segmentation_recall_combined.Rmd` is the main analysis script. An updated version working with R 4.1.1 is `segmentation_recall_combined_for_R_version_4.1.1.Rmd`

#### `segmentation_recall_combined_for_revision.Rmd`

Version working with R 4.1.1 and containing updated figures/tables. This file contains the following changes:
		+ Tables report GLMM results include odds ratios
		+ Some multiple panel figures have been separated into separated figures
		+ There is a single figure for all conditions in Experiment 2 (i.e., the non-recall conditions)
		+ Analyses after removing outliers.
		
#### `segmentation_recall_combined_for_revision2.Rmd`

Identical to `segmentation_recall_combined_for_revision.Rmd`, with the following changes. 

- Added the analysis for the simulations with PARSER. (The data for the simulations with Endress & Johnson's model is not included here.). The analyses are in subsection "Can a chunking model account for these results? Simulations with PARSER (added for revision for Cognitive Psychology)"
- Added analyses of correlations between recall and recognition performance. These analyses are given in subsection "Correlations (added for revision for Cognitive Psychology)"

### data

- `recall_city`: Lab-based version of Experiment 1
- `recall_testable`: Online version of Experiment 1. Participants have been anonymized
- `oversegmentation_city`: Experiment 2
- `oversegmentation_bcn`: pilot Experiments

### helper
Various helper functions

		
## simulations

Files for the various simulations.

	- parser: Simulations with PARSER. See Readme inside the folder for details
	- TP_model_subunits: Simulations with Endress & Johnson's model. Main file is tp_model_subunits.Rmd		
	



