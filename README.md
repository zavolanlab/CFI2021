Traverse to the desired path in your file system, then clone the repository and
move into it with:
```bash
git clone https://github.com/zavolanlab/CFI2021.git
cd CFI2021
Scripts for running KSEA analysis are located in the folder `KSEA_Phospho_Proteomics_Analysis`.
```
To reproduce Fig. S3, run the following `R` script:
```bash
Rscript Figure_S3/CDFs_proteomics_transcriptomics_targets.R
```
Files with panels of Fig. S3 `Figure_S3/Cumulative_distribution_proteomics_targets_nontargets.pdf` and 'Figure_S3/Cumulative_distribution_transcriptomics_targets_nontargets.pdf' will be generated.
