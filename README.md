# CFI 2021

Traverse to the desired path in your file system, then clone the repository and
move into it with:

```bash
git clone https://github.com/zavolanlab/CFI2021.git
cd CFI2021
```

Different analyses are described in different directories:

* Exact commands related to RNA-Seq data processing are described in
  `Preprocessing/README.md`.

* Reproducing the differential gene expression analyses is described in
  `Differential_Gene_Expression_Analysis/README.md`.

* Information on creating coverage profile plots were collected in
  `Coverage_Plots/plot-coverages.R`.

* All the analyses related to transcripts' 3' ends are gathered in the
  following notebook: `Terminal_Exons_Analysis/report.ipynb`.

  > Please follow the code there to reproduce plots presented in Fig. 1.

* Scripts for running KSEA analysis are located in the folder
  `KSEA_Phospho_Proteomics_Analysis`.

* To reproduce Fig. S3, run the following `R` script:

  ```bash
  Rscript Figure_S3/CDFs_proteomics_transcriptomics_targets.R
  ```

  Files with panels of Fig. S3
  `Figure_S3/Cumulative_distribution_proteomics_targets_nontargets.pdf` and
  `Figure_S3/Cumulative_distribution_transcriptomics_targets_nontargets.pdf`
  will be generated.
