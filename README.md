# CRC miRNA Random Forest Signature

Computational workflows supporting an accepted article on the development and
preliminary external validation of a three-miRNA random forest signature for
colorectal cancer tissue discrimination.

The article has been accepted by *ARO—The Scientific Journal of Koya
University* and is currently in production. The final citation and DOI will be
added when available.

## Locked miRNA signature

- hsa-miR-135b-5p
- hsa-miR-378a-5p
- hsa-miR-549a-5p

## Public datasets

| Role | Accession | Design |
|---|---|---|
| Discovery | PRJNA387172 | 104 matched tumour–normal pairs |
| External validation | GSE125961 | 6 matched pairs |
| External validation | GSE183437 | 5 matched pairs |
| External validation | GSE236991 | 5 matched pairs |

## Repository structure

```text
CRC-miRNA-RF-signature/
├── README.md
├── CITATION.cff
├── .gitignore
├── data/
│   └── README.md
├── original-scripts/
│   ├── original submitted script files
│   └── SHA256SUMS.txt
└── scripts/
    ├── preprocessing/
    │   ├── 01_discovery_PRJNA387172_pipeline.sh
    │   ├── 02_GSE125961_pipeline.sh
    │   ├── 03_GSE183437_pipeline.sh
    │   └── 04_GSE236991_pipeline.sh
    └── analysis/
        ├── 05_machine_learning_pipeline.R
        ├── 06_uncertainty_calibration_permutation.R
        └── 07_target_GO_KEGG_analysis.R
```

`original-scripts/` preserves the supplied files without analytical
modification. `scripts/` provides the same workflows under descriptive,
execution-ordered filenames and appropriate file extensions.

## Recommended execution order

1. Run the discovery Ubuntu pipeline.
2. Run each external-cohort Ubuntu pipeline.
3. Run the machine-learning pipeline.
4. Run the uncertainty, calibration and permutation analyses.
5. Run the target-integration and GO/KEGG enrichment analysis.

## Main computational environment

Primary read processing was performed under Ubuntu 24.04.2 LTS. Downstream
analysis was performed in R 4.5.2 on Windows 11 x64. Principal software and R
packages are documented within the scripts and the accepted manuscript.

The workflows include machine-learning feature selection with Elastic Net,
Boruta and RF-RFE; a random forest classifier with 1,000 trees; repeated
five-fold cross-validation with 10 repeats; confidence intervals, calibration,
Brier scores and permutation analyses; and multiMiR-based target integration
with GO and KEGG enrichment.

## Reproducibility notes

- The scripts reflect the computational workflow used for the accepted article.
- Several scripts contain workstation-specific paths. Update only the user
  settings and path variables before execution.
- Raw sequencing data are not redistributed because they are available from
  public SRA/GEO accessions.
- Generated intermediate and result files are excluded from version control by
  `.gitignore`.

## Citation

Until the final bibliographic record is available, cite this repository together
with the accepted article. Update `CITATION.cff` after the DOI is assigned.

