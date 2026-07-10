## Overview of nf-pipelines/
Holds nextflow pipelines originally from [philippinespire/nf-pipelines](https://github.com/philippinespire/nf-pipelines/).

## Content

```
nf-pipelines/
├── README.md
├── environments/
│   ├── nf-angsd-diversity.yml      # conda environment for nf-angsd
├── nf-angsd-diversity/             # angsd, PCA, gen div for all individuals
├── nf-angsd-diversity-cvi-only/    # same, but trimming out C. atripectoralis.
├── nf-trim-generode/               # replacement for nf-trim-merged-unmerged (done later)
└── nf-trim-merged-unmerged/        # read trimming and mapping and amber (done originally)
```
