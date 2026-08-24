## Overview of nf-pipelines/
Holds nextflow pipelines originally from [philippinespire/nf-pipelines](https://github.com/philippinespire/nf-pipelines/).

## Contents

```
nf-pipelines/
├── README.md
├── environments/
│   ├── nf-angsd-diversity.yml      # conda environment for nf-angsd
├── nf-angsd-diversity/             # angsd, PCA, gen div for all individuals
│                                   # replaced by nf-angsd-selection
├── nf-angsd-diversity-cvi-only/    # trimmed out C. atripectoralis.
│                                   # replaced by nf-angsd-selection
├── nf-angsd-selection/             # refined pipeline that included selection and pruning
│                                   # replaced by -1x
├── nf-angsd-selection-1x/          # finalized diversity analysis
├── nf-trim-generode/               # replacement for nf-trim-merged-unmerged (done later)
│                                   # adds repeat-masking, only trims modern reads
└── nf-trim-merged-unmerged/        # read trimming and mapping and amber
│                                   # replaced by -generode
```
