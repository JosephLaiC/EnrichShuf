# EnrichShuf

## How to install 

R script

```{R}
branch <- ""
token  <- "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"
repo   <- "JosephLaiC/EnrichShuf"

devtools::install_github(
  repo=repo, ref=branch,
  auth_token=token, upgrade=FALSE)
```

## Dependency (for version 0.1.0)

```
## regular package
readr
data.table
dplyr
stringr

## Bioconductor package
GenomicRanges
ChIPseeker
```

To build this package under conda environment

```
conda install -c conda-forge -c defaults -c bioconda\
  bioconductor-chipseeker\
  bioconductor-genomicranges\
  r-devtools r-readr r-stringr r-data.table
```



