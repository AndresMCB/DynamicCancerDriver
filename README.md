# DynamicCancerDriver
DynamicCancerDriver package contains functions to identify genes driving one or more bio-pathological transitions during cancer progression. Formally, we name a gene driving one or more core processes over cancer progression as dynamic cancer driver. 
### Note
Please find the results of our experiments in the folder [inst/Supplementary/](inst/Supplementary/). You can also access them from R by using the name of the file. For example:

```R
aux <- system.file("Supplementary/"
                   , "supplementary table 8 - dynamic cancer drivers HER2time(Bulk).csv"
                   , package = "DynamicCancerDriver")
DCD.HER2_bulk <- read.csv(aux)
```

## Introduction 
Our method takes gene expression data from cross-sectional studies,
as well as a **covariate** that reasonably modulates (in the sense described
by [Campbell and Yau, 2018](https://www.nature.com/articles/s41467-018-04696-6)) the pseudotemporal progression of one
relevant process occurring during cancer development. If no pseudotime
is provided, our method relies on [PhenoPath](https://www.bioconductor.org/packages/release/bioc/html/phenopath.html) package 
to find a pseudotime score to order the samples following
the trajectory encoded by the **covariate**. 

We use the pseudotime (either provided or inferred by using PhenoPath) and the covariate provided to
find a critical turning point in the trajectory along the pseudotime. We
name this critical point as the "event".

We applied our dynamic cancer drivers approach to a single cell RNA
sequencing dataset (NCBI GEO database, accession [GSE75688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688)) ([Chung
et al., 2017](https://www.nature.com/articles/ncomms15081)), and the cancer genome atlas breast cancer dataset ([TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA)).
Experiments implemented in our paper can be found as follows:
1. [demo/Test_DynamicCancerDriver(SC).R](demo/Test_DynamicCancerDriver(SC).R): Drivers inferred from a pre-processed Single Cell data, (GSE75688)
2. [demo/Test_DynamicCancerDriver(Bulk).R](demo/Test_DynamicCancerDriver(Bulk).R): Drivers inferred from the TCGA-BRCA dataset.

## Installation 
DynamicCancerDriver runs in the R statistical computing environment.

R (>=4.1.0), devtools(>=2.4.2), Rtools (>=4.0), Bioconductor (>=3.14), phenopath (tested with v. 1.18.0), CausalImpact(>= 1.2.7), and
 tidyverse(>= 1.3.1) are  required.
We also use some utilities from another of our packages ([AMCBGeneUtils](https://github.com/AndresMCB/AMCBGeneUtils)).

1. Please download and install Rtools from https://cran.r-project.org/bin/windows/Rtools (step only required for windows).

2. Please install devtools 

```R
install.packages("devtools")
```

3. Install DynamicCancerDriverKm package from github repository 
```R
devtools::install_github('AndresMCB/DynamicCancerDriverKM')
```

4. Please install DynamicCancerDriverKM additional packages required as follows: 
```R
devtools::install_github('AndresMCB/AMCBGeneUtils')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phenopath")

```
## Documentation 
Detailed information about the functions implemented can be found in the [user manual]()

Please find the datasets employed in our paper in the folder [data](data/)
