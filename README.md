# DynamicCancerDriver
DynamicCancerDriver package contains functions to identify genes driving one or more bio-pathological transitions during cancer progression. Formally, we name a gene driving one or more core processes over cancer progression as dynamic cancer driver. 

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
