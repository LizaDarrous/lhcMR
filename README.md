# lhcMR


:grey\_exclamation: The `lhcMR` R package is a work in progress :grey\_exclamation:

## Overview

`lhcMR` is an R package that performs bi-directional causal estimation between a pair of traits, while accounting for the presence of a potential heritable confounder acting on the pair. The method termed `LHC-MR` aims to overcome some limitations seen in standard two-sample Mendelian Randomisation (MR) methods such as potential sample overlap, under-exploitation of genome-wide markers, and sensitivity to the presence of a heritable confounder of the exposure-outcome relationship. 
To do so, our approach builds on the typical MR framework to incorporate the presence of latent heritable confounder in a structural equation model (SEM). 
Given starting points for the parameters we wish to estimate (including bidirectional causal effect, confounder effects, direct heritabilities and more), LHC-MR optimizes the likelihood function to obtain maximum likelihood estimates for the parameters. LHC-MR then runs a block jackknife approach in order to calculate the standard error of the estimated parameters.  

The data munging section of `lhcMR` is heavily inspired by the [`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM/) R package, and uses it as well to calculate the cross-trait LDSC and obtain individual trait intercepts. `lhcMR` also uses the [`TwoSampleMR`](https://github.com/MRCIEU/TwoSampleMR/) R package for inverse-variance weighted (IVW) MR causal effect estimation. 

## Installation

You can install the current version of `lhcMR` using the following code snippet:    
```r
#install.packages("devtools")
devtools::install_github("LizaDarrous/lhcMR")
library(lhcMR)
```

## Usage

#### There are several input files (in specific formats) needed to run `lhcMR`:

##### 1. GWAS summary statistics files of the two traits (exposure & outcome):
These files can be regular tab/space/comma- separated files or gzipped files [.tsv.bgz /.gz]), and they should be read into `R` as data frames. They must contain the following columns in these possible variations of their names (not case sensitive):

* SNP: `SNP`,`SNPID`,`RS_NUMBER`,`RS_NUMBERS`, `MARKERNAME`, `ID`,`PREDICTOR`,`SNP_ID`,`RS`
* Effect-size: `B`,`BETA`,`LOG_ODDS`,`EFFECTS`,`EFFECT`,`SIGNED_SUMSTAT`,`EST`, `BETA1`, `LOGOR`,`BETA_WF`
* Standard error: `STANDARD_ERROR`,`STD_ERR`,`STDERR`,`SE_BETA_WF`,`SE_BETA`
* T-statistic: `TSTAT`,`Z`,`ZSCORE`,`ZSTAT`,`ZSTATISTIC`
    If the T-statistic column is not present, it will be automatically calculated from the Effect-size and Standard error columns as such : `Effect-size/Standard error`  
* P-value: `P`,`PVALUE`,`PVAL`,`P_VALUE`,`P-VALUE`,`P.VALUE`,`P_VAL`,`GC_PVALUE`,`WALD_P`,`P_BETA_WF`
    If the P-value column is not present, it will be automatically calculated from the Effect-size and Standard error columns as such : `2*pnorm(-abs(T-statistic))`
* Alternate (Effect) allele:  `A1`, `ALLELE1`,`EFFECT_ALLELE`,`ALT`,`EA`    
* Reference allele: `A2`,`ALLELE2`,`ALLELE0`,`OTHER_ALLELE`,`REF`,`NON_EFFECT_ALLELE`,`DEC_ALLELE`,`OA`,`NEA`    
* Sample size: `N`,`NCOMPLETESAMPLES`, `TOTALSAMPLESIZE`, `TOTALN`, `TOTAL_N`,`N_COMPLETE_SAMPLES`, `SAMPLESIZE`,`N_REG`  
* Chromosome: `HG18CHR`,`CHR`, `CHROM`, `CHROMOSOME`    
* Position: `BP`,`POS`,`POSITION`   

##### 2. LD files paths (ld & rho):
The paths for two files should be provided first in the data merging file to select a subset of SNPs with LD information and create the data frame used in the actual analysis. These files are:

- LDfile: 

> To generate our own LDscore, we first took 4,773,627 SNPs with info (imputation certainty measure) ≥ 0.99 present in the association summary files from the second round of GWAS by the Neale lab. This set was restricted to 4,650,107 common, high-quality SNPs, defined as being present in both UK10K and UK Biobank, having MAF > 1% in both data sets, non-significant (Pdiff > 0.05) allele frequency difference between UK Biobank and UK10K and residing outside the HLA region (chr6:28.5-33.5Mb). For these SNPs, LD scores and regression weights were computed based on 3,781 individuals from the UK10K study. These LDscores can be downloaded [here](https://drive.google.com/file/d/1_B8DsGejRpeiW3XxNERWpmSQH2a8NWqS/view?usp=sharing). Users can also use the classical LD scores from the [Broad Insititue](https://data.broadinstitute.org/alkesgroup/LDSCORE/).

- Local LD pattern (`Rho`) for each SNP

> To estimate the local LD distribution for each SNP (k), characterised by (`π_k`,`σ^2_k`), we fitted a two-component Gaussian mixture distribution to the observed local correlations (focal SNP +/− 2’500 markers with MAF≥ 0.5% in the UK10K): (1) one Gaussian component corresponding to zero correlations, reflecting only measurement noise (whose variance is proportional to the inverse of the reference panel size) and (2) a second component with zero mean and a larger variance than the first component (encompassing measurement noise plus non-zero LD). These local LDscores can be downloaded [here](https://drive.google.com/file/d/1PwsJ2KiCR2VDVOodEXnQu6Tec54a0HC1/view?usp=sharing).


##### 3. Input files for LDSC (ld & hm3):
These are needed by the [`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM/) R-package to run the `GenomicSEM::ldsc` function. Their description is quoted from the manual:   

- hm3: 

> The name of the reference file. Here we use Hapmap 3 SNPs. This file can be obtained from https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v.

- ld:   

> A folder of LD scores used as the independent variable in LDSC (ld) and LDSC weights (wld). These are typically the same folder, and in the original LD score package is called "eur_w_ld_chr". We use the same LD scores and weights for our application, though the user can supply their own if desired. Weights for the european population used here can be obtained by downloading the `eur_w_ld_chr` folder in the link below (Note that these are the same weights provided by the original developers of LDSC): https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v 

* * *

#### There are 3 main functions in `lhcMR`: 

-   `merge_sumstats()`  
reads in summary statistics of two traits as well as file paths of LD scores and merges the data into a single data frame to be used as input for later functions. During the merge, this function harmonises the SNPs and filters out any that fall in the HLA region.
-   `calculate_SP()` 
uses the previously generated data frame to smartly generate starting points (SP) that will be used for the parameter estimation in the pair trait analysis done in the next step. 
-   `lhc_mr()` 
main function that uses the input data frame and the stating points to optimize the likelihood function and estimate the parameters (most notably the bidirectional causal effect, confounder effect and total trait heritability), as well as their standard error (SE) using block jackknife.

On overview of the usage of each is provided in the [manual](doc/lhcMR_manual.pdf), and their description and input parameters are detailed in the steps below:  
**Note1**: we advise to run the analysis on a master node instead of a worker node, especially if you choose to run the parallelisation using the `lapply` option in Step 3 in order to have access to cores.  
**Note2**: if running the analysis on a worker node with large summary statistics files (such as UKBiobank), it's best to allocate at least 15GB of memory.


##### Step 1: Reading in and merging Data using `merge_sumstats()`:
`merge_sumstats()` takes in a `input.files` list containing the two summary statistics files with all their required columns detailed above, a vector `trait.names` containing two strings corresponding to the trait names in the order they were placed in the list. As well as the **paths** for the LDscore files mentioned above in `LD.filepath` and `rho.filepath`. This function first checks for the presence of the needed columns in all files, then joins the different files based on chromosome and position, checks for swapped alleles and corrects their effects and finally returns a large dataframe with all the columns needed to be used in all the subsequent steps.

##### Step 2: Calculating smart starting points for the likelihood optimisation using `calculate_SP()`:
Taking in the resulting data frame from the previous method as `input.df`, as well as the vector `trait.names`, `calculate_SP()` first aims to generate starting points for the `iXY, axy, ayx` parameters, which are the cross-trait intercept and the bidirectional causal effect between the exposure `X` and the outcome `Y`. If `run_ldsc=TRUE`, then `iXY` will be calculated using the `GenomicSEM::ldsc()` function, where the input paths `hm3` and `ld` are needed, otherwise a random value will be generated instead. Similarly, if `run_MR=TRUE`, then both `axy` and `ayx` will be calculated from IVW-MR using `TwoSampleMR::mr()`, otherwise random values will be generated.

After obtaining these values to use as starting points for their parameter estimation, others can be obtained from a random distribution or from a single trait analysis. Parameters like the confounding effect `tX, tY` and direct trait heritabilities `h2X, h2Y` are obtained randomly. Whereas parameters such as the trait intercepts or SNP polygenicity are obtained from the single trait analysis where `LHC-MR` is run as a simplified model given only the SNP effects of a single trait, without any input from a second trait or a confounder.

To speed up the single trait (or pair-trait analysis in step 3), the data input is reduced by every n-th SNP using the `SNP_filter` parameter (default=10). The number of random starting points sets for the single trait analysis should be set using the `SP_single` parameter, we advise a value between 3-5.  
`nStep` is the parameter used to indicate how many steps should the `LHC-MR` analysis take. Values can be either `nStep=1`, where the analysis is done such that all 9 parameters are estimated in the second step simultaneously and only `iX,iY` are fixed, **or** `nStep=2` where 7 parameters are estimated and `piX,piY` are fixed in addition to `iX,iY` from the single-trait analysis.

`SP_pair` indicates how many sets of starting points should be generated from all these parameters (generated with small noise added or randomly). If `nStep=1`, we advise to have more starting points `80-100`, whereas if `nStep=2` the value can be 50 or higher.
`saveRFiles=TRUE` will write the results of the previous functions as `.csv` files in your working directory.  
Note1: if you are using your own LDscore, please specify the number of SNPs used to generate the file in `M`.  
Note2: if you do not wish to parallelise over cores, set `nCores=1`.  

This function returns a list containing `iX,iY,piX,piY` to be fixed in the optimisation of the next step, as well as the n-th SNP filtered data frame `input.df_filtered` and the generated starting points matrix `sp_mat`. On top of that, the values of the IVW estimated bidirectional effect `axy_MR,ayx_MR` and the cross-trait intercept `iXY` are reported. 

##### Step 3: Running the likelihood optimisation to estimate the parameters, followed by a block-jackknife procedure to calculate parameter-SE using `lhc_mr()`:
`lhc_mr()` takes the input list from the previous functions, the vector of `trait.names` to run the pair trait optimisation. The parallelisation of the optimisation can be done in two ways, depending on the `paral_method` chosen. If `paral_method="rslurm"`, then the Slurm workload manager will parallelise the optimisation over the sets of starting points (faster), and the partition you wish to use should be specified in `partition` as a string.  
Otherwise if `paral_method="lapply"` then the `parallel::mclapply` function is used to split the parallelisation over the available cores. This way is error-prone due to the availablity or lack thereof of cores to be used.  
`nBlock` indicates how many blocks should be created for the block jackknife analysis, by default this value is `200`.

The results of the optimisation over the starting points, as well as the block jackknife optimisation are written as `.csv` files. The summary calculations of the block jackknife procedure is also reported as a `.csv` file, and finally the parameter estimates, their SE and p-values are printed out and reported as a `.csv` file. 

**Note1**: if you are using your own LDscore, please specify the number of SNPs used to generate the file in `M`.  
**Note2:** if you do not wish to parallelise over cores, set `nCores=1`.  

## Example script
```
library(lhcMR)

## Set up a working directory where all the analysis files will be printed
setwd("~/BWeightDM/") 

## Read in summary statistics file, here we use Neale UKBB data 
X = data.table::fread(cmd = "zcat < ~/data/20022_irnt.gwas.imputed_v3.both_sexes.tsv.bgz")
Y = data.table::fread(cmd = "zcat < ~/data/2443.gwas.imputed_v3.both_sexes.tsv.bgz")

## To make sure we have all the columns needed (including allele data), we merge with Neale-provided variants file
VAR = data.table::fread(cmd = "zcat < ~/data/variants.tsv.bgz")
VAR$chr=as.numeric(VAR$chr)
X = dplyr::inner_join(X,VAR[,c(1:6)])
Y = dplyr::inner_join(Y,VAR[,c(1:6)])

## File paths needed for the analysis
LD.filepath = "~/data/LDscores_filtered.csv" # LD scores
rho.filepath = "~/data/LD_GM2_2prm.csv" # local/SNP-specfic LD scores

ld = "~/genomicSEM/data/eur_w_ld_chr/"
hm3 = "~/genomicSEM/data/w_hm3.noMHC.snplist"

rm(VAR) # to free up space

## Step 1
trait.names=c("BWeight","DM")
input.files = list(X,Y)
df = merge_sumstats(input.files,trait.names,LD.filepath,rho.filepath)

## Step 2
SP_list = calculate_SP(df,trait.names,run_ldsc=TRUE,run_MR=TRUE,hm3=hm3,ld=ld,nStep = 2,
                       SP_single=3,SP_pair=50,SNP_filter=10)

## Step 3
res = lhc_mr(SP_list, trait.names, paral_method="rslurm", nBlock=200)

```

## Citation
Our manuscript is under review, but you can read our pre-print [here](https://www.medrxiv.org/content/10.1101/2020.01.27.20018929v3).

## Contact
Liza Darrous <darrous.liza@gmail.com>
