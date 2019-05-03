# SclerodermaMethylation

We conducted a pilot study on DNA methylation in Systemic Sclerosis (SSc). In this cohort, we recruited nine SSc cases and nine controls. We performed whole-genome bisulfite sequencing (WGBS) and profiled genome-wide DNA cytosine methylation landscape for each sample.

In this repository, we harbor de-identified CpG methylation beta values of all samples genome-wide (MethylationBetaValues_chr*.RData). We separated the large dataset into per-chromosome subsets. Each of the subset contains information of all valid CpG dinucleotides in all samples passing quality control of read depth on the corresponding chromosome. Due to the excessively large size of CHG/CHH methylation data, we do not store them here. However, requisition for academic purpose is possible by contacting us.

We also harbor the de-identified demographic features (Covariates.RData). Samples are in the same order as all CpG methylation data.

We provide an R program which enables visualization of methylation pattern across samples in any user-specified regions (methylationShowcase.R). This program requires a recent version of R and R packages `argparse`, `dplyr` and `ComplexHeatmap`.

After downloading the script, demographic covariates and WGBS data of chromosomes of interest, or cloning this repository locally, users can retrieve usage of this program by executing from command line: 

    Rscript methylationShowcase.R -h
```ruby
usage: ./methylationShowcase.R [-h] [-c chromosome] [-s startcoord]
                               [-e endcoord] [-f SUBSET]
                               [-t CATEGORY [CATEGORY ...]]
                               [-m allowedDiseaseMissing]
                               [-n allowedControlMissing]
                               [-r REGRESSOUT [REGRESSOUT ...]] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -c chromosome, --chr chromosome
                        Chromosome number: 1-22 or X; Default: 1
  -s startcoord, --start startcoord
                        Coordinate for region start site; Default: 10000
  -e endcoord, --end endcoord
                        Coordinate for region stop site; Default: 20000
  -f SUBSET, --subset SUBSET
                        Use All samples or only Female samples; Default: All
  -t CATEGORY [CATEGORY ...], --category CATEGORY [CATEGORY ...]
                        Use all samples or disease subtype(s): choose from
                        Diffuse/Limited/Healthy; Default: All
  -m allowedDiseaseMissing, --diseaseCoverage allowedDiseaseMissing
                        Allowed maximum missingness in cases; Default: 0
  -n allowedControlMissing, --controlCoverage allowedControlMissing
                        Allowed maximum missingness in controls; Default: 0
  -r REGRESSOUT [REGRESSOUT ...], --regressOut REGRESSOUT [REGRESSOUT ...]
                        Regress out one or more fixed effects: choose from
                        Age/Male/Duration/Smoke/Ethnicity/Status where
                        Duration represents the progression duration since SSc
                        onset, Smoke represents smoking history (0/1);
                        Default: null
  -o OUTPUT, --output OUTPUT
                        Output pdf file prefix; Default output: Pattern.pdf
 ```

Users may refer to differentially methylated regions (DMRs) and visualize by supplying genomic coordinates accordingly. For example, a DMR was identified in gene *FNB3* with range chr19: 8137960 - 8138105. This DMR can then be displayed by executing from command line:

    Rscript ./methylationShowcase.R -c 19 -s 8137960 -e 8138105 -f Female -m 3 -n 1 -r Age -o Example

    # When making multiple selections, arguments should be separated by space. E.g. ... -t Diffuse Limited -r Age Male Ethnicity ... 

![ExamplePlot](Example.pdf)

In this plot, each column represents one sample and each row represents one CpG dinucleotide within the provided boundary (inclusive). Samples are ordered based on hierarchical clustering. CpG dinucleotides are in their original sequential order. Fixed effect (of age) has been regressed out from the raw methylation beta values, leaving the residual methylation levels colored in a gradient from blue (relatively hypomethylated) to red (relatively hypermethylated). CpG dinucleotides with low read depth are colored grey. Note that if no fixed effect were to be regressed out, the residual methylation levels would remain the same as the initial methylation beta values.

For more details, please refer to our publication.

Besides, we also provide an R script (bumphunterAnalysis.R) for reproducing our results presented in the publication using the R package `bumphunter`. This script analyzes only female samples (nine cases vs. four controls) and adjusts for age, while modifications should be easily achievable to incorporate more covariates and/or compare between different SSc subtypes. This script can be executed in RStudio and the results would be stored in a list object `GenomeSSCres` whose 23 factor levels correspond to the 23 (1-22 and X) chromosomes in order. q values and averaged methylation level difference can be obtained subsequently.
