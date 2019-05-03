# SclerodermaMethylation

We conducted a pilot study on DNA methylation in Systemic Sclerosis (SSc). In this cohort, we recruited nine SSc cases and nine controls. We performed whole-genome bisulfite sequencing (WGBS) and profiled genome-wide DNA cytosine methylation landscape for each sample.

In this repository, we harbor de-identified CpG methylation beta values of all samples genome-wide (MethylationBetaValues_chr*.RData). We separated the large dataset into per-chromosome subsets. Each of the subset contains information of all valid CpG dinucleotides in all samples passing quality control of read depth on the corresponding chromosome.

We also harbor the de-identified demographic features (Covariates.RData). Samples are in the same order as all methylation profile data.

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
                        Age/Male/Duration/Smoke/Ethnicity/Status; Default:
                        null
  -o OUTPUT, --output OUTPUT
                        Output pdf file prefix; Default output: Pattern.pdf
 ```
