#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library(argparse))
suppressMessages(library(dplyr))
suppressMessages(library(ComplexHeatmap))

parser <- ArgumentParser()
parser$add_argument("-c", "--chr", type="character", default=1, 
                    help="Chromosome number: 1-22 or X; Default: 1",
                    metavar="chromosome")
parser$add_argument("-s", "--start", type="integer", default=10000,
                    help="Coordinate for region start site; Default: 10000",
                    metavar="startcoord")
parser$add_argument("-e", "--end", type="integer", default=20000,
                    help="Coordinate for region stop site; Default: 20000",
                    metavar="endcoord")
parser$add_argument("-f", "--subset", type="character", default="All",
                    help="Use All samples or only Female samples; Default: All")
parser$add_argument("-t", "--category", type="character", default="All", nargs='+',
                    help="Use all samples or disease subtype(s): choose from Diffuse/Limited/Healthy; Default: All")
parser$add_argument("-m", "--diseaseCoverage", type="integer", default=0,
                    help="Allowed maximum missingness in cases; Default: 0",
                    metavar="allowedDiseaseMissing")
parser$add_argument("-n", "--controlCoverage", type="integer", default=0,
                    help="Allowed maximum missingness in controls; Default: 0",
                    metavar="allowedControlMissing")
parser$add_argument("-r", "--regressOut", type="character", nargs='+',
                    help="Regress out one or more fixed effects: choose from Age/Male/Duration/Smoke/Ethnicity/Status where Duration represents the progression duration since SSc onset, Smoke represents smoking history (0/1); Default: null")
parser$add_argument("-o", "--output", type="character", default="Pattern",
                    help="Output pdf file prefix; Default output: Pattern.pdf")

args <- parser$parse_args()

args$chr -> chr
args$start -> start
args$end -> end
args$subset -> subset
args$category -> category
args$diseaseCoverage -> diseaseCoverage
args$controlCoverage -> controlCoverage
args$regressOut -> regressOut
args$output -> filename

cat("The following arguments were used to generate methylation pattern:", collapse=("\n"))
cat("Chromosome: ", chr, collapse=("\n"))
cat("Start position: ", start, collapse=("\n"))
cat("End position: ", end, collapse=("\n"))
cat("Subset: ", subset, collapse=("\n"))
cat("Disease category: ", category, collapse=("\n"))
cat("Allowed maximum missingness in cases: ", diseaseCoverage, collapse=("\n"))
cat("Allowed maximum missingness in controls: ", controlCoverage, collapse=("\n"))
cat("Fixed effect(s) to be regressed out: ", regressOut, collapse=("\n"))
cat("Output: ", filename, ".pdf", collapse=("\n"))

imposeFilter <- function(vec, range1, range2, thres1, thres2){
  Na1 <- length(which(vec[range1]=="NaN"))
  Na2 <- length(which(vec[range2]=="NaN"))
  if (Na1 <= thres1 & Na2 <= thres2) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

linearRegression <- function(response, predictor) {
  NAcoord <- which(response=="NaN")
  response[NAcoord] <- NaN
  predictor$response <- response
  linearMod <- lm(response ~., data = predictor, na.action = na.exclude)
  residual <- residuals(linearMod)
  return(residual)
}

### Internal test
#chr <- 1
#start <- 150000
#end <- 2000000
#subset <- "Female"
#category <- c("Diffuse", "Limited")
#diseaseCoverage <- 3
#controlCoverage <- 1
#regressOut <- c("Age","Ethnicity","Male")

if ("Ethnicity" %in% regressOut) {
  regressOut <- c(regressOut,"Caucasian","Asian")
}

load("Covariates.RData")
load(paste0("MethylationBetaValues_chr",chr,".RData"))

if (chr == "X" & subset == "All") {
  subset == "Female"
  cat("Only display females for X chromosome.", collapse=("\n"))
}

if (subset == "Female") {
  datMeth <- datMeth[,covar$Male == 0]
  covar <- covar[covar$Male == 0,]
}

if (category != "All") {
  datMeth <- datMeth[,covar$Status %in% category]
  covar <- covar[covar$Status %in% category,]
}

pos <- unlist(strsplit(rownames(datMeth),"-"))[seq(2,length(unlist(strsplit(rownames(datMeth),"-"))),2)]
rownames(datMeth) <- pos

subMeth <- datMeth[as.numeric(rownames(datMeth)) >= start & as.numeric(rownames(datMeth)) <= end,]

subMeth <- subMeth[apply(subMeth, 1, function(x) imposeFilter(x,which(covar$Disease==1),which(covar$Disease==0),diseaseCoverage,controlCoverage)),]
if (is.null(nrow(subMeth))) {
  covariate <- covar[,which(colnames(covar) %in% regressOut)]
  covariate <- as.data.frame(covariate)
  cat("Only one CpG dinucleotide passed filtering. Plotted without clustering.", collapse=("\n"))
  regressResidual <- linearRegression(subMeth,covariate)
  regressResidual <- t(regressResidual)
  colnames(regressResidual) <- names(subMeth)
  if (ncol(covariate)==0) {
    regressResidual <- subMeth
  }
  pdf(paste0(filename,".pdf"))
  Heatmap(as.matrix(regressResidual), cluster_rows = F,
          cluster_columns = F,
          show_row_names = F,
          column_title = paste0("Chromosome ", chr, ": ", start, "-", end),
          name = "Residual methylation level")
  dev.off()
}
if (!is.null(nrow(subMeth))) {
  if (nrow(subMeth)==0) {
    cat("No CpG passed filtering. Consider alternative cutoff or region.", collapse=("\n"))
  }
  if (nrow(subMeth) > 1) {
    covariate <- covar[,which(colnames(covar) %in% regressOut)]
    covariate <- as.data.frame(covariate)
    cat(paste0(nrow(subMeth)," CpG dinucleotide(s) passed filtering."))
    regressResidual <- apply(subMeth, 1, function(x) linearRegression(x,covariate))
    regressResidual <- t(regressResidual)
    colnames(regressResidual) <- colnames(subMeth)
    if (ncol(covariate)==0) {
      regressResidual <- subMeth
    }
    plt <- tryCatch(
      Heatmap(as.matrix(regressResidual), cluster_rows = F,
              show_row_names = F,
              column_title = paste0("Chromosome ", chr, ": ", start, "-", end),
              name = "Residual methylation level"),
      error = function(err) {
        cat("Hierarchical clustering failed due to large volume of NAs. Plotted without clustering.", collapse=("\n"))
        return(Heatmap(as.matrix(regressResidual), cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        column_title = paste0("Chromosome ", chr, ": ", start, "-", end),
        name = "Residual methylation level"))
        })
    pdf(paste0(filename,".pdf"))
    print(plt)
    dev.off()
  }
}
