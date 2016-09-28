##
## Copyright (c) 2016 Broad Institute, Inc. and Massachusetts Institute of Technology.  All rights reserved.
##

suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(DESeq2)))

sessionInfo()

arguments <- commandArgs(trailingOnly=TRUE)

libdir <- arguments[1]

option_list <- list(
  make_option("--input.file", dest="input.file"),
  make_option("--cls.file", dest="cls.file"),
  make_option("--confounding.variable.cls.file", dest="confounding.variable.cls.file", default=NULL),
  make_option("--output.file.base", dest="output.file.base"),
  make_option("--qc.plot.format", dest="qc.plot.format"),
# Hiding advanced feature: not allowing >2 classes in Beta
#  make_option("--phenotype.test", dest="phenotype.test"),
  make_option("--fdrThreshold", dest="fdr.threshold", type="numeric", default=0.1),
  make_option("--topNCount", dest="top.N.count", type="integer", default=20),
  make_option("--randomSeed", dest="random.seed", type="integer", default=779948241)
  )

opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=arguments)
print(opt)
opts <- opt$options
source(file.path(libdir, "common.R"))
source(file.path(libdir, "gp_deseq2.R"))

check.output.format(opts$qc.plot.format)
# Hiding advanced feature: not allowing >2 classes in Beta
#check.phenotype.test(opts$phenotype.test)
#phenotype.test -> opts$phenotype.test
phenotype.test -> NULL

# Load the GCT and CLS.
gct <- read.gct(opts$input.file)
cls <- read.cls(opts$cls.file)
if (!is.null(opts$confounding.variable.cls.file)) {
   confounding.var.cls <- read.cls(opts$confounding.variable.cls.file)
} else {
   confounding.var.cls <- NULL
}

GP.deseq2(gct, cls, confounding.var.cls, opts$qc.plot.format, 
          phenotype.test, opts$output.file.base,
          opts$random.seed, opts$fdr.threshold, opts$top.N.count)

sessionInfo()