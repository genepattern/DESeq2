##
## Copyright (c) 2016 Broad Institute, Inc. and Massachusetts Institute of Technology.  All rights reserved.
##

GP.deseq2 <- function(gct, cls, confounding.var.cls, qc.plot.format, phenotype.test,
                      output.file.base, random.seed, fdrThreshold, topNCount) {

    # Make sure there are at least two classes for comparison
    if (NROW(cls$names) < 2) {
        stop("The cls.file must be categorical with at least two classes")
    }

    if (!is.null(random.seed) && !is.na(random.seed)) {
       set.seed(random.seed)
    }

    # Build a conditions data structure based on the CLS file(s), to be used as column metadata in the DESeqDataSet.
    # Note that our comparison order is always based on the order of classes in the CLS (e.g. cls$names)
    if (is.null(confounding.var.cls)) {
        colData <- data.frame(conditions=cls$labels)
        design = "~ conditions"
    } else if (NROW(cls$labels) != NROW(confounding.var.cls$labels)) {  # Check if the label counts for the CLS files match
        stop("Label mismatch detected between cls.file and confounding.var.cls.  Make sure both map to the same number of samples.")
    } else {
        colData <- data.frame(conditions=cls$labels, confoundingVar=confounding.var.cls$labels)
        design = "~ confoundingVar + conditions"
    }

    # Check if the CLS file maps to all columns of gct data
    if (NROW(cls$labels) != NCOL(gct$data)) {
        stop("Label mismatch detected between cls.file and gct.file.  Make sure the CLS maps to all GCT samples.")
    }
    
    rownames(colData) <- colnames(gct$data)
    
    # Filter out any low expression rows (with only 0 or 1 read) for better memory usage & speed of transformation.
    data_filtered <- round(gct$data[rowSums(gct$data) > 1,])
    
    # Build the DESeqDataSet to be used for the Differential Expression analyses.  We will specify the condition order level
    # later, at the time of each DESeq() call; we may also manipulate the colData as well, depending on the number of classes
    # and phenotype.test setting.
    # Note that we suppressMessages here as the call creates some spurious noise on stderr, causing issues in a GP setting.
    suppressMessages(
        ddsFullCountTable <- DESeqDataSetFromMatrix(countData=data_filtered, colData=colData, design=as.formula(design))
    )
    
    # Perform differential expression and build reports
    if (NROW(cls$names) == 2) {
        # With only two classes there is only one way to set up 
        ddsFullCountTable$conditions <- factor(ddsFullCountTable$conditions, levels=cls$names)
        dds <- DESeq(ddsFullCountTable, quiet=TRUE)
       
        write.reports(dds, output.file.base, cls$names[1], cls$names[2], fdrThreshold, topNCount, qc.plot.format)
    } else if (phenotype.test == "one-versus-all") {
        # We'll iterate through the the cls$names one by one, comparing samples of that class against the rest of the samples
        # (all grouped together into one "Rest" pseudo-class).
        for (currClass in cls$names) {
            # Set up a new conditions list with only this class (as the reference) and "Rest" labels for the other samples.
            newCond <- rep("Rest", NROW(cls$labels))
            classIndicies <- grep(paste0("^", currClass, "$"), cls$labels)
            newCond <- replace(newCond, classIndicies, currClass)
            colData$conditions <- factor(newCond, levels=c(currClass, "Rest"))
            ddsFullCountTable$conditions <- colData$conditions

            dds <- DESeq(ddsFullCountTable, quiet=TRUE)
            write.reports(dds, output.file.base, currClass, "Rest", fdrThreshold, topNCount, qc.plot.format)
        }
    } else { # Otherwise we're doing all-pairs
        # We'll iterate through the the cls$names one by one, comparing samples of that class to the other classes one at a time.
        for (i in 1:(NROW(cls$names) - 1) ) {  # Goes up to but does not include the last class
            # Set the currClass as the reference level
            currClass <- cls$names[i]
            ddsFullCountTable$conditions <- relevel(ddsFullCountTable$conditions, ref=currClass)

            # We only need to run the DESeq() call once per iteration to get the pairwise comparisons.
            dds <- DESeq(ddsFullCountTable, quiet=TRUE)
            for (j in (i+1):NROW(cls$names)) {
                otherClass <- cls$names[j]
                write.reports(dds, output.file.base, currClass, otherClass, fdrThreshold, topNCount, qc.plot.format)
            }
        }
    }
}

write.reports <- function(dds, output.file.base, class0, class1, fdrThreshold, topNCount, qc.plot.format) {
    # Note that this means we are comparing "class1 vs. class0" and thus class0 represents the control. 
    res <- results(dds, contrast=c("conditions", class1, class0))
    
    # Adjust NA padj values;  This is courtesy of Brian Haas' code.  DESeq2 sets these to NA if they do not pass the
    # significance threshold, but leaving them as NA can give misleading output; setting them to 1 means they will 
    # not pass the FDR threshold.
    res$padj[is.na(res$padj)]  <- 1
    
    # Based on the original code from Brian Haas' example.  Instead, we're going to give the class means in one
    # file and then the standard DESeq2 report separately
    #means_report <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
    #means_report <- cbind(id=rownames(means_report), as.data.frame(means_report))

    # Create some reports, ordered by adjusted p-value
    resOrdered <- res[order(res$padj),]

    # Write the basic results report from DESeq2
    output.file.base <- paste(output.file.base, class1, "vs", class0, sep=".")
    write.table(cbind(id=rownames(resOrdered), as.data.frame(resOrdered)), sep='\t', quote=FALSE, row.names=FALSE, 
                append=TRUE, paste0(output.file.base, ".DESeq2_results_report.txt"))

    # Extra info, inspired by Brian Haas' example code.  Write out the row-level mean values by class.
    baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == class0])
    baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == class1])
    means_report <- cbind(baseMeanA, baseMeanB)
    colnames(means_report) <- c(class0, class1)
    means_report <- cbind(id=rownames(res), means_report)
    
    # Re-order the means_report to match the basic report and then write it as well.
    means_report <- means_report[rownames(resOrdered),]
    write.table(as.data.frame(means_report), sep='\t', quote=FALSE, row.names=FALSE, append=TRUE,
               file=paste0(output.file.base, ".mean_values_by_class.txt"))
        
    # Print out some QC images. Plot nothing if the user chose "skip"
    device.open <- get.device.open(qc.plot.format)
    if (!is.null(device.open)) {
        print(paste0("Generating QC plots: ", output.file.base))
        print.MAplot(res, device.open, output.file.base)
        print.DispEsts(dds, device.open, output.file.base)
    }

    # Get the subset consisting only of the significant genes, as determined by an FDR threshold test.
    resSig <- subset(res, padj < fdrThreshold)
    
    # Check the size; if everything has been filtered out then we let the user know and skip these reports.
    if (NROW(resSig) == 0) {
       message <- paste0("Cannot report top up-regulated & down-regulated genes for ", class1, " vs. ", class0,
                    " .  No results pass the FDR filter threshold of ", fdrThreshold, ".")
       write(message, file=paste0(output.file.base, ".top", topNCount, "_upregulated_genes_report.txt"))
       write(message, file=paste0(output.file.base, ".top", topNCount, "_downregulated_genes_report.txt"))
    } else {
        # Top-N down-regulated & down-regulated genes in the dataset.
        upRegGenes <- head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ], n=topNCount)
        write.table(cbind(id=rownames(upRegGenes), as.data.frame(upRegGenes)), sep='\t', quote=FALSE, row.names=FALSE,
                   file=paste0(output.file.base, ".top", topNCount, "_upregulated_genes_report.txt"))
        downRegGenes <- head(resSig[ order(resSig$log2FoldChange), ], n=topNCount)
        write.table(cbind(id=rownames(downRegGenes), as.data.frame(downRegGenes)), sep='\t', quote=FALSE, row.names=TRUE,
                   file=paste0(output.file.base, ".top", topNCount, "_downregulated_genes_report.txt"))
    }
}

check.output.format <- function(output.format) {
   if (!(output.format %in% c("pdf", "svg", "png", "skip"))) {
      stop(paste0("Unrecognized output format '", output.format, "'"))
   }
}

check.phenotype.test <- function(phenotype.test) {
   if (!(phenotype.test %in% c("one-versus-all", "all-pairs"))) {
      stop(paste0("Unrecognized phenotype.test '", phenotype.test, "'"))
   }
}

get.device.open <- function(extension) {
   if (extension == "skip") { return(NULL) }
   if (extension == "pdf") {
      return(function(filename_base) {
         pdf(paste0(filename_base, ".pdf"))
      })
   }
   if (extension == "svg") {
      return(function(filename_base) {
         svg(paste0(filename_base, ".svg"))
      })
   }
   if (extension == "png") {
      return(function(filename_base) {
         png(paste0(filename_base, ".png"))
      })
   }
   stop(paste0("Unhandled plot file format '", extension, "'"))
}

print.MAplot <- function(res, device.open, output.file.base) {
   plotname <- paste0(output.file.base, ".QC.MAplot")
   tryCatch({
      device.open(plotname)
      plotMA(res, ylim=c(-5,5))
      dev.off()
   },
   error = function(err) {
      print(paste0("Error printing the ", plotname, " plot - skipping"))
      print(conditionMessage(err))
   })
}

print.DispEsts <- function(dds, device.open, output.file.base) {
   plotname <- paste0(output.file.base, ".QC.DispEsts")
   tryCatch({
      device.open(plotname)
      plotDispEsts(dds)
      dev.off()
   },
   error = function(err) {
      print(paste0("Error printing the ", plotname, " plot - skipping"))
      print(conditionMessage(err))
   })
}