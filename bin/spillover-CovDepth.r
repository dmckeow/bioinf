#!/usr/bin/env Rscript
### setwd("/panfs/jay/groups/27/dcschroe/dmckeow/data/Spillover_FINAL/tmp.mapping")

### for parallel (all other packages must be loaded within the parallel loop)
library(foreach)
library(doParallel)
library(here)

# Register parallel backend
# Adjust 'n_cores' to the number of cores you want to utilize
n_cores <- 12
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# List all files you want to process
files <- list.files(path=here(), pattern=".*.samtoolsdepth", full.names=TRUE, recursive=FALSE)

# Initialize an empty list to store the outputs
output_list <- list()

# Loop over each file - this is the parallel equivalent of for loop
output_list <- foreach(file = files, .combine = rbind) %dopar% {
    library(dplyr)
    library(zoo)
    library(data.table)

    print(paste("Processing:", file))
    data <- data.table::fread(file, header = FALSE, sep = "\t")
    colnames(data) <- c("ReadSource", "Contig", "Locus", "Depth")
    
    # Perform operations on data
    MappingCov200bpWindow <- data %>% 
      group_by(Contig,ReadSource) %>% 
        do(
          data.frame(
            start = rollapply(.$Locus, width=200, by=200, FUN=min, align="left"),
            end = rollapply(.$Locus, width=200, by=200, FUN=max, align="left"),
            coverage = rollapply(.$Depth, width=200, by=200, FUN=mean, align="left")
    )
  ) %>% ungroup()
    
    # Store the output in the list
    output_list[[file]] <- MappingCov200bpWindow
}

# Stop the parallel backend
stopCluster(cl)

data.table::fwrite(output_list, file = "AllSamtoolsDepthCovMappingByWindow.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

