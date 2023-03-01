library(PharmacoGx)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
processed_dir <- paste0(args[1], "processed")
sens_version <- args[2]

processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/processed"
sens_version <- "8.2"

slices_dir <- file.path(processed_dir, paste0('slices', str_replace(sens_version, '\\.', '_')))
dir.create(slices_dir)
unzip(file.path(processed_dir, paste0("raw_sense_slices_", sens_version, ".zip")), exdir = slices_dir, junkpaths = TRUE)

files <- list.files(slices_dir, full.names = TRUE)
recomp_dir <- file.path(processed_dir, paste0('slices_recomp', str_replace(sens_version, '\\.', '_')))
dir.create(recomp_dir)
for (file in files) {
  print(file)
  mybasenm <- basename(file)
  
  slice <- readRDS(file)
  
  res <- PharmacoGx:::.calculateFromRaw(slice)
  
  saveRDS(res, file = file.path(recomp_dir, gsub(mybasenm, pattern = ".rds", replacement = "_recomp.rds", fixed = TRUE)))
}

###assemble###
myfn <- list.files(path=recomp_dir, full.names = TRUE)
slices <- list()

for(fn in myfn){
  temp <- readRDS(fn)
  parTable <- do.call(rbind,temp[[3]])
  # print(head(rownames(parTable)))
  # print(str(temp[[3]]))
  n <- cbind("aac_recomputed" = as.numeric(unlist(temp[[1]]))/100, 
             "ic50_recomputed" = as.numeric(unlist(temp[[2]])), 
             "HS" = as.numeric(unlist(parTable[,1])),
             "E_inf" = as.numeric(unlist(parTable[,2])),
             "EC50" = as.numeric(unlist(parTable[,3]))) 
  rownames(n) <- names(temp[[3]])
  slices[[fn]] <- n
}

res <- do.call(rbind, slices)

save(res, file=file.path(processed_dir, paste0("profiles", sens_version, ".RData")))

unlink(slices_dir, recursive=TRUE)
unlink(recomp_dir, recursive=TRUE)
