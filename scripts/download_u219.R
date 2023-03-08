library("downloader")

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/download"

tmpdir <- file.path(download_dir, "dwl")
dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(download_dir, "array"), showWarnings = FALSE, recursive = TRUE)

ftpdir <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3610/"

## phenodata
dwl.status <- download(url = sprintf("%s/E-MTAB-3610.sdrf.txt", ftpdir), destfile = file.path(download_dir, "E-MTAB-3610.sdrf.txt"), quiet = TRUE)
sampleinfo <- read.csv(file.path(download_dir, "E-MTAB-3610.sdrf.txt"), sep = "\t", stringsAsFactors = FALSE)
rownames(sampleinfo) <- sampleinfo[, "Assay.Name"]
sampleinfo[, "Array.Data.File"] <- gsub("[.]cel$", ".CEL.gz", sampleinfo[, "Array.Data.File"])
uarchive <- sort(unique(sampleinfo[, "Comment..ArrayExpress.FTP.file."]))


message("Download genomic data")

require(R.utils) || stop("Library R.utils is not available!")

## download and compress CEL files
celfile.timestamp <- celfn <- NULL
i <- 1
while (i <= length(uarchive)) {
  ## assuming there are only 9 zip archives (need to check if the update version has more)
  dwl.status <- download(url = uarchive[i], destfile = file.path(tmpdir, basename(uarchive)[i]), quiet = TRUE)
  if (dwl.status != 0) {
    message("\t-> download failed, let's try again ...")
    file.remove(file.path(tmpdir, basename(uarchive)[i]))
    i <- i - 1
  } else {
    ## unzip archive
    fff <- unzip(zipfile = file.path(tmpdir, basename(uarchive)[i]), list = TRUE)
    celfile.timestamp <- c(celfile.timestamp, as.character(fff[, "Date"]))
    celfn <- c(celfn, as.character(fff[, "Name"]))
    res <- unzip(zipfile = file.path(tmpdir, basename(uarchive)[i]), exdir = file.path(download_dir, "array"))
    ## compress each CEL file individually using gzip
    library(R.utils)
    sapply(file.path(download_dir, "array", as.character(fff[, "Name"])), R.utils::gzip, overwrite = TRUE)
    i <- i + 1
  }
}
celfile.timestamp <- t(sapply(strsplit(celfile.timestamp, split = " "), function(x) {
  return(x)
}))
dimnames(celfile.timestamp) <- list(gsub("[.]CEL$", "", celfn), c("file.day", "file.hour"))
write.csv(celfile.timestamp, file = file.path(download_dir, "celfile_timestamp_u219.csv"))
save(list = c("celfile.timestamp"), compress = TRUE, file = file.path(download_dir, "celfile_timestamp_u219.RData"))

ftpdir <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3610/"
## phenodata
dwl.status <- download(url = sprintf("%s/E-MTAB-3610.sdrf.txt", ftpdir), destfile = file.path(download_dir, "E-MTAB-3610.sdrf.txt"), quiet = TRUE)

zip(zipfile = file.path(download_dir, "gdsc_array_u219"), files = list.files(path = file.path(download_dir, "array"), full.names = TRUE), extras = "-j")

unlink(file.path(download_dir, "dwl"), recursive = T)
unlink(file.path(download_dir, "array"), recursive = T)
