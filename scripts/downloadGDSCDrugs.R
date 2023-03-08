options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/download"

matchToIDTable <- function(ids, tbl, column, returnColumn = "unique.cellid") {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)", Hmisc::escapeRegex(x), "((///)|$)"), tbl[, column])
    if (length(myx) > 1) {
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if (length(myx) == 0) {
      return(NA_character_)
    }
    return(tbl[myx, returnColumn])
  })
}

require(downloader)
badchars <- "[\xb5]|[]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
sessionInfo()

drugFileURL <- "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.0/"
drugFileName <- "screened_compounds_rel_8.0.csv"

drugFileURL_8.2 <- "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.2/"
drugFileName_8.2 <- "screened_compunds_rel_8.2.csv"

## download sample information
message("Download drug info")
myfn <- file.path(download_dir, "screened_compounds_rel_8.0.csv")
myfn2 <- file.path(download_dir, "screened_compunds_rel_8.2.csv")

dwl.status <- download.file(url = sprintf("%s/%s", drugFileURL, drugFileName), destfile = myfn, quiet = TRUE)
if (dwl.status != 0) {
  stop("Download failed, please rerun the pipeline!")
}

dwl.status <- download.file(url = sprintf("%s/%s", drugFileURL_8.2, drugFileName_8.2), destfile = myfn2, quiet = TRUE)
if (dwl.status != 0) {
  stop("Download failed, please rerun the pipeline!")
}
