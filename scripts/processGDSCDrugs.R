library(PharmacoGx)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")
processed_dir <- paste0(args[1], "processed")

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/processed"

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

myfn <- file.path(download_dir, "screened_compounds_rel_8.0.csv")
myfn2 <- file.path(download_dir, "screened_compunds_rel_8.2.csv")

# require(gdata)
drug.info <- read.csv(myfn)
# drug.info <- drug.info[-nrow(drug.info),]
## Last row is a total summation row

drug.all <- read.csv(file.path(download_dir, "drugs_with_ids.csv"), na.strings = c("", " ", "NA"))

# ver 8.0 drug annotations (July 2019)
drug.info$unique.drugid <- matchToIDTable(ids = drug.info[, "DRUG_NAME"], tbl = drug.all, column = "GDSC2019.drugid", returnColumn = "unique.drugid")
save(drug.info, file = file.path(processed_dir, "drugInfo_8.0.RData"))

# ver 8.2 drug annotations (Feb 2020)
drug.info <- read.csv(myfn2)
drug.info$unique.drugid <- matchToIDTable(ids = drug.info[, "DRUG_NAME"], tbl = drug.all, column = "GDSC2019.drugid", returnColumn = "unique.drugid")
save(drug.info, file = file.path(processed_dir, "drugInfo_8.2.RData"))
