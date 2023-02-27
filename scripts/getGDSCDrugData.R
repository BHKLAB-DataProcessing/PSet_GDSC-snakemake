library(PharmacoGx)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")
processed_dir <- paste0(args[1], "processed")

download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/download"
processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/processed"

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


# require(gdata)
drug.info <- read.csv(myfn)
# drug.info <- drug.info[-nrow(drug.info),]
## Last row is a total summation row

drug.all <- read.csv(file.path(download_dir, "drugs_with_ids.csv"), na.strings = c("", " ", "NA"))

# ver 8.0 drug annotations (July 2019)
drug.info$unique.drugid <- matchToIDTable(ids = drug.info[, "DRUG_NAME"], tbl = drug.all, column = "GDSC2019.drugid", returnColumn = "unique.drugid")
save(drug.info, file = file.path(processed_dir, "drugInfo.RData"))

# ver 8.2 drug annotations (Feb 2020)
drug.info <- read.csv(myfn2)
drug.info$unique.drugid <- matchToIDTable(ids = drug.info[, "DRUG_NAME"], tbl = drug.all, column = "GDSC2019.drugid", returnColumn = "unique.drugid")
save(drug.info, file = file.path(processed_dir, "drugInfo_8.2.RData"))

# downloadGDSCdrugs <- function(path.data="/pfs/out",  path.drug = path.data){
# 	if(!file.exists(path.drug)){
# 		dir.create(path.drug)
# 	}
# 	## download drug sensitivity
# 	message("Download drug sensitivity measurements")
# 	myfn <- file.path(path.drug, "gdsc_drug_sensitivity.csv")
# 	if (!file.exists(myfn)) {
# 		dir.create(file.path(path.drug, "dwl"), showWarnings=FALSE, recursive=TRUE)
# 		dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_manova_input_w5.csv", destfile=file.path(path.drug, "dwl", "gdsc_manova_input_w5.csv"), quiet=TRUE)
# 		if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
# 		file.copy(from=file.path(path.drug, "dwl", "gdsc_manova_input_w5.csv"), to=myfn)
# 	}

# 	## download drug concentration
# 	message("Download screening drug concentrations")
# 	myfn <- file.path(path.drug, "gdsc_drug_concentration.csv")
# 	if (!file.exists(myfn)) {
# 		dir.create(file.path(path.drug, "dwl"), showWarnings=FALSE, recursive=TRUE)
# 		dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_compounds_conc_w5.csv", destfile=file.path(path.drug, "dwl", "gdsc_compounds_conc_w5.csv"), quiet=TRUE)
# 		if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
# 		file.copy(from=file.path(path.drug, "dwl", "gdsc_compounds_conc_w5.csv"), to=myfn)
# 	}

# 	## download drug information
# 	message("Download drug information")
# 	myfn <- file.path(path.drug, "gdsc_drug_information.csv")
# 	if (!file.exists(myfn)) {
# 		dir.create(file.path(path.drug, "dwl"), showWarnings=FALSE, recursive=TRUE)
# 	  # dwl.status <- download.file(url="http://www.cancerrxgene.org/action/ExportJsonTable/CSV", destfile=file.path(path.drug, "dwl", "export-Automatically_generated_table_data.csv"), quiet=TRUE)
# 	  # if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
# 	  drugs <- read.csv("https://www.cancerrxgene.org/translation/drug_list?list=all&export=csv")
# 	  # drugs <- tables[1][[1]]
# 	  write.csv(drugs, row.names=FALSE, file=file.path(path.drug, "dwl", "export.csv"))
# 	  file.copy(from=file.path(path.drug, "dwl", "export.csv"), to=myfn)
# 	}
# 	myfn <- file.path(path.drug, "nature_supplementary_information.xlsx")
# 	if (!file.exists(myfn)) {
# 		dir.create(file.path(path.drug, "dwl"), showWarnings=FALSE, recursive=TRUE)
# 		dwl.status <- download.file(url="http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip", destfile=file.path(path.drug, "dwl", "nature11005-s2.zip"), quiet=TRUE)
# 		ff <- as.character(unzip(zipfile=file.path(path.drug, "dwl", "nature11005-s2.zip"), list=TRUE)[1, 1])
# 		unzip(zipfile=file.path(path.drug, "dwl", "nature11005-s2.zip"), exdir=file.path(path.drug, "dwl"))
# 		file.copy(from=file.path(path.drug, "dwl", ff), to=myfn)
# 	}

# }

# downloadGDSCdrugs()
