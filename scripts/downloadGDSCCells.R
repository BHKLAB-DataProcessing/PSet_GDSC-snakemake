# library("readxl")
require(downloader)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")
# processsed_dir <- paste0(args[1], "processed")

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/processed"

badchars <- "[\xb5]|[]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
sessionInfo()

cellFileURL <- "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.0/"
cellFileName <- "Cell_Lines_Details.xlsx"

cellFileURL_8.2 <- "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.2/"
cellFileName_8.2 <- "Cell_Lines_Details.xlsx"

## download sample information
message("Download cell info")
myfn <- file.path(download_dir, "gdsc_cellinfo.xlsx")
myfn2 <- file.path(download_dir, "gdsc_cellinfo_8.2.xlsx")

dwl.status <- download.file(url = sprintf("%s/%s", cellFileURL, cellFileName), destfile = myfn, quiet = TRUE)
if (dwl.status != 0) {
  stop("Download failed, please rerun the pipeline!")
}

dwl.status <- download.file(url = sprintf("%s/%s", cellFileURL_8.2, cellFileName_8.2), destfile = myfn2, quiet = TRUE)
if (dwl.status != 0) {
  stop("Download failed, please rerun the pipeline!")
}


# # require(gdata)
# # ver 8.0 cell data (July 2019)
# cell.info <- as.data.frame(read_excel(myfn, sheet = 1, .name_repair = make.names))
# cell.info <- cell.info[-nrow(cell.info), ]
# ## Last row is a total summation row

# cell.all <- read.csv(file.path(download_dir, "cell_annotation_all.csv"), na.strings = c("", " ", "NA"))


# cell.info$unique.cellid <- cell.all[match(cell.info[, "Sample.Name"], cell.all[, "GDSC1000.cellid"]), "unique.cellid"]
# save(cell.info, file = file.path(processed_dir, "cellInfo.RData"))



# # ver 8.2 cell data (Feb 2020)
# cell.info <- as.data.frame(read_excel(myfn2, sheet = 1, .name_repair = make.names))
# cell.info <- cell.info[-nrow(cell.info), ]
# cell.info$unique.cellid <- cell.all[match(cell.info[, "Sample.Name"], cell.all[, "GDSC1000.cellid"]), "unique.cellid"]
# save(cell.info, file = file.path(processed_dir, "cellInfo_8.2.RData"))

# cellcuration <- cell_all[,c("CGP.cellid", "GDSC.SNP.cellid", "CGP_EMTAB3610.cellid", "unique.cellid")]
# EMTAB3610_matches <- match(toupper(gsub(pattern=badchars, "", x=cell.info$Sample.Name)), toupper(gsub(pattern=badchars, "", x=cellcuration[,"CGP_EMTAB3610.cellid"])))
# SNP_matches <- match(toupper(gsub(pattern=badchars, "", x=cell.info$Sample.Name)), toupper(gsub(pattern=badchars, "", x=cellcuration[,"GDSC.SNP.cellid"])))
# CGP_matches <- match(toupper(gsub(pattern=badchars, "", x=cell.info$Sample.Name)), toupper(gsub(pattern=badchars, "", x=cellcuration[,"CGP.cellid"])))
# unique_matches <- match(toupper(gsub(pattern=badchars, "", x=cell.info$Sample.Name)), toupper(gsub(pattern=badchars, "", x=cellcuration[,"unique.cellid"])))

# matches <- cbind(EMTAB3610_matches, SNP_matches, CGP_matches, unique_matches)

# doubleMatchedCells <- apply(matches,1,function(x){
#   return(ifelse(!all(is.na(x)), !all(c(x[1]==x[2], x[1]==x[3], x[1]==x[4]), na.rm=T), FALSE))
# })

# matchesNames <- cbind(cell_all$unique.cellid[EMTAB3610_matches], cell_all$unique.cellid[SNP_matches], cell_all$unique.cellid[CGP_matches], cell_all$unique.cellid[unique_matches])

# ##### for whatever reason, the SNP_matches are correct while the others are wrong except for JUKRAT which the first column in correct
# matches <- apply(matches, 1, function(x){
#   if(!is.na(x[4])){
#     return(x[4])
#   } else {
#     return(ifelse(length(unique(na.omit(x[1:3]))), unique(na.omit(x[1:3])), NA))
#   }
# })
# cell.info$cellid <- cell.info$Sample.Name
# cell.info$cellid[!is.na(matches)] <- cell_all$unique.cellid[na.omit(matches)]
# cell.info$cellid[grep("KM-H2", cell.info$cellid)] <- cell.info$Sample.Name[grep("KM-H2", cell.info$cellid)]
# cell.info$cellid[grep("^T[-]T$", cell.info$Sample.Name)] <- "T.T"
# cell.info$cellid[grep("^TT$", cell.info$Sample.Name)] <- "TT"
# save(cell.info, file="/pfs/out/cellInfo.RData")
