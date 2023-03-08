library("readxl")
require(downloader)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")
processsed_dir <- paste0(args[1], "processed")

myfn <- file.path(download_dir, "gdsc_cellinfo.xlsx")
myfn2 <- file.path(download_dir, "gdsc_cellinfo_8.2.xlsx")

# require(gdata)
# ver 8.0 cell data (July 2019)
cell.info <- as.data.frame(read_excel(myfn, sheet = 1, .name_repair = make.names))
cell.info <- cell.info[-nrow(cell.info), ]
## Last row is a total summation row

cell.all <- read.csv(file.path(download_dir, "cell_annotation_all.csv"), na.strings = c("", " ", "NA"))


cell.info$unique.cellid <- cell.all[match(cell.info[, "Sample.Name"], cell.all[, "GDSC1000.cellid"]), "unique.cellid"]
save(cell.info, file = file.path(processed_dir, "cellInfo.RData"))

# ver 8.2 cell data (Feb 2020)
cell.info <- as.data.frame(read_excel(myfn2, sheet = 1, .name_repair = make.names))
cell.info <- cell.info[-nrow(cell.info), ]
cell.info$unique.cellid <- cell.all[match(cell.info[, "Sample.Name"], cell.all[, "GDSC1000.cellid"]), "unique.cellid"]
save(cell.info, file = file.path(processed_dir, "cellInfo_8.2.RData"))
