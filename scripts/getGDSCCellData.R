library(PharmacoGx)
library(RCurl)
options(stringsAsFactors = FALSE)

# Using the celline.gdsc.RData file obtained from a previous commit when the download link was valid.
# This is a temporary solution until the updated download link from COSMIC is found.
load(file.path("/pfs/downloadGDSCcells/celline.gdsc.RData"))
save(celline.gdsc, file = "/pfs/out/celline.gdsc.RData")

# getCosmic <- function(em, passw, directory="/pfs/out") {
#   if (missing(em)) { stop ("Email must be provided") }
#   if (missing(passw)) { stop ("Password must be provided") }
#
#
#   auth_key <- base64(paste(em, passw, sep=":"))
#
#   json_resp <- getURL("https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cell_lines/v87/CosmicCLP_MutantExport.tsv.gz",
#                       httpheader=c(Authorization=paste0("Basic ", auth_key)))
#   json_resp <- gsub(pattern="\\{\"url\":\"", rep="", json_resp)
#   json_resp <- gsub(pattern="\"\\}", rep="", json_resp)
#   resp <- RCurl::getBinaryURL(json_resp)
#   if (!file.exists(file.path(directory))) { dir.create(file.path(directory), showWarnings=FALSE, recursive=TRUE) }
#   myfl <- file(file.path(directory, "CosmicCLP_MutantExport.tsv.gz"), "wb")
#   writeBin(resp, myfl)
#   close(myfl)
#   unlink(file.path(directory, "cookies.txt"))
#   return (0)
# }
#
#
#
# generateGDSCCell.lines <- function(path.data="/pfs/out",  path.cell = file.path(path.data, "celline"), saveres=file.path("/pfs/out/saveres")){
#   # ftpdir <- "ftp://ftp.ebi.ac.uk//pub/databases/microarray/data/experiment/MTAB/E-MTAB-783/"
#   options(stringsAsFactors=FALSE)
#   if(!file.exists(path.cell)){
#     dir.create(path.cell)
#   }
#   if(!file.exists(saveres)){
#     dir.create(saveres)
#   }
#   ## download cell line annotations and COSMIC IDs
#   ## annotations from COSMIC cell line project
#   em <- "benjamin.haibe.kains@utoronto.ca"
#   passw <- "princessmargaret"
#   ## sshpass -p 'princessmargaret' scp -r 'benjamin.haibe.kains@utoronto.ca'@sftp-cancer.sanger.ac.uk:/files/grch38/cosmic/v73/CosmicCLP_MutantExport.tsv.gz dwl
# 
#   # system(sprintf("scp -r '%s':'%s'@sftp-cancer.sanger.ac.uk:/files/grch38/cosmic/v73/CosmicCLP_MutantExport.tsv.gz dwl", em, passw))
# 
#   ## download cell line annotations and COSMIC IDs
#   ## annotations from COSMIC cell line project
#   myfn <- file.path(path.cell, "cosmic_annotations.RData")
#   if(!file.exists(myfn)) {
#     message("Download COSMIC annotations for cell lines")
#     myfn2 <- file.path("cosmic_cell_line_collection.txt")
#     if(!file.exists(myfn2)) {
#       dir.create(file.path("dwl"), showWarnings=FALSE, recursive=TRUE)
#       dwl.status <- getCosmic(em=em, passw=passw, directory=file.path("dwl"))
#       if(dwl.status != 0) { stop("Download failed, please rerun the pipeline") }
#       ## untar
#       res <- R.utils::gunzip(filename=file.path("dwl", "CosmicCLP_MutantExport.tsv.gz"), overwrite=TRUE)
#       file.copy(from=file.path("dwl", "CosmicCLP_MutantExport.tsv"), to=myfn2)
#     }
#     message("Process COSMIC annotations")
#     cosmic.celline <- read.csv(file=myfn2, sep="\t")
#     cosmic.celline <- cosmic.celline[complete.cases(cosmic.celline[ , c("Sample.name", "Sample.source")]) & cosmic.celline[ , "Sample.source"] == "cell-line", , drop=FALSE]
#     cosmic.celline[cosmic.celline == "NS" | cosmic.celline == "" | cosmic.celline == " " | cosmic.celline == "  "] <- NA
#     ## merge the gene targets
#     dupln <- sort(unique(cosmic.celline[ , "Sample.name"][duplicated(cosmic.celline[ , "Sample.name"])]))
#     tt <- cosmic.celline
#     ## select unique cell lines
#     iix.rm <- NULL
#     for(i in 1:length(dupln)) {
#       duplix <- cosmic.celline[ ,"Sample.name"] == dupln[i]
#       iix <- sort((which(duplix)), decreasing=FALSE)[1]
#       iix.rm <- c(iix.rm, setdiff(which(duplix), iix))
#       ## get the most frequent tissue type
#       tissuet <- table(cosmic.celline[duplix, "Primary.site"])
#       if (length(tissuet) == 0) {
#         tt[iix, "Primary.site"] <- NA
#       } else {
#         tt[iix, "Primary.site"] <- names(sort(tissuet, decreasing=TRUE))[1]
#       }
#     }
#     tt <- tt[-iix.rm, , drop=FALSE]
#     tt <- tt[!is.na(tt[ , "Sample.name"]), , drop=FALSE]
#     rownames(tt) <- tt[ , "Sample.name"]
#     ## remove unnecessary annotations
#     tt <- tt[ , c("Sample.name", "ID_sample", "ID_tumour", "Primary.site", "Site.subtype.1", "Site.subtype.2", "Site.subtype.3", "Primary.histology", "Histology.subtype.1", "Histology.subtype.2", "Histology.subtype.3", "Sample.source", "Tumour.origin"), drop=FALSE]
#     cosmic.celline <- tt
#     save(list=c("cosmic.celline"), compress=TRUE, file=myfn)
#   } else { load(myfn) }
#   ## annotations from GDSC (Genomics of Drug Sensitivity in Cancer)
#   myfn <- file.path(saveres, "gdsc_annotations.RData")
#   if(!file.exists(myfn)) {
#     message("Download GDSC annotations for cell liness")
#     myfn2 <- file.path(path.cell, "gdsc_celline_collection.csv")
#     if(!file.exists(myfn2)) {
#       dir.create(file.path(path.cell, "dwl"), showWarnings=FALSE, recursive=TRUE)
#       dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_cell_lines_w5.csv", destfile=file.path(path.cell, "dwl", "gdsc_cell_lines_w5.csv"), quiet=TRUE)
#       if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
#       file.copy(from=file.path(path.cell, "dwl", "gdsc_cell_lines_w5.csv"), to=myfn2)
#     }
#     gdsc.celline <- read.csv(file=file.path(path.cell, "gdsc_celline_collection.csv"))
#     gdsc.celline[gdsc.celline == "" | gdsc.celline == " " | gdsc.celline == "  "] <- NA
#     gdsc.celline <- gdsc.celline[!is.na(gdsc.celline[ , "CELL_LINE_NAME"]), , drop=FALSE]
#     dupln <- unique(gdsc.celline[ , "CELL_LINE_NAME"][duplicated(gdsc.celline[ , "CELL_LINE_NAME"])])
#     gdsc.celline <- gdsc.celline[!duplicated(gdsc.celline[ , "CELL_LINE_NAME"]), , drop=FALSE]
#     rownames(gdsc.celline) <- gdsc.celline[ , "CELL_LINE_NAME"]
#     save(list=c("gdsc.celline"), compress=TRUE, file=myfn)
#   } else { load(myfn) }
# 
#   ## merge GDSC and COSMIC annotations through COSMIC_ID
#   message("Merge COSMIC and GDSC annotations for cell lines")
#   iix <- which(complete.cases(gdsc.celline[ , c("CELL_LINE_NAME", "COSMIC_ID")]) & !is.element(gdsc.celline[ , "COSMIC_ID"], cosmic.celline[ , "ID_sample"]) & !is.element(gdsc.celline[ , "CELL_LINE_NAME"], cosmic.celline[ , "Sample.name"]))
#   tt <- data.frame(matrix(NA, nrow=nrow(cosmic.celline) + length(iix), ncol=ncol(cosmic.celline), dimnames=list(c(rownames(cosmic.celline), rownames(gdsc.celline)[iix]), colnames(cosmic.celline))))
#   tt[rownames(cosmic.celline), ] <- cosmic.celline
#   tt[rownames(gdsc.celline)[iix], "Sample.name"] <- gdsc.celline[iix, "CELL_LINE_NAME"]
#   tt[rownames(gdsc.celline)[iix], "ID_sample"] <- gdsc.celline[iix, "COSMIC_ID"]
#   celline.gdsc <- tt
# 
#   return(celline.gdsc)
# 
# }
# 
# celline.gdsc <- generateGDSCCell.lines()
# 
# save(celline.gdsc, file="/pfs/out/celline.gdsc.RData")
