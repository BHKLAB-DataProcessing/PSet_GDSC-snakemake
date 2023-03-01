library(affy)
library(downloader)
library(PharmacoGx)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")
processsed_dir <- paste0(args[1], "processed")

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/processed"

load(file.path(processed_dir, "celline.gdsc.RData"))
celline <- celline.gdsc

celfile.timestamp <- read.csv(file.path(download_dir, "celfile_timestamp.csv"), row.names = 1)

file.paths <- file.path(download_dir, c(
  list.files(pattern = "^hthgu133ahsensg*", path = download_dir),
  list.files(pattern = "^pd.hthgu133a.hs.ensg*", path = download_dir)
))

print(file.paths)

install.packages(file.paths, repos = NULL, type = "source")

`celfileChip` <-
  function(filename) {
    h <- affyio::read.celfile.header(filename, info = "full")
    return(as.character(h$cdfName))
  }

`celfileDateHour` <-
  function(filename) {
    h <- affyio::read.celfile.header(filename, info = "full")
    # ddate <- grep("/", strsplit(h$DatHeader, " ")[[1]], value=TRUE)
    # ddate <- strsplit(ddate, split="/")[[1]]
    # CC <- ifelse(substr(ddate[3],1,1)=="9", "19", "20")
    if (length(h$ScanDate) > 0) {
      h$ScanDate <- gsub(pattern = "T", replacement = " ", x = h$ScanDate)
      ddate <- strsplit(h$ScanDate, " ")[[1]]
    } else {
      ddate <- rep(NA, 2)
    }
    names(ddate) <- c("day", "hour")
    return(ddate)
  }

dir.create(file.path(download_dir, 'gdsc_array'))
unzip(file.path(download_dir, 'gdsc_array.zip'), exdir = file.path(download_dir, 'gdsc_array'))
celfn <- list.celfiles(file.path(download_dir, 'gdsc_array'), full.names = TRUE)
celfns <- list.celfiles(file.path(download_dir, 'gdsc_array'), full.names = FALSE)
## experiments' names
names(celfn) <- names(celfns) <- gsub(".CEL.gz", "", celfns)
## chip type and date
chipt <- sapply(celfn, celfileChip)
chipd <- t(sapply(celfn, celfileDateHour))

message("Read sample information")
sampleinfo <- read.csv(file.path(download_dir, "gdsc_ge_sampleinfo.txt"), sep = "\t")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
## curate cell line names
sampleinfo[sampleinfo[, "Source.Name"] == "MZ2-MEL.", "Source.Name"] <- "MZ2-MEL"
sampleinfo[sampleinfo[, "Source.Name"] == "KMS-12-BM", "Source.Name"] <- "KMS12-BM"
iix <- which(!duplicated(sampleinfo[, "Source.Name"]) & !is.element(sampleinfo[, "Source.Name"], celline[, "Sample.name"]))
if (length(iix) > 0) {
  ## enrich the list of cell lines
  tt <- matrix(NA, nrow = length(iix), ncol = ncol(celline), dimnames = list(sampleinfo[iix, "Source.Name"], colnames(celline)))
  tt[, "Sample.name"] <- sampleinfo[iix, "Source.Name"]
  celline <- rbind(celline, tt)
}
fn <- gsub(patter = "[.]CEL", replacement = "", x = sampleinfo[, "Array.Data.File"])
if (any(!is.element(fn[!is.na(fn)], names(celfns)))) {
  stop("some CEL files are missing for the GDSC project")
}
rownames(sampleinfo) <- fn
sampleinfo <- sampleinfo[names(celfn), , drop = FALSE]
sampleinfo <- data.frame(
  "samplename" = names(celfns), 
  "filename" = celfns, 
  "chiptype" = chipt, 
  "hybridization.date" = chipd[, "day"], 
  "hybridization.hour" = chipd[, "hour"], 
  "file.day" = celfile.timestamp[, "file.day"], 
  "file.hour" = celfile.timestamp[, "file.hour"], 
  "batchid" = NA, 
  "cellid" = sampleinfo[, "Source.Name"], 
  sampleinfo
)




## phenodata
# load("/pfs/gdscU133a/celfile_timestamp.RData")

# rownames(celfile.timestamp) <- basename(rownames(celfile.timestamp))

# celfn <- list.files(pattern="*.CEL.gz", path="/pfs/gdscU133a/", full.names=TRUE)

cgp.u133a <- just.rma(filenames = celfn, cdfname = "hthgu133ahsensgcdf")
save(cgp.u133a, compress = TRUE, file = file.path(processed_dir, "GDSC_u133a_ENSG_RAW.RData"))
print(head(rownames(pData(cgp.u133a))))
colnames(cgp.u133a) <- gsub(colnames(cgp.u133a), pat = ".gz", rep = "", fixed = TRUE)
pData(cgp.u133a) <- data.frame(
  pData(cgp.u133a), 
  sampleinfo[match(colnames(exprs(cgp.u133a)), sampleinfo[, "Array.Data.File"]),, drop = FALSE], 
  celfile.timestamp[rownames(pData(cgp.u133a)), , drop = FALSE]
)
colnames(exprs(cgp.u133a)) <- rownames(pData(cgp.u133a)) <- colnames(exprs(cgp.u133a))
fData(cgp.u133a) <- data.frame(
  "PROBE" = rownames(exprs(cgp.u133a)), 
  "GENEID" = sapply(strsplit(rownames(exprs(cgp.u133a)), "_"), function(x) {
      return(x[[1]])
  }), 
"BEST" = TRUE
)
rownames(fData(cgp.u133a)) <- rownames(exprs(cgp.u133a))
cgp.u133a.ensg <- cgp.u133a


load(file.path(download_dir, "Ensembl.v99.annotation.RData"))
eset <- cgp.u133a.ensg
controls <- rownames(exprs(eset))[grep("AFFX", rownames(exprs(eset)))]
ensemblIds <- sapply(strsplit(rownames(exprs(eset)), "_at"), function(x) {
  return(x[[1]])
})
fData(eset) <- data.frame(
  "Probe" = rownames(exprs(eset)),
  "EnsemblGeneId" = ensemblIds,
  "Symbol" = features_gene[ensemblIds, "gene_name"],
  "GeneBioType" = features_gene[ensemblIds, "gene_biotype"]
)

rownames(fData(eset)) <- rownames(exprs(eset))
# rownames(fData(eset)) <- eset@featureData@data$EnsemblGeneId
# rownames(eset) <-  eset@featureData@data$EnsemblGeneId
# eset@featureData@data$EnsemblGeneId[1:68] <- NA
pData(eset)[, "batchid"] <- NA
annotation(eset) <- "rna"

cgp.u133a.ensg <- eset
save(cgp.u133a.ensg, compress = TRUE, file = file.path(processed_dir, "GDSC_U133a_ENSG.RData"))

unlink(file.path(download_dir, 'gdsc_array'), recursive = T)
