library(downloader)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")

download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGDSC/download"

basePath <- "http://mbni.org/customcdf/20.0.0/ensg.download"
download(file.path(basePath, "hgu219hsensgcdf_20.0.0.tar.gz"), destfile = file.path(download_dir, "hgu219hsensgcdf_20.0.0.tar.gz"))
download(file.path(basePath, "hgu219hsensgprobe_20.0.0.tar.gz"), destfile = file.path(download_dir, "hgu219hsensgprobe_20.0.0.tar.gz"))
download(file.path(basePath, "pd.hgu219.hs.ensg_20.0.0.tar.gz"), destfile = file.path(download_dir, "pd.hgu219.hs.ensg_20.0.0.tar.gz"))
download(file.path(basePath, "hthgu133ahsensgcdf_20.0.0.tar.gz"), destfile = file.path(download_dir, "hthgu133ahsensgcdf_20.0.0.tar.gz"))
download(file.path(basePath, "hthgu133ahsensgprobe_20.0.0.tar.gz"), destfile = file.path(download_dir, "hthgu133ahsensgprobe_20.0.0.tar.gz"))
download(file.path(basePath, "pd.hthgu133a.hs.ensg_20.0.0.tar.gz"), destfile = file.path(download_dir, "pd.hthgu133a.hs.ensg_20.0.0.tar.gz"))
download(file.path(basePath, "hgu133ahsensgcdf_20.0.0.tar.gz"), destfile = file.path(download_dir, "hgu133ahsensgcdf_20.0.0.tar.gz"))
download(file.path(basePath, "hgu133ahsensgprobe_20.0.0.tar.gz"), destfile = file.path(download_dir, "hgu133ahsensgprobe_20.0.0.tar.gz"))
