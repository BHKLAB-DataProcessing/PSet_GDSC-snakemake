library(data.table)
library(RaggedExperiment)
library(GenomicRanges)
# library(rtracklayer)

# source("~/Code/Github/PharmacoGx-private/R/matchToIDTable.R")

inputDir <- "/pfs/input"
outDir <- "/pfs/out"

cellData <- read.csv(file.path(inputDir, "sample_info.csv"))


allMuts <- fread(file.path(inputDir, "CCLE_mutations.csv"))



gdsc.muts <- allMuts[SangerRecalibWES_AC != ""]
gdsc.muts[, CGA_WES_AC := NULL]
gdsc.muts[, RNAseq_AC := NULL]
gdsc.muts[, HC_AC := NULL]
gdsc.muts[, RD_AC := NULL]
gdsc.muts[, WGS_AC := NULL]


gdscGRanges <- makeGRangesFromDataFrame(gdsc.muts, start.field = "Start_position", end.field = "End_position", keep.extra.columns = TRUE)

gdscGRangesSplit <- split(gdscGRanges, gdscGRanges$DepMap_ID)

gdscRag <- RaggedExperiment(gdscGRangesSplit)

colData(gdscRag)$DepMap_ID <- colnames(gdscRag)
colData(gdscRag)$stripped_cell_line_name <- cellData[match(colnames(gdscRag), cellData$DepMap_ID), "stripped_cell_line_name"]

saveRDS(gdscRag, file = file.path(outDir, "gdscMutationExtended.rds"))
