#!/usr/bin/env Rscript

library(PharmacoGx)
# library(PharmacoGxPrivate)
library(data.table)
library(Biobase)
library(reshape2)
library(dplyr)
options(stringsAsFactors=FALSE)

myDirPrefix <- "/pfs"

args = commandArgs(trailingOnly=TRUE)

version <- args[[1]]

matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
	sapply(ids, function(x) {
                          myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
                          if(length(myx) > 1){
                            stop("Something went wrong in curating ids, we have multiple matches")
                          }
			  if(length(myx) == 0){return(NA_character_)}
                          return(tbl[myx, returnColumn])
                        })
}


switch(version, v1 = {
	myOutFile <- "GDSC_v1.RData"
	myInPrefix <- "gdscv1"
  sensFolder <- "GDSC2019v1Normalize"
  profFolder <- "gdscProfilesV1"

	}, v2 = {
	myOutFile <- "GDSC_v2.RData"
	myInPrefix <- "gdscv2"
  sensFolder <- "GDSC2019v2Normalize"
  profFolder <- "gdscprofilesV2"

	})

message("Loading Sensitivity Data")

sens.info <- readRDS(file=file.path(myDirPrefix, sensFolder, paste0(myInPrefix, "_sens_info.rds")))
sens.raw <- readRDS(file=file.path(myDirPrefix, sensFolder, paste0(myInPrefix, "_sens_raw.rds")))
rownames(sens.info) <- sens.info$exp_id
rownames(sens.raw) <- sens.info$exp_id

# sens.recalc <- PharmacoGx:::.calculateFromRaw(sens.raw, nthread=10)

# saveRDS(sens.recalc, file=paste0(myInPrefix, "_sens_recalc.rds"))

# sens.recalc <- readRDS(paste0(myInPrefix, "_sens_recalc.rds"))
# sens.recalc$pars <- lapply(sens.recalc$pars, unlist)

# sens.pars <- do.call(rbind, sens.recalc$pars)

# sens.profiles <- cbind(data.frame("AAC" = sens.recalc$AUC, "IC50" = sens.recalc$IC50), sens.pars)

load(file.path(myDirPrefix, profFolder, "profiles.RData"))

sens.profiles <- res

sens.profiles <- sens.profiles[rownames(sens.info),]

message("Loading RNA Data")

load(file.path(myDirPrefix, "gdscU219normalized/GDSC_U219_ENSG.RData"))

cell.all <- read.csv(file.path(myDirPrefix, "downAnnotations/cell_annotation_all.csv"))


rna.cellid <- matchToIDTable(ids=phenoData(cgp.u219.ensg)$Characteristics.cell.line., tbl=cell.all, column = "CGP_EMTAB3610.cellid", returnColumn="unique.cellid")


message("Loading CNV Data")


#load(file.path(myDirPrefix, "gdscCNA/GDSC_eset.RData"))
cl.eset <- readRDS("/pfs/gdsc_cnv_new/GDSC_CN.gene.RDS")

cl.eset$GDSC.cellid <- as.character(cl.eset$`Sample Name`)
cnv.cellid <- as.character(matchToIDTable(ids=cl.eset$GDSC.cellid, tbl=cell.all, column="GDSC.SNP.cellid", returnColumn = "unique.cellid"))


#myx <- which(is.na(cl.eset$GDSC.cellid))

#toRep <- cl.eset$GDSC.cellid
#toRep[myx] <- sapply(strsplit(rownames(phenoData(cl.eset))[myx], split="_"), `[`, 1)

#cl.eset$GDSC.cellid <- toRep

# phenoData(cl.eset)[myx,"GDSC.cellid"] <- sapply(strsplit(rownames(phenoData(cl.eset))[myx], split="_"), `[`, 1)

#cnv.cellid <- matchToIDTable(ids=cl.eset$GDSC.cellid, tbl=cell.all, column="GDSC.SNP.cellid", returnColumn = "unique.cellid")

message("Loading Mutation/fusion Data")


mut.matrix <- read.csv(file.path(myDirPrefix, "gdscMutPanel/gdsc_mutation_w5.csv"))

mut.cellid <- matchToIDTable(ids=mut.matrix[,1], tbl=cell.all, column="CGP.cellid", returnColumn = "unique.cellid")


rangeg <- which(colnames(mut.matrix) == "AKT2"):which(colnames(mut.matrix) == "VHL")

mutation <- as.matrix(mut.matrix[ , rangeg, drop=FALSE])
mutation <- apply(X=mutation, MARGIN=c(1, 2), FUN=function(x) {
    x <- unlist(strsplit(x, split="::"))
    if(length(x) == 2) {
      if(!is.na(x[[1]]) && (x[[1]] == "na")) {
        x <- NA
      } else {
        x <- x[[1]]
      }
    } else { x <- NA }
    return(x)
  })

MutationEset <- ExpressionSet(t(mutation)) 

colnames(MutationEset) <- mut.cellid

geneMap <- read.csv(file.path(myDirPrefix, "downAnnotations/annot_ensembl_all_genes.csv"))
geneInfoM <- geneMap[na.omit(match(rownames(MutationEset),geneMap[ , "gene_name"]) ), c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")] 
rownames(geneInfoM) <- geneInfoM[ , "gene_name"]     
geneInfoM <- geneInfoM[rownames(MutationEset),]
colnames(geneInfoM) <- c("EnsemblGeneId", "EntrezGeneId", "Symbol", "GeneBioType")
rownames(geneInfoM) <- rownames(MutationEset)
fData(MutationEset) <- geneInfoM 
tttt <- data.frame(row.names=colnames(MutationEset), colnames(MutationEset))
colnames(tttt) <- 'cellid'
pData(MutationEset) <- tttt
pData(MutationEset)[, "batchid"] <- NA
annotation(MutationEset) <- "mutation"



rangefus <- which(colnames(mut.matrix) == "BCR_ABL"):which(colnames(mut.matrix) == "MLL_AFF1")

fusion <- as.matrix(mut.matrix[ , rangefus, drop=FALSE])
fusion <- apply(X=fusion, MARGIN=c(1, 2), FUN=function(x) {
    # x <- unlist(strsplit(x, split="::"))
    if(x == "") {
        x <- NA
      } else if (x == "0"){
        x <- "wt"
      }
    return(x)
  })



FusionEset <- ExpressionSet(t(fusion)) 

colnames(FusionEset) <- mut.cellid
tttt <- data.frame(row.names=colnames(FusionEset), colnames(FusionEset))
colnames(tttt) <- 'cellid'
pData(FusionEset) <- tttt
annotation(FusionEset) <- "fusion"
pData(FusionEset)[, "batchid"] <- NA


message("Loading Cell and Drug Info")


load(file.path(myDirPrefix, "gdsc1000CellInfo/cellInfo.RData"))
load(file.path(myDirPrefix, "gdscDrugInfo/drugInfo.RData"))


rownames(cell.info) <- cell.info$unique.cellid



message("Summarizing Processed RNAseq")

summarizeRnaSeq <- function (dir, 
                                 tool=c("kallisto", "stringtie", "cufflinks", "rsem", "salmon"), 
                                 features_annotation,
                                 samples_annotation) {
      library(Biobase)
      library(readr)
      library(tximport)
      
      load(features_annotation)
      tx2gene <- as.data.frame(cbind("transcript"=toil.transcripts$transcript_id, "gene"=toil.transcripts$gene_id))
      
      files <- list.files(dir, recursive = TRUE, full.names = T)
      resFiles <- grep("abundance.h5", files)
      resFiles <- files[resFiles]
      length(resFiles)
      names(resFiles) <- basename(dirname(resFiles))
      
      txi <- tximport(resFiles, type="kallisto", tx2gene=tx2gene)
      head(txi$counts[,1:5])
      dim(txi$counts)
      
      xx <- txi$abundance
      gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
      fData(gene.exp) <- toil.genes[featureNames(gene.exp),]
      pData(gene.exp) <- samples_annotation[sampleNames(gene.exp),]
      annotation(gene.exp) <- "rnaseq"
      
      xx <- txi$counts
      gene.count <- Biobase::ExpressionSet(log2(xx + 1))
      fData(gene.count) <- toil.genes[featureNames(gene.count),]
      pData(gene.count) <- samples_annotation[sampleNames(gene.count),]
      annotation(gene.count) <- "rnaseq"
      
      txii <- tximport(resFiles, type="kallisto", txOut=T)
      
      xx <- txii$abundance
      transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
      fData(transcript.exp) <- toil.transcripts[featureNames(transcript.exp),]
      pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp),]
      annotation(transcript.exp) <- "isoforms"
      
      xx <- txii$counts
      transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
      fData(transcript.count) <- toil.transcripts[featureNames(transcript.count),]
      pData(transcript.count) <- samples_annotation[sampleNames(transcript.count),]
      annotation(transcript.count) <- "isoforms"
      
      return(list("rnaseq"=gene.exp, 
                  "rnaseq.counts"=gene.count, 
                  "isoforms"=transcript.exp, 
                  "isoforms.counts"=transcript.count))
    }


rnaseq.sampleinfo <- read.csv("/pfs/downAnnotations/E-MTAB-3983.sdrf.txt", sep="\t")
rownames(rnaseq.sampleinfo) <- rnaseq.sampleinfo$Comment.EGA_RUN.
rnaseq.sampleinfo$cellid <- matchToIDTable(ids=rnaseq.sampleinfo$Source.Name, tbl=cell.all, column = "GDSC_rnaseq.cellid", returnColumn = "unique.cellid")
rnaseq.sampleinfo <- rnaseq.sampleinfo[,c("cellid","Characteristics.organism.part.","Characteristics.disease.","Characteristics.sex.","Scan.Name","Comment.EGA_RUN.")]
   
rnaseq <- summarizeRnaSeq(dir="/pfs/downloadgdscrnaseq/KallistoGDSC_hg38/KallistoGDSC_hg38", 
                                tool="kallisto", 
                                features_annotation="/pfs/downAnnotations/Gencode.v23.annotation.RData",
                                samples_annotation=rnaseq.sampleinfo)



message("Compile All GDSC Mutation Data")

######### Loading ALL exome data #########
mutation_raw <- read.csv("/pfs/gdscmutation_all/mutations_latest.csv", na.strings=c("", " ", "NA"))
mutation_raw <- mutation_raw[,c("gene_symbol","protein_mutation","model_name","cancer_driver")]
mutation_raw <- mutation_raw[which(mutation_raw$cancer_driver=="True"),]
cells_matched <- matchToIDTable(ids = mutation_raw[,3], tbl = cell.all, column = "GDSC.SNP.cellid", returnColumn = "unique.cellid")
mutation_raw[,3] <- as.character(cells_matched)
mutation_raw$cancer_driver <- NULL

#concatenate cases where one cell line maps to the same gene twice ("///")
xx <- mutation_raw %>% group_by(gene_symbol, model_name) %>% 
  mutate(protein_mutation = paste(protein_mutation, collapse="///"))

xx_df <- data.frame(xx)

#removes duplicated concatenation
xx_df_2 <- distinct(xx_df,gene_symbol, model_name,.keep_all = TRUE)

#flatten data frame to gene name x cell line matrix
matrix_final <- reshape2::acast(xx_df_2, gene_symbol ~ model_name, value.var = "protein_mutation")
matrix_final[which(is.na(matrix_final))] <- "wt"


geneMap <- read.csv("/pfs/downAnnotations/annot_ensembl_all_genes.csv")
geneInfoM <- geneMap[na.omit(match(rownames(matrix_final),geneMap[ , "gene_name"])), c('gene_biotype','gene_name','EntrezGene.ID')] 
rownames(geneInfoM) <- geneInfoM[ , "gene_name"] 
missing_genes <- rownames(matrix_final)[which(!rownames(matrix_final) %in% geneMap$gene_name)]
all_genes <- c(rownames(geneInfoM), missing_genes)
geneInfoM[nrow(geneInfoM)+ length(missing_genes),] <- NA
rownames(geneInfoM) <- all_genes

geneInfoM <- geneInfoM[rownames(matrix_final),] 

MutationAll <- Biobase::ExpressionSet(matrix_final)
tttt <- data.frame(row.names=colnames(MutationAll), colnames(MutationAll))
colnames(tttt) <- 'cellid'
pData(MutationAll) <- tttt
fData(MutationAll) <- geneInfoM 
annotation(MutationAll) <- "mutation"
pData(MutationAll)[, "batchid"] <- NA



cellnall <- unionList(rownames(cell.info), 
					  cnv.cellid, 
					  rna.cellid, 
					  mut.cellid,
		     			  rnaseq$rnaseq$cellid,
		     			  MutationAll$cellid)
newcells <- setdiff(cellnall, rownames(cell.info))
newRows <- matrix(NA_character_, nrow=length(newcells), ncol=ncol(cell.info))
# newRows <- cell.info[newcells,]

rownames(newRows) <- newcells
colnames(newRows) <- colnames(cell.info)
newRows[,"unique.cellid"] <- newcells

cell.info <- rbind(cell.info, newRows)

collapseRows2 <- function(x, rows){
    xNew <- lapply(x[rows, ], function(x) {
      xx <- na.omit(x)
      if (length(xx) == 0) {
        xx <- NA
      }
      if (length(unique(xx)) > 1) {
        xx <- paste(xx, collapse="///")
      } else {xx <- xx[1]}
      return(as.vector(xx))
      })
    xNew <- as.data.frame(xNew, as.is = TRUE)
    x[rows[1], ] <- xNew
    x <- x[-rows[-1], ]
    return(x)
  }

message("Deduplicating Drugs")


drugDupsIDs <- unique(drug.info$unique.drugid[duplicated(drug.info$unique.drugid)])

for(dupID in drugDupsIDs){
	myx <- which(drug.info$unique.drugid == dupID)
	drug.info <- collapseRows2(drug.info, myx)
}
rownames(drug.info) <- drug.info$unique.drugid

message("Making Curation Tables")

curationCell <- data.frame(unique.cellid = rownames(cell.info),
						   GDSC2019.cellid = cell.info$Sample.Name,
						   CGP.cellid = NA_character_,
						   GDSC.SNP.cellid = NA_character_,
						   CGP_EMTAB3610.cellid = NA_character_,
			  			   GDSC_rnaseq.cellid = NA_character_,
			  		           GDSC1000.cellid = NA_character_)
rownames(curationCell) <- curationCell$unique.cellid

myx <- match(rownames(curationCell),cell.all$unique.cellid)

curationCell$CGP.cellid <- cell.all[myx, "CGP.cellid"]
curationCell$GDSC.SNP.cellid <- cell.all[myx, "GDSC.SNP.cellid"]
curationCell$CGP_EMTAB3610.cellid <- cell.all[myx, "CGP_EMTAB3610.cellid"]
curationCell$GDSC_rnaseq.cellid <- cell.all[myx, "GDSC_rnaseq.cellid"]
curationCell$GDSC1000.cellid <- cell.all[myx, "GDSC1000.cellid"]

cell.info$tissueid <- cell.all[myx, "unique.tissueid"]

curationTissue <- data.frame("unique.tissueid" = cell.info$tissueid, "GDSC2019.tissueid" = cell.info$GDSC..Tissue.descriptor.1)
rownames(curationTissue) <- rownames(cell.info)

curationDrug <- data.frame(unique.drugid = drug.info$unique.drugid,
						   "GDSC2019.drugid" = drug.info$DRUG_NAME)
rownames(curationDrug) <- rownames(drug.info)

annot <- geneMap
rownames(annot) <- annot$gene_id
gdsc.u219.ensg <- cgp.u219.ensg
annotation(gdsc.u219.ensg) <- "rna"
ensemblIds <- sapply(strsplit(rownames(exprs(gdsc.u219.ensg)), "_"), function (x) { return (x[[1]]) }) 
fData(gdsc.u219.ensg) <- data.frame("Probe"=rownames(exprs(gdsc.u219.ensg)), 
                          "EnsemblGeneId"=ensemblIds,
                          "EntrezGeneId"=annot[ensemblIds, "EntrezGene.ID"],
                          "Symbol"=annot[ensemblIds, "gene_name"],
                          "GeneBioType"=annot[ensemblIds, "gene_biotype"],
                          "BEST"=TRUE)
rownames(fData(gdsc.u219.ensg)) <- rownames(exprs(gdsc.u219.ensg))
pData(gdsc.u219.ensg)[,"batchid"] <- NA
pData(gdsc.u219.ensg)[,"cellid"] <- rna.cellid

#Compile CNV data

tt <- rownames(pData(cl.eset))
pData(cl.eset) <- as.data.frame(apply(pData(cl.eset), MARGIN=2, as.character), stringsAsFactors=FALSE)
rownames(pData(cl.eset)) <- tt
pData(cl.eset)[,"batchid"] <- NA
pData(cl.eset)[,"cellid"] <- cnv.cellid
tt <- annot[match(rownames(fData(cl.eset)), annot$gene_name), c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")]
rownames(tt) <- rownames(fData(cl.eset))
colnames(tt) <- c("EnsemblGeneId", "EntrezGeneId", "Symbol", "GeneBioType")
fData(cl.eset) <- tt
annotation(cl.eset) <- "cnv"


cellsPresent <- sort(unionList(sens.info$cellid, 
					  pData(gdsc.u219.ensg)$cellid, 
					  pData(MutationEset)$cellid,
					  pData(FusionEset)$cellid,
					  pData(cl.eset)$cellid,
		    			  rnaseq$rnaseq$cellid,
			      		  MutationAll$cellid))
cell.info <- cell.info[cellsPresent,]

drugsPresent <- sort(unique(sens.info$drugid))

drug.info <- drug.info[drugsPresent,]


drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
drug_all <- drug_all[which(!is.na(drug_all[ , "GDSC2019.drugid"])),]
drug_all <- drug_all[ , c("unique.drugid", "GDSC2019.drugid","smiles","inchikey","cid","FDA")]
rownames(drug_all) <- drug_all[ , "unique.drugid"]

drug_all <- drug_all[rownames(drug.info),]
drug.info[,c("smiles","inchikey","cid","FDA")] <- drug_all[,c("smiles","inchikey","cid","FDA")]


message("Making PSet")



GDSC <- PharmacoSet(molecularProfiles=list("rna"=gdsc.u219.ensg, "mutation"=MutationEset, "mutation_exome"=MutationAll ,"fusion"=FusionEset, "cnv"=cl.eset, "rnaseq"=rnaseq$rnaseq, "rnaseq.counts"=rnaseq$rnaseq.counts, "isoforms"=rnaseq$isoforms, "isoforms.counts"=rnaseq$isoforms.counts),
                      name=paste("GDSC", version, sep="_"), 
                      cell=cell.info, 
                      drug=drug.info, 
                      sensitivityInfo=sens.info, 
                      sensitivityRaw=sens.raw, 
                      sensitivityProfiles=sens.profiles, 
                      sensitivityN=NULL, 
                      curationCell=curationCell, 
                      curationDrug=curationDrug, 
                      curationTissue=curationTissue, 
                      datasetType="sensitivity")

message("Saving")


save(GDSC, file=paste0("/pfs/out/GDSC", version, ".RData"), version=2)
saveRDS(GDSC, file=paste0("/pfs/out/GDSC", version, ".rds"), version=2)




