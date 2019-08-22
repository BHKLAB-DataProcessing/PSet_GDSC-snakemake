#################################################
## Re-Implementation of getGDSC for 2019 changes
## And pachyderm support.
## 
#################################################
library(PharmacoGxPrivate)
options(stringsAsFactors=FALSE)

getGDSC <- function(tmpdir=tempdir(),
										verbose=TRUE){

  if(!file.exists(tmpdir)) { dir.create(tmpdir, showWarnings=FALSE, recursive=TRUE) }
  badchars <- "[\xb5]|[]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  
  celline <- get(load("~/Data/getGDSCTest/gdscCellInfo/celline.gdsc.RData"))
  ## COSMIC has fixed this error in the cell line name, but the data still has the old name
  ## The error in the pSet will be fixed when we remap the names
  celline[celline[, "Sample.name"] == "NTERA-2_cl_D1", "Sample.name"] <- "NTERA-S-cl-D1"

  # TODO:: make sure this is actually needed
  cosmic <- celline 

  ## format column names
  message("Read drug sensitivity measurements")
  drugpheno <- read.csv(file.path("~/Data/getGDSCTest/gdscDrugInfo", "gdsc_drug_sensitivity.csv"), na.strings=c("NA", "", " "))
  coln2 <- unlist(drugpheno[1, ])
  coln2[coln2 == ""] <- NA  
  drugpheno <- drugpheno[!is.na(drugpheno[ , "Cell.Line"]), ,drop=FALSE]
  coln <- colnames(drugpheno)
  coln2[is.na(coln2)] <- coln[is.na(coln2)]
  coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
  myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
  coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
  colnames(drugpheno) <- coln2

  ## drug identifiers and names
  dn <- toupper(gsub(badchars, "", sapply(strsplit(coln, "_"), function(x) { return(x[[1]]) })))
  ## manual curation for drug names starting with a digit
  dn[!is.na(dn) & dn == "X17AAG"] <- "17AAG"
  dn[!is.na(dn) & dn == "X681640"] <- "681640"
  did <- sapply(strsplit(coln2, "_"), function(x) { if(x[[1]] == "drugid") { return(x[[2]]) } else { return(NA) } })
  drugnid <- cbind("drug.name"=dn, "drug.id"=did)[!is.na(did) & !duplicated(did), ]
  rownames(drugnid) <- paste("drugid", drugnid[ , "drug.id"], sep="_")

  # Check that all cell line names still exist in cosmic
	if(any(!is.element(drugpheno[ , "Cell.Line"], celline[ , "Sample.name"]))) { 
		stop("Some cell line names are not included in the COSMIC database") }
  #cellnames 
  celln <- drugpheno[ , "Cell.Line"]
  drugpheno <- data.frame("cellid"=celln, drugpheno)
  rownames(drugpheno) <- celln

  #assumption: Mutation columns range from AKT2 to MLL_AFF1
  rangeg <- which(colnames(drugpheno) == "AKT2"):which(colnames(drugpheno) == "MLL_AFF1")
  mutation <- as.matrix(drugpheno[ , rangeg, drop=FALSE])
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

  # Formating cosmic info
  celline <- data.frame("GDSC.cellid"=as.character(celline[ , "Sample.name"]), celline)
  celline[ , "GDSC.cellid"] <- as.character(celline[ , "GDSC.cellid"])
  ## add url based on COSMIC IDs
  uurl <- paste("http://cancer.sanger.ac.uk/cell_lines/sample/overview?id=", celline[ , "ID_sample"], sep="")
  uurl[is.na(celline[ , "ID_sample"])] <- NA
  celline <- data.frame("GDSC.cellid"=celline[ , "GDSC.cellid"], "link"=uurl, celline[ , !is.element(colnames(celline), "GDSC.cellid")])

  ## TODO: what are the lines below actually doing?
	## drugpheno
  cellnall <- sort(unique(c(as.character(drugpheno[ , "cellid"]))))
  # dd <- data.frame(matrix(NA, ncol=ncol(drugpheno), nrow=length(cellnall), dimnames=list(cellnall, colnames(drugpheno))))
  # newlev <- sapply(drugpheno, levels)
  # newlev$cellid <- cellnall
  # dd <- genefu::setcolclass.df(df=dd, colclass=sapply(drugpheno, class), factor.levels=newlev)
  # dd[rownames(drugpheno), colnames(drugpheno)] <- drugpheno
  # dd[ , "cellid"] <- cellnall
  # drugpheno <- dd
  
  ## mutation
  dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
  dd[rownames(mutation), colnames(mutation)] <- mutation
  rownames(dd) <- cellnall
  mutation <- dd  

   
	MutationEset <- ExpressionSet(t(mutation)) 
	geneMap <- read.csv("~/Data/getGDSCTest/downAnnotations/annot_ensembl_all_genes.csv")
	geneInfoM <- geneMap[na.omit(match(rownames(MutationEset),geneMap[ , "gene_name"]) ), c('gene_biotype','gene_name','EntrezGene.ID')] 
	rownames(geneInfoM) <- geneInfoM[ , "gene_name"]     
	geneInfoM <- geneInfoM[rownames(MutationEset),]      
	rownames(geneInfoM) <- rownames(MutationEset)
	fData(MutationEset) <- geneInfoM 
	tttt <- data.frame(row.names=colnames(MutationEset), colnames(MutationEset))
	colnames(tttt) <- 'cellid'
	pData(MutationEset) <- tttt
	annotation(MutationEset) <- "mutation"
	pData(MutationEset)[, "batchid"] <- NA


	## drug information
  message("Read drug information")
  druginfo <- read.csv(file.path("~/Data/getGDSCTest/gdscDrugInfo", "gdsc_drug_information.csv"), na.strings=c("NA", "", " "))
  druginfo <- data.frame("drug.name"=toupper(gsub(badchars, "", druginfo[ , "Name"])), druginfo)
  druginfo <- druginfo[,!colnames(druginfo) %in% c("Sample.Size", "Count")]
	## Fixing up CAMPTOTHECIN
  druginfo[druginfo[,"drug.name"]=="SN38",c("drug.name", "drug_id")] <- c("CAMPTOTHECIN", drugnid[drugnid[,"drug.name"]=="CAMPTOTHECIN","drug.id"])
  myx <- match(drugnid[ , "drug.id"], druginfo[ , "drug_id"])

  if (any(is.na(myx))) { warning ("Some drugs have missing annotations") }
  druginfo <- data.frame("drugid"=rownames(drugnid), drugnid, druginfo[myx, , drop=FALSE])
  rownames(druginfo) <- as.character(druginfo[ , "drugid"])
  ## complement drug infomration with the supplementary infomration from the Nature website
  druginfo.nature <- as.data.frame(readxl::read_xlsx(file.path("~/Data/getGDSCTest/gdscDrugInfo", "nature_supplementary_information.xlsx"), sheet=4,  .name_repair=make.names))
  druginfo.nature[druginfo.nature == "" | druginfo.nature == " "] <- NA
  rownames(druginfo.nature) <- paste("drugid", druginfo.nature[ , "Drug.ID"], sep="_")
  druginfo <- data.frame(druginfo, druginfo.nature[match(druginfo[,"drug_id"],druginfo.nature[ , "Drug.ID"]), c("Brand.name", "Site.of.screening", "Drug.type", "Drug.class.I", "Drug.class.II", "Target.family", "Effector.pathway.biological.process", "Clinical.trials", "Source")])


  ## drug concentration
  message("Read drug concentration")
  drugconc <- read.csv(file.path("~/Data/getGDSCTest/gdscDrugInfo", "gdsc_drug_concentration.csv"))
  drugconc[!is.na(drugconc) & (drugconc == "" | drugconc == " ")] <- NA
  drugconc <- data.frame("drug.name"=toupper(gsub(badchars, "", drugconc[ , "Compound.Name"])), drugconc)
  
  if(all(!is.element(drugconc[ , "drug.name"], drugnid[ , "drug.name"]))) { stop("Screening concentration for drugs without identifiers!") }
  
  myx <- match(drugconc[ , "drug.name"], drugnid[ , "drug.name"])
  drugid <- rownames(drugnid)[myx]
  drugconc <- data.frame("drugid"=drugid, drugconc)

  ## correct ambiguity for AZD6482: drugid_156 corresponds to the first occurence of AZD6482 while drugid_1066 corresponds to the second
	drugconc[drugconc$drug.name == "AZD6482" & drugconc$Max.Concentration.micromolar. == 5, "drugid"] <- "drugid_1066"

  rownames(drugconc) <- drugconc$drugid


  #Make tables to contain all drugs 
  drugpheno.drugs <- paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_")
  if(!all(unionList(c(drugpheno.drugs, rownames(drugconc))) %in% rownames(druginfo))){
  	stop("Some drugs in drugpheno or drugconc not in druginfo")
  }

 #  dix <- sort(unique(c(rownames(druginfo), rownames(drugconc), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
 #  druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))))
	# newlev <- sapply(druginfo, levels)
 #  newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
 #  druginfo2 <- genefu::setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
 #  druginfo2[match(rownames(druginfo), dix), colnames(druginfo)] <- druginfo
 #  druginfo2[ , "drugid"] <- newlev$drugid
 #  druginfo <- druginfo2
 #  ## update drugconc
 #  drugconc2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(drugconc), dimnames=list(dix, colnames(drugconc))))
 #  newlev <- sapply(drugconc, levels)
 #  newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
 #  drugconc2 <- genefu::setcolclass.df(df=drugconc2, colclass=sapply(drugconc, class), factor.levels=newlev)
 #  drugconc2[match(rownames(drugconc), dix), colnames(drugconc)] <- drugconc
 #  drugconc2[ , "drugid"] <- newlev$drugid
 #  drugconc <- drugconc2

  ## report concentrations per cell line and per drug
	drugconc2 <- data.frame(matrix(NA, nrow=nrow(drugconc) * length(cellnall), ncol=6, dimnames=list(paste(rep(rownames(drugconc), times=length(cellnall)), rep(cellnall, each=nrow(drugconc)), sep="_"), c("cellid", "drugid", "drug.name", "nbr.conc.tested", "min.Dose.uM", "max.Dose.uM"))))
  drugconc2[ , "drugid"] <- rep(rownames(drugconc), times=length(cellnall))
  drugconc2[ , "cellid"] <- rep(cellnall, each=nrow(drugconc))
  drugconc2[ , "drug.name"] <- rep(as.character(drugconc[ , "drug.name"]), times=length(cellnall))
  ## as mentioned in the supplementary information of Garnett et al., a single cell line is used on each plate and treated with 28 different drugs over a 9-pt, 256-fold concentration range
  drugconc2[ , "nbr.conc.tested"] <- 9
  drugconc2[ , "min.Dose.uM"] <- rep(drugconc[ , "Min.Concentration.micromolar."], times=length(cellnall))
  drugconc2[ , "max.Dose.uM"] <- rep(drugconc[ , "Max.Concentration.micromolar."], times=length(cellnall))
  drugconc <- drugconc2

  drugnall <- rownames(druginfo)

  ## update drugpheno
  ## IC50 in µM 
  iix <- grep("_IC_50$", colnames(drugpheno))
  iixn <- gsub("_IC_50$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.ic50 <- exp(dd) # XXX: What is going on here? do we have logIC50


	iix <- grep("_AUC$", colnames(drugpheno))
  iixn <- gsub("_AUC$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.auc <- 1 - dd

  # Slope not currently included in PSet
  if(FALSE){
	  iix <- grep("_BETA$", colnames(drugpheno))
	  iixn <- gsub("_BETA$", "", colnames(drugpheno)[iix])
	  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
	  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
	  dd[dd < 0] <- 0
	  dd[dd > 1] <- 1
	  drugpheno.slope <- dd
  }

  nms <- sapply(colnames(drugpheno.ic50), paste, rownames(drugpheno.ic50), sep="_")
  ic50 <- c(drugpheno.ic50)
  names(ic50) <- nms
  auc <- c(drugpheno.auc)
  names(auc) <- nms

  profiles <- cbind("ic50_published"=ic50, "auc_published"=auc)#, "slope_published"=slope)

  load("~/Data/getGDSCTest/gdscRawSensitivity/GDSC_sens_raw.RData")
  # con_tested <- con_tested
  
  if(!all(rownames(raw.sensitivity) %in% rownames(drugconc))){
  	stop("Raw sensitivity data for which we are missing drugconc (sensInfo)")
  }
  if(!all(rownames(raw.sensitivity) %in% rownames(profiles))){
  	stop("Raw sensitivity data for which we are missing published data")
  }
  ## NOTE: there are an equal number complete cases of profiles as rows of raw.sensitivity
  
  if(!sum(complete.cases(profiles)) == nrow(raw.sensitivity)){
  	stop("Unequal number of published and recomputed profiles")
  }

  drugconc <- drugconc[rownames(raw.sensitivity),]
  duration <- rep(x=72, length=nrow(drugconc))
  sensitivityInfo <- cbind(drugconc, "duration_h"=duration)
  profiles <- profiles[rownames(sensitivityInfo),]
  # recomputed <- .calculateFromRaw(raw.sensitivity,dose.range=c(log10(2), log10(1000)), cap=100)
  # myfn <- file.path(saveres, "GDSC_sens_recomputed.RData")
  # if(!file.exists(myfn)){
  # 	 recomputed <- .calculateFromRaw(raw.sensitivity, cap=100)
  #    save(recomputed, file=myfn)
  # } else {
	 #  load(myfn, verbose=TRUE)
  # }

  recomputed <- get(load("~/Data/getGDSCTest/gdscProfiles/gdscProfile.RData"))
	if(!setequal(rownames(recomputed), rownames(profiles)){
		stop("Recomputed values for which we dont have raw values???")
	}
  profiles <- cbind(profiles, recomputed[rownames(profiles),])    

	## TODO::: Continue from here!!

	## Load in molecular Data
  annot <- read.csv(file.path("~/Data/getGDSCTest/downAnnotations", "annot_ensembl_all_genes.csv"), stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)
  
	# U219 array
  myf <- "GDSC_U219_ENSG.RData"

  load(file.path("~/Data/getGDSCTest/gdscU219normalized", myf), verbose=TRUE)

  ### the new microarray data loaded to gdsc.u219.ensg
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
	  
  # Old U133a array
  myfn2 <- file.path("~/Data/getGDSCTest/gdscU133anormalized", "GDSC_U133a_ENSG.RData")
  load(myfn2)
  #eset <- just.rma(filenames=celfn, verbose=TRUE, cdfname="hgu133ahsensgcdf")
  # eset <- just.rma(filenames=celfn, verbose=TRUE, cdfname="hthgu133ahsensgcdf")
  eset <- cgp.u133a.ensg
  # pData(eset) <- as.data.frame(sampleinfo[match(gsub("[.]CEL[.]gz$", "", rownames(pData(eset))), rownames(sampleinfo)), , drop=FALSE])
  # colnames(exprs(eset)) <- rownames(pData(eset)) <- gsub("[.]CEL[.]gz$", "", colnames(exprs(eset)))
  controls <- rownames(exprs(eset))[grep("AFFX", rownames(exprs(eset)))]
  # fData(eset) <- fData(eset)[which(!rownames(fData(eset)) %in% controls), , drop=FALSE]
  # exprs(eset) <- exprs(eset)[which(!rownames(exprs(eset)) %in% controls), , drop=FALSE]
  eset <- eset[!rownames(eset) %in% controls,]
  ensemblIds <- sapply(strsplit(rownames(exprs(eset)), "_"), function (x) { return (x[[1]]) }) 
  fData(eset) <- data.frame("Probe"=rownames(exprs(eset)), 
                            "EnsemblGeneId"=ensemblIds,
                            "EntrezGeneId"=annot[ensemblIds, "EntrezGene.ID"],
                            "Symbol"=annot[ensemblIds, "gene_name"],
                            "GeneBioType"=annot[ensemblIds, "gene_biotype"],
                            "BEST"=TRUE)
  rownames(fData(eset)) <- rownames(exprs(eset))
  pData(eset)[,"batchid"] <- NA
  annotation(eset) <- "rna"

  ########################
  ## Fix up all cellids across all data types 
  ########################

  cell_all <- read.csv("~/Data/getGDSCTest/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
  rownames(cell_all) <- cell_all[, "unique.cellid"]

  ## First u133a

  mymatches <- matchToIDTable(cell_all, ids=pData(eset)[,"cellid"], column="CGP.cellid")["unique.cellid",]
  pData(eset)[,"cellid"] <- as.character(mymatches)

  # U219
  mymatches <- matchToIDTable(cell_all, ids=pData(gdsc.u219.ensg)[,"Factor.Value.cell_line."], column="CGP_EMTAB3610.cellid")["unique.cellid",]
  pData(gdsc.u219.ensg)[,"cellid"] <- as.character(mymatches)

  # Mutation
  mymatches <- matchToIDTable(cell_all, ids=pData(MutationEset)[,"cellid"], column="CGP.cellid")["unique.cellid",]
	pData(MutationEset)[,"cellid"] <- as.character(mymatches)

  # Sensitivity Info
  mymatches <- matchToIDTable(cell_all, ids=sensitivityInfo[ , "cellid"], column="CGP.cellid")["unique.cellid",]
	sensitivityInfo[ , "cellid"] <- as.character(mymatches)

	cellMatch1 <- matchToIDTable(cell_all, ids=celline[,"GDSC.cellid"], column="CGP.cellid")["unique.cellid",]
	cellMatch2 <- matchToIDTable(cell_all, ids=celline[,"GDSC.cellid"], column="CGP_EMTAB3610.cellid")["unique.cellid",]
	cellMatch3 <- matchToIDTable(cell_all, ids=celline[,"GDSC.cellid"], column="GDSC1000.cellid")["unique.cellid",]
	
	cellMatchAll <- do.call(rbind, list(cellMatch1, cellMatch2, cellMatch3))	

	mode(cellMatchAll) <- "character"

	cellMatchAll[cellMatchAll == "character(0)"] <- NA_character_

	cellMatchAll <- t(cellMatchAll)

	mybadx <- which(apply(cellMatchAll, 1, function(x) length(unique(na.omit(x)))>1))

  ########################
  ## merge data and put NA when the data are anot avaialble
  ########################
	message("Merge data")

	## cell lines and drugs
	# cellnall <- sort(unique(c(as.character(sampleinfo[ , "cellid"]), as.character(drugpheno[ , "cellid"]), rownames(genexprs), rownames(mutation))))
	sampleinfo.u113a <- pData(eset)
	sampleinfo.u219 <- pData(gdsc.u219.ensg)
	sampleinfo.mutation <- pData(MutationEset)

	cellnall <- sort(unique(c(as.character(sampleinfo.u113a[ , "cellid"]), 
										as.character(sensitivityInfo[ , "cellid"]), 
										as.character(sampleinfo.u219[ , "cellid"]),
										as.character(sampleinfo.mutation[, "cellid"]))))

	drugnall <- sort(unique(c(rownames(druginfo), as.character(drugconc[ , "drugid"]), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))

 





  ## Deduplicate AD6482

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
  druginfo <- collapseRows2(druginfo, which(druginfo$Name %in% "AZD6482"))



}

























