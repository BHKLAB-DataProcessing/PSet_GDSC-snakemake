#################################################
## Get fRMA normalized GDSC data from InSilicoDB
##
## 
#################################################
library(PharmacoGxPrivate)
getGDSC <- 
function (#gene=TRUE, 
			tmpdir=tempdir(),
  		verbose=TRUE, 
			nthread=1) {
  
  ## TODO:: DECIDE add another parameter to ask user if Pset should be created with standardized ids for the cell lines, drugs, tissues or the raw ones
  options(stringsAsFactors=FALSE)
  ## maximum number of CEL files to normalize at once
  
  
  ## create directories for temporary files
  if(!file.exists(tmpdir)) { dir.create(tmpdir, showWarnings=FALSE, recursive=TRUE) }
  
  badchars <- "[\xb5]|[]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  
  
  # ## create directories
  # if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
  # if(!file.exists(path.ge)) { dir.create(path.ge, showWarnings=FALSE, recursive=TRUE) }
  # if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
  # if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }
  # if(!file.exists(path.mut)) { dir.create(path.mut, showWarnings=FALSE, recursive=TRUE) }
  # if(!file.exists(path.fus)) { dir.create(path.fus, showWarnings=FALSE, recursive=TRUE) }
  
  # if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE, recursive=TRUE) }
  
  
  # celfile.timestamp <- downloadGDSCCel(path.data=path.data, path.ge=path.ge)
  
  # downloadGDSCsamples(path.data=path.data, path.ge=path.ge)
  
  # downloadGDSCdrugs(path.data=path.data, path.drug=path.drug)
  
  # celline <- generateGDSCCell.lines(path.data=path.data, path.cell=path.cell, saveres=saveres)
  celline <- get(load("/pfs/gdscCellInfo/celline.gdsc.RData"))
  celline[celline[, "Sample.name"] == "NTERA-2_cl_D1", "Sample.name"] <- "NTERA-S-cl-D1"

  cosmic <- celline 
  
  ## profiles for the drugs
  message("Read drug sensitivity measurements")
  # myfn2 <- file.path("/pfs/gdscDrugInfo", "gdsc_drug_sensitivity.RData")
  # if(!file.exists(myfn2)) {
  drugpheno <- read.csv(file.path("/pfs/gdscDrugInfo", "gdsc_drug_sensitivity.csv"))
  drugpheno[drugpheno == "" | drugpheno == " "] <- NA
    # save(list="drugpheno", compress=TRUE, file=myfn2)
  # } else { load(myfn2) }
  ## format column names
  coln2 <- unlist(drugpheno[1, ,drop=TRUE])
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
  ## manual curation for drug names starting with a figure
  dn[!is.na(dn) & dn == "X17AAG"] <- "17AAG"
  dn[!is.na(dn) & dn == "X681640"] <- "681640"
  did <- sapply(strsplit(coln2, "_"), function(x) { if(x[[1]] == "drugid") { return(x[[2]]) } else { return(NA) } })
  drugnid <- cbind("drug.name"=dn, "drug.id"=did)[!is.na(did) & !duplicated(did), ]
  rownames(drugnid) <- paste("drugid", drugnid[ , "drug.id"], sep="_")
  
  # drugpheno[drugpheno[ , "Cell.Line"] == "NTERA-S-cl-D1" , "Cell.Line"] <- "NTERA-2_cl_D1"
  if(any(!is.element(drugpheno[ , "Cell.Line"], celline[ , "Sample.name"]))) { stop("Some cell line names are not included in the COSMIC database") }
  celln <- drugpheno[ , "Cell.Line"]
  drugpheno <- data.frame("cellid"=celln, drugpheno)
  rownames(drugpheno) <- celln
  
  ## protein coding variants
  ## Genetic mutation data for cancer genes. Includes MSI status (1=unstable and 0=stable) and gene-fusions. A binary code 'x::y' description is used for each gene where 'x' identifies a coding variant and 'y' indicates copy number information from SNP6.0 data. For gene fusions, cell lines are identified as fusion not-detected (0) or the identified fusion is given. The following abbreviations are used: not analysed (na), not detected or wild-type (wt), no copy number information (nci).
  ## we assume that AKT2 and WT1 are the first and last genes in the file
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
  # write.csv(mutation, file=file.path(path.mut, "gdsc_mutation.csv"))
  
  
  ### info about each experiment
  # message("Read sample information")
  # sampleinfo <- read.csv(file.path("/pfs/gdscU133a", "gdsc_ge_sampleinfo.txt"), sep="\t")
  # sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
  # ## curate cell line names
  # sampleinfo[sampleinfo[ , "Source.Name"] == "MZ2-MEL.", "Source.Name"] <- "MZ2-MEL"
  # sampleinfo[sampleinfo[ , "Source.Name"] == "KMS-12-BM", "Source.Name"] <- "KMS12-BM"
  # iix <- which(!duplicated(sampleinfo[ , "Source.Name"]) & !is.element(sampleinfo[ , "Source.Name"], celline[ , "Sample.name"]))
  # if(length(iix) > 0) {
  #   ## enrich the list of cell lines
  #   tt <- matrix(NA, nrow=length(iix), ncol=ncol(celline), dimnames=list(sampleinfo[iix, "Source.Name"], colnames(celline)))
  #   tt[ , "Sample.name"] <- sampleinfo[iix, "Source.Name"]
  #   celline <- rbind(celline, tt)
  # }
  # fn <- gsub(patter="[.]CEL", replacement="", x=sampleinfo[ , "Array.Data.File"])
  # if(any(!is.element(fn[!is.na(fn)], names(celfns)))) { stop("some CEL files are missing for the GDSC project") }
  # rownames(sampleinfo) <- fn
  # sampleinfo <- sampleinfo[names(celfn), , drop=FALSE]
  # sampleinfo <- data.frame("samplename"=names(celfns), "filename"=celfns, "chiptype"=chipt, "hybridization.date"=chipd[ , "day"], "hybridization.hour"=chipd[ , "hour"], "file.day"=celfile.timestamp[ , "file.day"], "file.hour"=celfile.timestamp[ , "file.hour"], "batchid"=NA, "cellid"=sampleinfo[ , "Source.Name"], sampleinfo)
  # sampleinfo2 <- sampleinfo
  ## remove duplcated cell line hybridization
  
  
  ## TODO:: DECIDE what to do with the ALL option and rownames
  
  ## update of gdsc cell line collection
  celline <- data.frame("GDSC.cellid"=as.character(celline[ , "Sample.name"]), celline)
  celline[ , "GDSC.cellid"] <- as.character(celline[ , "GDSC.cellid"])
  ## add url based on COSMIC IDs
  uurl <- paste("http://cancer.sanger.ac.uk/cell_lines/sample/overview?id=", celline[ , "ID_sample"], sep="")
  uurl[is.na(celline[ , "ID_sample"])] <- NA
  celline <- data.frame("GDSC.cellid"=celline[ , "GDSC.cellid"], "link"=uurl, celline[ , !is.element(colnames(celline), "GDSC.cellid")])
  
  ## drugpheno
  #cellnall <- sort(unique(c(as.character(sampleinfo[ , "cellid"]), as.character(drugpheno[ , "cellid"]))))
  cellnall <- sort(unique(c(as.character(drugpheno[ , "cellid"]))))
  dd <- data.frame(matrix(NA, ncol=ncol(drugpheno), nrow=length(cellnall), dimnames=list(cellnall, colnames(drugpheno))))
  newlev <- sapply(drugpheno, levels)
  newlev$cellid <- cellnall
  dd <- genefu::setcolclass.df(df=dd, colclass=sapply(drugpheno, class), factor.levels=newlev)
  dd[rownames(drugpheno), colnames(drugpheno)] <- drugpheno
  dd[ , "cellid"] <- cellnall
  drugpheno <- dd
  
  ## mutation
  dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
  dd[rownames(mutation), colnames(mutation)] <- mutation
  rownames(dd) <- cellnall
  mutation <- dd
  
  
  
  
  # internal check only
  
  # ## reproducibility between different screening sites
  # ## camptothecin was screened at MGH (drug id 195) and WTSI (drug id 1003)
  # ## data only available in the supplementary infomration of the Nature website
  # myfn2 <- file.path(saveres, "nature_supplinfo_drugpheno_gdsc.RData")
  # if(!file.exists(myfn2)) {
  #   drugpheno.nature <- gdata::read.xls(xls=file.path(path.drug, "nature_supplementary_information.xls"), sheet=2)
  #   drugpheno.nature[drugpheno.nature == "" | drugpheno.nature == " "] <- NA
  #   save(list="drugpheno.nature", compress=TRUE, file=myfn2)
  # } else { load(myfn2) }
  # ## format column names
  # coln2 <- gsub(" ", "", sapply(drugpheno.nature[1,], as.character))
  # coln2[coln2 == ""] <- NA
  # drugpheno.nature <- drugpheno.nature[-1, ,drop=FALSE]
  # coln <- colnames(drugpheno.nature)
  # coln2[is.na(coln2)] <- coln[is.na(coln2)]
  # coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
  # myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
  # coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
  # colnames(drugpheno.nature) <- coln2
  # rownames(drugpheno.nature) <- as.character(drugpheno.nature[ , "Cell.Line"])
  # myx <- sapply(strsplit(colnames(drugpheno.nature), "_"), function(x) { return(all(x[c(length(x)-1, length(x))] == c("IC", "50"))) })
  # ic50 <- drugpheno.nature[ , myx, drop=FALSE]
  # nn <- dimnames(ic50)
  # nn[[2]] <- gsub("_IC_50", "", nn[[2]])
  # ic50 <- apply(ic50, 2, as.numeric)
  # dimnames(ic50) <- nn
  # ic50 <- exp(ic50) / 10^6
  # ## camptothecin
  # pdf(file.path(saveres, "gdsc_camptothecin_mgh_wtsi_paper.pdf"))
  # yy <- -log10(ic50[ , "drugid_195", drop=FALSE])
  # xx <- -log10(ic50[ , "drugid_1003", drop=FALSE])
  # ccix <- complete.cases(xx, yy)
  # nnn <- sum(ccix)
  # cc <- stats::cor.test(x=xx, y=yy, method="spearman", use="complete.obs", alternative="greater")
  # cci <- spearmanCI(x=cc$estimate, n=sum(ccix))
  # par(mar=c(4, 4, 3, 1) + 0.1)
  # llim <- round(range(c(xx, yy), na.rm=TRUE) * 10) / 10
  # myScatterPlot(x=xx, y=yy, xlab="-log10 IC50 (WTSI)", ylab="-log10 IC50 (MGH)", main="CAMPTOTHECIN", pch=16, method="transparent", transparency=0.75)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2)
  # dev.off()
  #
  
  
  ## drug information
  message("Read drug information")
  druginfo <- read.csv(file.path("/pfs/gdscDrugInfo", "gdsc_drug_information.csv"))
  druginfo[!is.na(druginfo) & (druginfo == " " | druginfo == " ")] <- NA
  druginfo <- data.frame("drug.name"=toupper(gsub(badchars, "", druginfo[ , "Name"])), druginfo)
  myx <- match(drugnid[ , "drug.id"], druginfo[ , "drug_id"])
  myx[drugnid[,"drug.name"]=="CAMPTOTHECIN"] <- grep("SN-38", druginfo[,"Name"]) ## Fixing up CAMPTOTHECIN
  if (any(is.na(myx))) { warning ("Some drugs have missing annotations") }
  ## correct ambiguity for AZD6482: drugid_156 corresponds to the first occurence of AZD6482 while drugid_1066 corresponds to the second
  ## table(!is.na(drugpheno[ , "drugid_156_AUC"]))
  ## table(!is.na(drugpheno[ , "drugid_1066_AUC"]))
  # myx[drugnid[ , "drug.name"] == "AZD6482"][2] <- which(druginfo[ , "drug.name"] == "AZD6482")[2]
  druginfo <- data.frame("drugid"=rownames(drugnid), drugnid, druginfo[myx, , drop=FALSE])
  rownames(druginfo) <- as.character(druginfo[ , "drugid"])
  ## complement drug infomration with the supplementary infomration from the Nature website
  # myfn2 <- file.path(saveres, "nature_supplinfo_druginfo_gdsc.RData")
  # if(!file.exists(myfn2)) {
    druginfo.nature <- as.data.frame(readxl::read_xlsx(file.path("/pfs/gdscDrugInfo", "nature_supplementary_information.xls"), sheet=4,  .name_repair=make.names))
    druginfo.nature[druginfo.nature == "" | druginfo.nature == " "] <- NA
    # save(list="druginfo.nature", compress=TRUE, file=myfn2)
  # } else { load(myfn2) }
  rownames(druginfo.nature) <- paste("drugid", druginfo.nature[ , "Drug.ID"], sep="_")
  druginfo <- data.frame(druginfo, druginfo.nature[match(druginfo[,"drug_id"],druginfo.nature[ , "Drug.ID"]), c("Brand.name", "Site.of.screening", "Drug.type", "Drug.class.I", "Drug.class.II", "Target.family", "Effector.pathway.biological.process", "Clinical.trials", "Source")])

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
  ## drug concentration
  message("Read drug concentration")
  drugconc <- read.csv(file.path("/pfs/gdscDrugInfo", "gdsc_drug_concentration.csv"))
  drugconc[!is.na(drugconc) & (drugconc == "" | drugconc == " ")] <- NA
  drugconc <- data.frame("drug.name"=toupper(gsub(badchars, "", drugconc[ , "Compound.Name"])), drugconc)
  if(all(!is.element(drugconc[ , "drug.name"], drugnid[ , "drug.name"]))) { stop("Screening concentration for drugs without identifiers!") }
  myx <- match(drugconc[ , "drug.name"], drugnid[ , "drug.name"])
  ## correct ambiguity for AZD6482: drugid_156 corresponds to the first occurence of AZD6482 while drugid_1066 corresponds to the second
  myx[drugconc[ , "drug.name"] == "AZD6482"][2] <- which(drugnid[ , "drug.name"] == "AZD6482")[2]
  rownames(drugconc) <- rownames(drugnid)[myx]
  drugconc <- data.frame("drugid"=rownames(drugconc), drugconc)

  ## combine all drugs
  dix <- sort(unique(c(rownames(druginfo), rownames(drugconc), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
  ## update druginfo
  druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))))
  newlev <- sapply(druginfo, levels)
  newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
  druginfo2 <- genefu::setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
  druginfo2[match(rownames(druginfo), dix), colnames(druginfo)] <- druginfo
  druginfo2[ , "drugid"] <- newlev$drugid
  druginfo <- druginfo2
  ## update drugconc
  drugconc2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(drugconc), dimnames=list(dix, colnames(drugconc))))
  newlev <- sapply(drugconc, levels)
  newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
  drugconc2 <- genefu::setcolclass.df(df=drugconc2, colclass=sapply(drugconc, class), factor.levels=newlev)
  drugconc2[match(rownames(drugconc), dix), colnames(drugconc)] <- drugconc
  drugconc2[ , "drugid"] <- newlev$drugid
  drugconc <- drugconc2

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
  
  ##
  annot <- read.csv(file.path("/pfs/downAnnotations", "annot_ensembl_all_genes.csv"), stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)
  

  ## normalization of gene expression data
  myf <- "GDSC_U219_ENSG.RData"

  load(file.path("/pfs/gdscU219normalized", myf), verbose=TRUE)

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


##TODO: broken from here down
	  

  myfn2 <- file.path("/pfs/gdscU133anormalized", "GDSC_U133a_ENSG.RData")
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
        
  ## match the experiment labels
  # myx <- rownames(sampleinfo)[match(rownames(genexprs), gsub(".CEL.gz", "", as.character(sampleinfo[ , "filename"])))]
  # genexprs <- genexprs[!is.na(myx), , drop=FALSE]
  # myx <- myx[!is.na(myx)]
  # rownames(genexprs) <- myx

  ## build annotation matrix
  message("Build annotation matrix")
    
    
#   myfn2 <- file.path(saveres, "gdsc_ge_annot.RData")
#   if(!file.exists(myfn2)) {
#     ## select the best probe for a single gene
#     require(jetset) || stop("Library jetset is not available!")
#     js <- jetset::jscores(chip="hgu133a", probeset=colnames(genexprs))
#     js <- js[colnames(genexprs), , drop=FALSE]
#     ## identify the best probeset for each entrez gene id
#     geneid1 <- as.character(js[ , "EntrezID"])
#     names(geneid1) <- rownames(js)
#     geneid2 <- sort(unique(geneid1))
#     names(geneid2) <- paste("geneid", geneid2, sep=".")
#     gix1 <- !is.na(geneid1)
#     gix2 <- !is.na(geneid2)
#     geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
#     ## probes corresponding to common gene ids
#     gg <- names(geneid1)[is.element(geneid1, geneid.common)]
#     gid <- geneid1[is.element(geneid1, geneid.common)]
#     ## duplicated gene ids
#     gid.dupl <- unique(gid[duplicated(gid)])
#     gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
#     ## unique gene ids
#     gid.uniq <- gid[!is.element(gid, gid.dupl)]
#     gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
#     ## which are the best probe for each gene
#     js <- data.frame(js, "best"=FALSE)
#     js[gg.uniq, "best"] <- TRUE
#     ## data for duplicated gene ids
#     if(length(gid.dupl) > 0) {
#       library(jetset)
#       ## use jetset oevrall score to select the best probeset
#       myscore <- js[gg.dupl, "overall"]
#       myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
#       myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
#       myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
#       js[myscore[ , "probe"], "best"] <- TRUE
#     }
#     annot <- data.frame("probe"=rownames(js), "EntrezGene.ID"=js[ , "EntrezID"], js)
#     annot <- annot[colnames(genexprs), , drop=FALSE]
#     save(list=c("annot"), compress=TRUE, file=myfn2)
#   } else { load(myfn2) }
#
	

  ########################
  ## merge data and put NA when the data are anot avaialble
  ########################
   message("Merge data")
  
   ## cell lines and drugs
   # cellnall <- sort(unique(c(as.character(sampleinfo[ , "cellid"]), as.character(drugpheno[ , "cellid"]), rownames(genexprs), rownames(mutation))))
   sampleinfo <- pData(eset)

   cellnall <- sort(unique(c(as.character(sampleinfo[ , "cellid"]), as.character(drugpheno[ , "cellid"]), rownames(mutation))))
   
   drugnall <- sort(unique(c(rownames(druginfo), as.character(drugconc[ , "drugid"]), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
   
   
   
   # ## update sampleinfo
#    dd <- data.frame(matrix(NA, nrow=length(cellnall), ncol=ncol(sampleinfo), dimnames=list(rownames(sampleinfo), colnames(sampleinfo))))
#    newlev <- sapply(sampleinfo, levels)
#    newlev$cellid <- cellnall
#    dd <- genefu::setcolclass.df(df=dd, colclass=sapply(sampleinfo, class), factor.levels=newlev)
#    dd[rownames(sampleinfo), colnames(sampleinfo)] <- sampleinfo
#    dd[ , "cellid"] <- rownames(dd)
#    sampleinfo<- dd
  
   # ## update gene expression
#    dd <- matrix(NA, ncol=ncol(genexprs), nrow=length(cellnall), dimnames=list(cellnall, colnames(genexprs)))
#    dd[rownames(genexprs), colnames(genexprs)] <- genexprs
#    genexprs <- dd
  
   ## update mutation
   dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
   dd[rownames(mutation), colnames(mutation)] <- mutation
   mutation <- dd
   
   MutationEset <- ExpressionSet(t(mutation)) 
   geneMap <- read.csv("/pfs/downAnnotations/annot_ensembl_all_genes.csv")
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



  ## update drugpheno
  ## IC50 in microM
  iix <- grep("_IC_50$", colnames(drugpheno))
  iixn <- gsub("_IC_50$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.ic50 <- exp(dd)
  # ## IC25 in microM
#   iix <- grep("_IC_25$", colnames(drugpheno))
#   iixn <- gsub("_IC_25$", "", colnames(drugpheno)[iix])
#   dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
#   dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
#   drugpheno.ic25 <- exp(dd)
#   ## IC75 in microM
#   iix <- grep("_IC_75$", colnames(drugpheno))
#   iixn <- gsub("_IC_75$", "", colnames(drugpheno)[iix])
#   dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
#   dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
#   drugpheno.ic75 <- exp(dd)
#   ## IC90 in microM
#   iix <- grep("_IC_90$", colnames(drugpheno))
#   iixn <- gsub("_IC_90$", "", colnames(drugpheno)[iix])
#   dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
#   dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
#   drugpheno.ic90 <- exp(dd)
  
  ## AUC (higher value represents sensitivity)
  iix <- grep("_AUC$", colnames(drugpheno))
  iixn <- gsub("_AUC$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.auc <- 1 - dd
  # slope
  iix <- grep("_BETA$", colnames(drugpheno))
  iixn <- gsub("_BETA$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  dd[dd < 0] <- 0
  dd[dd > 1] <- 1
  drugpheno.slope <- dd
  
  # save(drugpheno.auc, drugpheno.slope, file=file.path(saveres, "GDSC_drugpheno.measures.RData"))
  
  nms <- sapply(colnames(drugpheno.ic50), paste, rownames(drugpheno.ic50), sep="_")
  ic50 <- c(drugpheno.ic50)
  names(ic50) <- nms
  auc <- c(drugpheno.auc)
  names(auc) <- nms
  #slope <-c(drugpheno.slope)
  #names(slope) <- nms

  profiles <- cbind("ic50_published"=ic50, "auc_published"=auc)#, "slope_published"=slope)
  # if(cdfname != "brainarray") {
  #   message("Make eset")
  #   fff <- file.path(saveres, sprintf("gdsc_%s.RData", normalization))
  #   if(!file.exists(fff)) {
  #     data.ge <- genexprs
  #     annot.ge <- annot
  #     data.ge <- data.ge[rownames(sampleinfo), rownames(annot.ge), drop=FALSE]
  #     save(list=c("data.ge", "annot.ge", "sampleinfo", "druginfo"), compress=TRUE, file=fff)
  #   } else { load(fff) }
  #   eset <- ExpressionSet(t(data.ge))
  #   pData(eset) <- sampleinfo
  #   fData(eset) <- annot.ge
  #   annotation(eset) <- "rna"
  # }
  
  
  ## update celline
  dd <- data.frame(matrix(NA, ncol=ncol(celline), nrow=length(cellnall), dimnames=list(cellnall, colnames(celline))), check.names=FALSE)
  iix <- intersect(rownames(celline), cellnall)
  dd[iix, colnames(celline)] <- celline[iix, , drop=FALSE]
  celline <- dd
  celline[ , "cell_id"] <- celline[ , "CELL_LINE_NAME"] <- rownames(celline)
  ## annotate cell lines with curated tissue type
  tissue.type <- read.csv(file.path("/pfs/downAnnotations", "cell_annotation_all.csv"), stringsAsFactors=FALSE)
  rownames(tissue.type) <- tissue.type[ , "unique.cellid"]
  celline <- cbind("tissue.type"=tissue.type[match(celline[ , "cell_id"], tissue.type[ , "CGP.cellid"]), "unique.tissueid"], celline)
  
  
  # ## update drugconc
#   drugconcnall <- as.vector(outer(drugnall, cellnall, paste, sep="..."))
#   dd <- data.frame(matrix(NA, nrow=length(drugconcnall), ncol=ncol(drugconc), dimnames=list(drugconcnall, colnames(drugconc))))
#   newlev <- sapply(drugconc, levels)
#   newlev$drugid <- drugnall
#   newlev$cellid <- cellnall
#   dd <- genefu::setcolclass.df(df=dd, colclass=sapply(drugconc, class), factor.levels=newlev)
#   dd[rownames(drugconc), colnames(drugconc)] <- drugconc
#   dd[ , "cellid"] <- sapply(strsplit(x=rownames(dd), split="[.][.][.]"), function (x) { return (x[2]) })
#   dd[ , "drugid"] <- sapply(strsplit(x=rownames(dd), split="[.][.][.]"), function (x) { return (x[1]) })
#   drugconc <- dd
  
  

  
  # myfn2 <- file.path(path.drug, "cpg_concentrations.RData")
#   if (!file.exists(myfn2)){
#     drugconc2 <- array(data=NA, dim=c(length(unique(drugconc[, "drugid"])),length(unique(drugconc[, "cellid"])),3))
#
#     dimnames(drugconc2)[[1]] <- unique(drugconc$drugid)
#     dimnames(drugconc2)[[2]] <- unique(drugconc$cellid)
#     dimnames(drugconc2)[[3]] <- colnames(drugconc)[-c(1:3)]
#
#     # ptm <- proc.time()
#
#     if(verbose){message("Format Drug Concentrations")}
#     for (drug in unique(drugconc$drugid)){
#     for (cell in unique(drugconc$cellid)){
#     drugconc2[drug, cell, 1]<- drugconc[paste(drug, cell, sep="..."),4]
#     drugconc2[drug, cell, 2]<- drugconc[paste(drug, cell, sep="..."),5]
#     drugconc2[drug, cell, 3]<- drugconc[paste(drug, cell, sep="..."),6]
#     }}
#     save(drugconc2, file=myfn2)
#   } else { load(myfn2) }
#
#   drugconc <- drugconc2
  # print(proc.time() - ptm)
  
  ## update druginfo
  dd <- data.frame(matrix(NA, nrow=length(drugnall), ncol=ncol(druginfo), dimnames=list(drugnall, colnames(druginfo))))
  newlev <- sapply(druginfo, levels)
  newlev$drugid <- sapply(strsplit(drugnall, split="_"), function(x) { return(x[2]) })
  dd <- genefu::setcolclass.df(df=dd, colclass=sapply(druginfo, class), factor.levels=newlev)
  dd[match(rownames(druginfo), drugnall), colnames(druginfo)] <- druginfo
  dd[ , "drugid"] <- newlev$drugid
  druginfo <- dd
  druginfo[, "drug.name"] <- druginfo[, "Name"]


  
  
  
  
 
  #raw.sensitivity <- read.csv(file.path(inst("PharmacoGx"), "extdata", "gdsc_sensitivity_detail.csv"))
  getGDSCrawData <-
    function(path.data=file.path("data", "GDSC"), result.type=c("array", "list")){
    require(stringr) || stop("Library stringr is not available!")
    if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
    ## download gdsc raw sensitivity data
    archivn <- "gdsc_drug_sensitivity_raw_data_w5"
    if(!file.exists(file.path(path.data, "dwl"))){dir.create(file.path(path.data, "dwl"), showWarnings=FALSE)}
    myfn <- file.path(path.data, "dwl", sprintf("%s.csv", archivn))
    if (!file.exists(myfn)) {
      dwl.status <- download.file(url=sprintf("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/%s", sprintf("%s.zip", archivn)), destfile=file.path(path.data, "dwl", sprintf("%s.zip", archivn)), quiet=TRUE)
      if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    }
    res <- unzip(zipfile=file.path(path.data, "dwl", sprintf("%s.zip", archivn)),exdir=file.path(path.data, "dwl"))
    gdsc.raw.drug.sensitivity <- read.csv(myfn, stringsAsFactors=FALSE)
    
    concentrations.no <- length(grep("raw", names(gdsc.raw.drug.sensitivity)))
    
    if(result.type == "array"){
      ## create the gdsc.drug.response object including information viablilities and concentrations for each cell/drug pair
      obj <- array(NA,  dim=c(length(unique(gdsc.raw.drug.sensitivity[ , "CELL_LINE_NAME"])),  length(unique(gdsc.raw.drug.sensitivity[ , "DRUG_ID"])),  2,  concentrations.no), dimnames=list(unique(gdsc.raw.drug.sensitivity[ , "CELL_LINE_NAME"]), unique(gdsc.raw.drug.sensitivity[ , "DRUG_ID"]), c("concentration", "viability"), 1:concentrations.no))
    }
    
    fnexperiment <- 
      function(values)  {
        cellline <- values["CELL_LINE_NAME"]
        drug <- stringr::str_trim(values["DRUG_ID"])
        fold_dillution <- as.numeric(values["FOLD_DILUTION"])
        max_dose <- as.numeric(values["MAX_CONC"])
        doses <- rev(as.numeric(c(max_dose, unlist(lapply(1:8, function(x){max_dose/(fold_dillution ^ x)}))))) # micro molar
        if(concentrations.no > length(doses)) {doses <- c(doses, rep(NA, concentrations.no - length(doses)))}
        
        responses <- as.numeric(values[grep("raw",names(values))])
        if(concentrations.no > length(responses)) {responses <- c(responses, rep(NA, concentrations.no - length(responses)))}
        controls <- values[grep("control", names(values))]#qc_fail
        controls <- as.numeric(controls[which(!is.na(controls) & controls != "qc_fail")])
        
        blanks <- values[grep("blank", names(values))]#qc_fail
        blanks <- as.numeric(blanks[which(!is.na(blanks) & blanks != "qc_fail")])
        
        responses <- rev((responses - mean(blanks))/(mean(controls) - mean(blanks))) * 100 #mean can be replaced by median         
        
        if(result.type == "array"){
          obj[cellline,drug, "concentration", 1:length(doses)] <<- doses
          obj[cellline,drug, "viability", 1:length(responses)] <<- responses
        }else{
          return(list(cell=cellline, drug=drug, doses=doses, responses=responses)) #paste(doses, collapse=", "), responses=paste(responses, collapse=", ")))
        }
      }
    
    gdsc.raw.drug.sensitivity.list <- do.call(c, apply(gdsc.raw.drug.sensitivity, 1, list))
    gdsc.raw.drug.sensitivity.res <- mapply(fnexperiment, values=gdsc.raw.drug.sensitivity.list)
    if(result.type == "array"){
      return(list(data=obj, concentrations.no=concentrations.no))
    }else{
      return(list(data =gdsc.raw.drug.sensitivity.res, concentrations.no=concentrations.no))
    }
  }

  load("/pfs/gdscRawSensitivity/GDSC_sens_raw.RData")
  # con_tested <- con_tested
  
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

  recomputed <- readRDS("/pfs/gdscProfiles/gdscProfiles.rds")

  profiles <- cbind(profiles, recomputed[rownames(profiles)])    


  
  	
  cell_all <- read.csv("/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
  rownames(cell_all) <- cell_all[, "unique.cellid"]

  curationCell <- cell_all[which(!is.na(cell_all[ , "CGP.cellid"]) | !is.na(cell_all[, "CGP_EMTAB3610.cellid"])),]
  curationTissue <- cell_all[which(!is.na(cell_all[ , "CGP.cellid"]) | !is.na(cell_all[, "CGP_EMTAB3610.cellid"])),]


	curationCell <- curationCell[ , c("unique.cellid", "CGP.cellid", "CGP_EMTAB3610.cellid")]
	colnames(curationCell) <- gsub("CGP", "GDSC", colnames(curationCell))
	curationTissue <- curationTissue[ , c("unique.tissueid", "CGP.tissueid", "CGP_EMTAB3610.tissueid")]
	colnames(curationTissue) <- gsub("CGP", "GDSC", colnames(curationTissue))
	
	rownames(curationTissue) <- curationCell[, "unique.cellid"]
	rownames(curationCell) <- curationCell[, "unique.cellid"]

	drug_all <- read.csv("/pfs/downAnnotations/drug_annotation_all.csv", na.strings=c("", " ", "NA"))
	curationDrug <- drug_all[ , c("unique.drugid", "CGP.drugid")]
	colnames(curationDrug) <- gsub("CGP", "GDSC", colnames(curationDrug))

	curationDrug <- curationDrug[match(druginfo[, "drug.name"], curationDrug[, "GDSC.drugid"]),]
	mynames <- rownames(druginfo)[match(curationDrug[, "GDSC.drugid"],druginfo[, "drug.name"])]
  	
	mynames[duplicated(mynames)] <- rownames(druginfo)[!rownames(druginfo) %in% mynames]
  rownames(curationDrug) <- mynames
  
  ##update esets pData 
  pData(MutationEset)[, "cellid"] <- rownames(curationCell)[match(rownames(pData(MutationEset)), curationCell[, "GDSC.cellid"])]
  pData(eset)[, "cellid"] <- rownames(curationCell)[match(pData(eset)[, "cellid"], curationCell[, "GDSC.cellid"])]
  pData(gdsc.u219.ensg)[, "cellid"] <- rownames(curationCell)[match(pData(gdsc.u219.ensg)[, "Characteristics.cell.line."], curationCell[, "GDSC_EMTAB3610.cellid"])]

  #integrate cell slot with all the extra new rna cells
  # celline[which(celline[, "cell_id"] == "KMS12-BM"), "cell_id"] <- "KMS-12-BM"
  rownames(celline) <- rownames(curationCell)[match(celline[, "cell_id"], curationCell[, "CGP.cellid"])]
  celline <- cbind(celline, "GDSC_EMTAB3610.cellid"=curationCell[rownames(celline), "GDSC_EMTAB3610.cellid"])
  ##after cbind NA is replaced by emty strings!!!
  celline[!is.na(celline) & celline == ""] <- NA
  tt <- as.data.frame(matrix(NA, ncol=ncol(celline), nrow=length(which(!(rownames(curationCell) %in% rownames(celline))))), stringsAsFactors=FALSE)
  colnames(tt) <- colnames(celline)
  tt[, "GDSC_EMTAB3610.cellid"] <- curationCell[which(!(rownames(curationCell) %in% rownames(celline))), "GDSC_EMTAB3610.cellid"]
  rownames(tt) <- rownames(curationCell)[match(tt[, "GDSC_EMTAB3610.cellid"], curationCell[, "GDSC_EMTAB3610.cellid"])]

  tt[, "cell_id"] <- tt[ , "GDSC_EMTAB3610.cellid"]
  tt[which(tt[ , "cell_id"] == "KMS-12-BM"), "cell_id"] <- "KMS12-BM"
  tt[which(tt[ , "cell_id"] == "IOSE-364-"), "cell_id"] <- "IOSE-364(-)"
  tt[which(tt[ , "cell_id"] == "IOSE-523-"), "cell_id"] <- "IOSE-523(-)"
  tt[which(tt[ , "cell_id"] == "NB-TU-1-10"), "cell_id"] <- "NB(TU)1-10"
  temp <- tt[which(!is.na(match(tt[ , "cell_id"], cosmic[ , "Sample.name"]))),]
  temp[ , colnames(cosmic)] <- cosmic[match(temp[, "cell_id"], cosmic[ , "Sample.name"]),]
  tt[rownames(temp), ] <- temp
  tt[, "tissue.type"] <- cell_all[rownames(tt), "CGP_EMTAB3610.tissueid"]
  celline <- rbind(celline,tt)
  celline[, "tissueid"] <- curationTissue[rownames(celline), "unique.tissueid"]
  celline[, "cellid"] <- curationTissue[rownames(celline), "unique.cellid"]


  ##update sensitivity info cellid to be unique.cellid
  sensitivityInfo[, "cellid"] <- rownames(curationCell)[match(sensitivityInfo[, "cellid"], curationCell[, "GDSC.cellid"])]

  ##update drugid to be unique dugid (there is problem with AZD6482)
  tt <- NULL
  for(drug in unique(druginfo[, "drug.name"])){
    hits <- which(druginfo[, "drug.name"] == drug)
    if(length(hits) > 1)
    {
      hit <- hits[which.min(apply(druginfo[hits,], 1, function(x){length(which(is.na(x)))}))]
      druginfo[hits, "drugid"] <- paste(druginfo[hits, "drugid"], collapse= "///")
      druginfo[hits, "Last.screening.date"] <- paste(druginfo[hits, "Last.screening.date"], collapse= "///")
    }else{
      hit <- hits
    }
    tt <- rbind(tt, druginfo[hit, ])
  }
  druginfo <- as.data.frame(tt, stringAsFactors=FALSE)
  rownames(druginfo) <- curationDrug[match(rownames(druginfo), rownames(curationDrug)), "unique.drugid"]
  sensitivityInfo[, "drugid"] <- curationDrug[match(sensitivityInfo[, "drugid"], rownames(curationDrug)), "unique.drugid"]
  curationDrug <- curationDrug[which(!duplicated(curationDrug[, "unique.drugid"])),]
  rownames(curationDrug) <- curationDrug[, "unique.drugid"]
  ###clear drug info columns
  druginfo <- druginfo[, -3]
  druginfo <- druginfo[, -3]

  curationTissue <- cell_all[rownames(celline), c("unique.tissueid", "CGP.tissueid", "CGP_EMTAB3610.tissueid")]
  rownames(curationTissue) <- rownames(celline)

  ## fusion genes 
  ## TODO:: is this not just drugpheno?
  myf <- "/pfs/gdscDrugInfo/dwl/gdsc_manova_input_w5.csv"

  gdsc_manova_input <- read.csv(myf, stringsAsFactors=FALSE, header=TRUE, na.strings=c("", " "))
  gdsc_manova_input <- gdsc_manova_input[7:nrow(gdsc_manova_input),]
  start.index <- which(colnames(gdsc_manova_input) == "BCR_ABL")
  gdsc.fusion <- gdsc_manova_input[,start.index:(start.index+2)]
  rownames(gdsc.fusion) <- gdsc_manova_input[, "Cell.Line"]
  rownames(gdsc.fusion) <- rownames(curationCell)[match(rownames(gdsc.fusion), curationCell$GDSC.cellid)]

  fusion.eset <- ExpressionSet(t(gdsc.fusion))
  pData(fusion.eset)[ ,"batchid"] <- NA
  pData(fusion.eset)[ ,"cellid"] <- rownames(pData(fusion.eset))
  annotation(fusion.eset) <- "fusion"

  ##manual curation//not obvious why it's missed!
  pData(eset)[which(pData(eset)[,"Characteristics.CellLine."] == "KMS-12-BM"),"cellid"] <- "KMS-12-BM"
  ####CNV
  cell_all <- read.csv(system.file("extdata", "cell_annotation_all.csv", package="PharmacoGx"), na.strings=c("", " ", "NA"))
  
  ## TODO: missing gdsc.cnv.eset.RData

  # load(file.path(path.cnv, "gdsc.cnv.eset.RData"))
  # tt <- rownames(pData(gdsc.eset))
  # pData(gdsc.eset) <- as.data.frame(apply(pData(gdsc.eset), MARGIN=2, as.character), stringsAsFactors=FALSE)
  # rownames(pData(gdsc.eset)) <- tt
  # pData(gdsc.eset)[ , "cellid"] <- cell_all$unique.cellid[match(pData(gdsc.eset)[, "Name"], cell_all$GDSC.SNP.cellid)]
  # pData(gdsc.eset)[,"batchid"] <- NA
  # annot <- read.csv(system.file("extdata", "annot_ensembl_all_genes.csv", package="PharmacoGx"), stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)
  # tt <- annot[match(rownames(fData(gdsc.eset)), annot$gene_name), c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")]
  # rownames(tt) <- rownames(fData(gdsc.eset))
  # colnames(tt) <- c("EnsemblGeneId", "EntrezGeneId", "Symbol", "GeneBioType")
  # fData(gdsc.eset) <- tt
  # annotation(gdsc.eset) <- "cnv"
  
#   myf2 <- "saveres/gdsc_csv_exprs.RData"
#   if(!file.exists(myf2)) {
#     tt <- exprs(gdsc.eset)    
#     tt <- apply(tt, MARGIN=c(1,2), function(x){as.character(max(as.numeric(unlist(strsplit(x,split="[/]")))))})
#     temp <- rownames(tt)
#     tt <- apply(tt, MARGIN=2, as.numeric)
#     rownames(tt) <- temp
#     save(tt, file=myf2)
#   }else{
#     load(myf2)
#   }
  # tt <- exprs(gdsc.eset)  
  # temp <- rownames(tt)
  # tt <- apply(tt, MARGIN=2, as.numeric)
  # rownames(tt) <- temp
  # exprs(gdsc.eset) <- tt

  # data <- exprs(gdsc.eset)
  # tt <- apply(data, MARGIN = 1, function(x){table(is.na(x))["FALSE"]})
  # features <- setdiff(rownames(exprs(gdsc.eset)), names(tt)[which(is.na(tt))])
  # gdsc.eset <- gdsc.eset[features, ]

  curationCell <- cell_all[which(!is.na(cell_all$CGP.cellid) | !is.na(cell_all$CGP_EMTAB3610.cellid) | !is.na(cell_all$GDSC.SNP.cellid)), c("unique.cellid", "CGP.cellid", "CGP_EMTAB3610.cellid", "GDSC.SNP.cellid")]
  rownames(curationCell) <- curationCell$unique.cellid
  colnames(curationCell) <- gsub("CGP", "GDSC", colnames(curationCell))
  
  curationTissue <- cell_all[which(!is.na(cell_all$CGP.cellid) | !is.na(cell_all$CGP_EMTAB3610.cellid) | !is.na(cell_all$GDSC.SNP.cellid)), c("unique.tissueid", "CGP.cellid", "CGP_EMTAB3610.cellid", "GDSC.SNP.tissueid")]
  rownames(curationTissue) <- curationCell$unique.cellid
  colnames(curationTissue) <- gsub("CGP", "GDSC", colnames(curationTissue))

  tt <- matrix(NA, ncol=ncol(celline), nrow=length(setdiff(curationCell$unique.cellid, rownames(celline))))
  colnames(tt) <- colnames(celline)
  rownames(tt) <- setdiff(curationCell$unique.cellid, rownames(celline))
  celline <- rbind(celline, tt)
  ##update celline
  celline <- celline[rownames(curationCell), ]
  celline[ ,"tissueid"] <- curationTissue[rownames(celline), "unique.tissueid"]
  celline[ ,"cellid"] <- rownames(celline)
  

  GDSC <- PharmacoSet(molecularProfiles=list("rna"=eset, "rna2"=gdsc.u219.ensg, "mutation"=MutationEset, "fusion"=fusion.eset),#, "cnv"=gdsc.eset),
                      name="GDSC", 
                      cell=celline, 
                      drug=druginfo, 
                      sensitivityInfo=sensitivityInfo, 
                      sensitivityRaw=raw.sensitivity, 
                      sensitivityProfiles=profiles, 
                      sensitivityN=NULL, 
                      curationCell=curationCell, 
                      curationDrug=curationDrug, 
                      curationTissue=curationTissue, 
                      datasetType="sensitivity")
  
  # myfn <- file.path(saveres, "GDSC_sens_tables.RData")
#   if(!file.exists(myfn)){
# 	  tables <- array(NA, c(nrow(celline), nrow(druginfo), ncol(profiles)), dimnames=list(rownames(celline), rownames(druginfo), colnames(profiles)))
# 	  for (measure in colnames(profiles)){
# 	  	  if(verbose){
# 			  print("Now generating table for: ")
# 			  print(measure)
# 	  	  }
# 		  tables[,,measure] <- summarizeSensitivityprofiles(gdsc, measure=measure, summary="median")
#
# 	  }
# 	  save(tables, file=myfn)
#   } else {
# 	  load(myfn)
#   }
#   gdsc@tables <- tables
#   gdsc@table.summary[colnames(profiles)] <- "median"
  
  # gdsc@exprs_cell_id <- "cellid"
#
#   gdsc@sens_cell_id <- "cellid"
#
#   gdsc@sens_drug_id <- "drugid"
#
 
 
  ##JetSet probe/gene mapping
  # save(GDSC, file="../PSets/GDSC.RData")
 return (GDSC)
  
}
GDSC <- getGDSC()
save(GDSC, file="/pfs/out/GDSC.RData")
## End
