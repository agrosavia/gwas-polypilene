#!/usr/bin/Rscript
# r4.1: Without structure file input, Working with G4,G2,PL,SH
# r2.41: Changed to parameters by config file. Beginning changes...
# r2.3: Removed annotations as it fails. Before config by file instead parameters
# r2.2: Write model type and inflation factor to significative QTLs
# r2.1: Support for multiple command line arguments (dynamic)

LOAD_DATA = T
args = commandArgs(trailingOnly = TRUE)
#args = c("in/config-gwaspoly-Kinship+PCs.config")

USAGE="USAGE: Rscript gwas-polypiline.R <config file>"
if (length (args) != 1) {
	message (USAGE)
	quit()
}
#-------------------------------------------------------------
#-------------------------------------------------------------

library (GWASpoly) #
suppressMessages(library (config))  # For read config file

setClass ("GWASpolyStruct", slots=c(params="list"),#testModels="character" snpModels="character", 
		  contains="GWASpoly.K")

options (width=300)
source ("gwas-formats.R")           # Module with functions to convert between different genotype formats 

#-------------------------------------------------------------
# Global configs
#-------------------------------------------------------------
#genotypeFormats = "numeric"|"AB"|"ACGT"
#gwasModel = "Naive"|"Kinship+PCs"
#-------------------------------------------------------------
data = data1 = data2 = data3 = data4 = NULL
phenotype = genotype = NULL
dt=NULL
#-------------------------------------------------------------
# Main
#-------------------------------------------------------------
main <- function (args) 
{
	# Read and check config file arguments
	configFile =  args [1]
	params = getConfigurationParameters (configFile)
	genotypeFile  = params$genotypeFile
	phenotypeFile = params$phenotypeFile

	msg (">>>>>>>>>>>>", params$gwasModel, "<<<<<<<<<<<") 

	# Read, filter, and check phenotype and genotype
	data <- dataPreprocessing (genotypeFile, phenotypeFile, params)
	genoGwas2File   = "out/filtered-gwasp2-genotype.tbl"
	phenoGwasp2File = "out/filtered-gwasp4-phenotype.tbl"

	runGwaspolyGwas (phenotypeFile, genotypeFile, 4, params, data)
	runGwaspolyGwas (phenoGwasp2File, genoGwas2File, 2, params, data)

	runPlinkGwas (genotypeFile, phenotypeFile, params$gwasModel)
	runShesisGwas ()
}

#-------------------------------------------------------------
#-------------------------------------------------------------
runGwaspolyGwas <- function (phenotypeFile, genotypeFile, ploidy,  params, data) 
{
	if (LOAD_DATA & file.exists ("gwas.RData")) load (file="gwas.RData")

	if (ploidy==4) {
		#snpModels=testModels = ("general")
		snpModels  = c("general","additive","1-dom", "2-dom")
		testModels = c("general", "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
	}else if (ploidy==2) {
		#snpModels=testModels = ("general")
		snpModels  = c("general","additive","1-dom")
		testModels = c("general", "additive","1-dom-alt","1-dom-ref")
	}else 
		stop ("unknown ploidy number (only 2 or 4)")

	params = append (params, list (snpModels=snpModels, testModels=testModels))

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
	data1 <- initGWAS (phenotypeFile, genotypeFile, ploidy, params$genotypeFormat, data1)

	# Set the kinship
	data2 <- setKinship (data1, params$gwasModel, data2)

	# Populations structure and kinship
	data3 <- setPopulationStructure (data2, params$gwasModel, data$phenotype, data3)

	# GWAS execution
	data4 <- runGwaspoly (data3, params$gwasModel, params$snpModels, data4)
	showResults (data4, params$testModels, data$trait, params$gwasModel, 
				 params$phenotypeFile, params$annotationsFile, ploidy)

	if (LOAD_DATA) save(data, data1, data2, data3, data4, file="gwas.RData") 
}

#-------------------------------------------------------------
# 
#-------------------------------------------------------------
runPlinkGwas <- function (genotypeFile, phenotypeFile, model) 
{
	msg ("Running plink gwas..", model)
	if (model=="Naive") {
		cmm = paste0 ("plink --file out/filtered-plink-genotype --linear --assoc --adjust --pheno out/filtered-plink-phenotype.tbl --all-pheno --allow-no-sex --out out/out-plink-", model)
		runCommand (cmm)
	}else if (model=="Kinship"|model=="Kinship+PCs") {
		runCommand ("plink --file out/filtered-plink-genotype --make-bed --out out/tmp-plinkb")
		runCommand ("plink2 --bfile out/tmp-plinkb --make-pgen --out out/tmp-plinkg")
		runCommand ("plink2 --pfile out/tmp-plinkg --king-cutoff 0.177 --out out/tmp-plinkg")

		msg ("Writing new filtered kinship genotype...")
		outKinship      = read.table ("out/tmp-plinkg.king.cutoff.out.id")
		plinkGeno       = read.csv (file="out/filtered-plink-genotype.ped", sep="\t", header=F)
		markers         = outKinship [,1]
		markersFiltered = plinkGeno %>% filter (V2 %in% markers)
		write.table (file="out/filtered-plink-genotype.ped", markersFiltered, quote=F, sep="\t", row.names=F, col.names=F)

		runCommand (paste0 ("plink --file out/filtered-plink-genotype --linear --assoc --adjust --pheno out/filtered-plink-phenotype.tbl --all-pheno --allow-no-sex --out out/out-plink-", model))
	}
}
	
runShesisGwas <- function () 
{
	inGenoPheno = "out/filtered-shesis-genopheno.tbl"
	inMarkers   = "out/filtered-shesis-markernames.tbl"
	outFile     = "out/out-shesis-naive"
	cmm=sprintf ("shesis --input %s --ploidy 4 --assoc --qtl --snpname-file %s --report-txt --adjust  --output %s", 
				 inGenoPheno, inMarkers, outFile)
	msg (cmm)
	system (cmm)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
runCommand <- function (command, logFile="gwas.log") {
	msg (">>>> ", command)
	if (logFile != "")
		system (paste0 (command, " >> ", logFile))
	else
		system (command)
}

#-------------------------------------------------------------
# Filters the genotype by different quality control filters
# Read and check files, sample samples
#-------------------------------------------------------------
dataPreprocessing <- function (genotypeFile, phenotypeFile, params) 
{
	# Format convertion from gwasp4 to plink2
	msg();msg ("Converting gwaspoly to plink formats...")
	markersIdsMap = gwasp2plinkGenoMap (genotypeFile)
	plinkFile = gwaspTetraGenoToPlinkPed (genotypeFile, markersIdsMap)

	# Apply filters to genotype (markers and samples) by calling external program
	filtered = genotypeFiltering (plinkFile, params)

	# Filter phenotype
	phenoAll = read.csv (phenotypeFile, header=T, sep=",", check.names=F)
	phenoFiltered  = phenoAll [phenoAll[,1] %in% filtered$Samples,]

	# Filter genotype
	genoAll  = read.csv (genotypeFile, header=T, sep=",", check.names=F)
	rownames (genoAll) = genoAll [,1]
	genoFiltered <- genoAll [filtered$Markers,]

	# Select common names geno and pheno
	samplesNamesGeno   <- colnames (genoFiltered)
	samplesNamesPheno  <- phenoFiltered[,1]
	commonSamples <- Reduce (intersect, list (samplesNamesGeno, samplesNamesPheno)) 

	# Filter geno and pheno by common names
	phenoCommon <- phenoFiltered [phenoFiltered[,1] %in% commonSamples,]
	genoColums  <- c(samplesNamesGeno[1:3], commonSamples)
	genoCommon  <- genoFiltered  [,colnames(genoFiltered) %in% genoColums]

	outGenoFile  <- "out/filtered-gwasp4-genotype.tbl"
	outPhenoFile <- "out/filtered-gwasp4-phenotype.tbl"
	write.table (file=outGenoFile, genoCommon, row.names=F, quote=F, sep=",")
	write.table (file=outPhenoFile, phenoCommon, row.names=F, quote=F, sep=",")

	trait  <- colnames (phenoCommon)[2]
	msg  (">>>> Evaluating trait ", trait)

	# Convert and write plink pheno filtered
	msg (outPhenoFile)
	gwasp2plinkPhenotype  (outPhenoFile,"out/filtered-plink-phenotype.tbl") 
	gwasp2tasselPhenotype (outPhenoFile,"out/filtered-tassel-phenotype.tbl") 
	gwasp2shesisGenoPheno ("out/filtered-gwasp4-genotype.tbl", "out/filtered-gwasp4-phenotype.tbl") 

	return (list (genotypeFile=outGenoFile, phenotypeFile=outPhenoFile, trait=trait))

}
#-------------------------------------------------------------
# Apply filters to genotype (markers and samples) by calling external program
#-------------------------------------------------------------
genotypeFiltering <- function (plinkFile, params) {
	# Filter missingness per sample (MIND)"
	MIND = params$MIND
	# Filter missingness per SNP    (GENO)
	GENO = params$GENO
	# Filter SNPs with a low minor allele frequency (MAF)
	MAF  = params$MAF
	# Filter SNPs which are not in Hardy-Weinberg equilibrium (HWE).
	HWE  = params$HWE
	# Filter individuals closely related (Cryptic relatedness)  : --genome --min 0.0001
	runCommand (sprintf ("plink --file %s --mind %s --geno %s --maf %s --hwe %s --out %s-QC --make-bed", 
						 plinkFile,MIND,GENO,MAF,HWE,plinkFile))
	runCommand (sprintf ("plink --bfile %s-QC --recode tab --out %s-QC", plinkFile, plinkFile))

	# Copy links of filtered plink files to main dir
	runCommand (sprintf ("ln -s %s/%s-QC.ped out/filtered-plink-genotype.ped", getwd(), plinkFile))
	runCommand (sprintf ("ln -s %s/%s-QC.map out/filtered-plink-genotype.map", getwd(), plinkFile))

	# Get final markers and individuals"
	sampleNames = read.table (paste0(plinkFile,"-QC.ped"))[,2]
	write.table (file="out/filtered-names-samples.tbl", sampleNames, quote=F,row.names=F,col.names=F)
	markerNames = read.table (paste0(plinkFile,"-QC.map"))[,2]
	write.table (file="out/filtered-names-markers.tbl", markerNames, quote=F,row.names=F,col.names=F)

	msg ("Reading filtered markers and individuals...")
	filteredMarkers <- read.table (file="out/filtered-names-markers.tbl", stringsAsFactors=F)[,1]
	filteredSamples <- read.table (file="out/filtered-names-samples.tbl", stringsAsFactors=F)[,1]

	return (list (Markers=filteredMarkers, Samples=filteredSamples))
}

#-------------------------------------------------------------
# Set the kinship matrix using the realized relationship matrix  
#-------------------------------------------------------------
setKinship <- function (data1,  gwasModel, data2) 
{
	msg();msg("Setting kinship...")
 	# Load data instead calculate it
	if (!is.null (data2)) {msg (">>>> Loading kinship..."); return (data2) }

	if (gwasModel %in% c("Kinship", "Kinship+Structure", "Kinship+PCs")) {
		msg ("    >>>> With default kinship... ")
		kinshipMatrix = NULL
		data2  = set.K (data1)
	}else { 
		msg ("    >>>> Without kinship...") 
		markerNames   = data1@pheno [,1]
		n             = length (markerNames)
		data2         = set.K (data1, K=NULL)
		#kinshipMatrix = matrix (diag (n), n, n, dimnames=list (markerNames, markerNames))
	}		
	return (data2)
}
		
#-------------------------------------------------------------
# Fix Populations structure and kinship
#-------------------------------------------------------------
setPopulationStructure <- function (data2, gwasModel, phenotype, data3) 
{
	st <<- structure
	msg();msg ("Setting population structure...")
	#setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")
	data3 <- new ("GWASpolyStruct", data2)
	
	nPCs=0
	structNames = structTypes = NULL

	if (gwasModel %in% c("PCs","Kinship+PCs")) {
		nPCs=5
		msg (">>>> Population structure with nPCS=", nPCs, "...")
	}else 
		msg (">>>> Without population structure")

	data3@params = set.params(n.PC=nPCs, fixed=NULL, fixed.type=NULL)

	return (data3)
}

#-------------------------------------------------------------
# GWAS execution
#-------------------------------------------------------------
runGwaspoly <- function (data3, gwasModel, snpModels, data4) 
{
	msg();msg("Running GWASpoly...")

	if (!is.null (data4)) { msg (">>>> Loading GWASpoly..."); return (data4) }

 	if (gwasModel %in% c("Naive","Kinship")) {
		msg (">>>> Without params")
		data4 = GWASpoly(data3, models=snpModels, traits=NULL, params=NULL, n.core=4)
	}else {
		msg (">>>> With params")
		data4 = GWASpoly(data3, models=snpModels, traits=NULL, params=data3@params)
	}
	
	return (data4)
}

#-------------------------------------------------------------
# Plot results
#-------------------------------------------------------------
showResults <- function (data4, testModels, trait, gwasModel, phenotypeFile, snpsAnnFile, ploidy) 
{
	msg();msg ("Writing results...", trait)

	msg (">>>> Plotting results...")
	phenoName = strsplit (phenotypeFile, split=".tbl")[[1]][1]
	plotName = sprintf("out/out-gwasp%s-%s-plots.pdf", ploidy, gwasModel)

	n = length (testModels)
	
	# QTL Detection
	data5 = set.threshold (data4, method="Bonferroni",level=0.05,n.permute=100,n.core=4)

	# Plots
	pdf (file=plotName, width=11, height=7)
	# QQ plot 
	op <- par(mfrow = c(2,n), oma=c(0,0,3,0))
	for (i in 1:length(testModels)) {
		#par (cex.main=0.5, cex.lab=0.5, cex.axis=0.5, ann=T)
		qqPlot(data4,trait=trait, model=testModels[i], cex=0.3)
	}

	# Manhattan plot 
	for (i in 1:length(testModels)) {
		#par (cex=1.5)
		#manhattan.plot (y.max=20,data5, trait=trait, model=testModels [i])
		manhattan.plot (data5, trait=trait, model=testModels [i])
	}
	plotTitle = sprintf ("%s gwas %s-ploidy for %s trait", gwasModel, ploidy, trait)  
	mtext(plotTitle, outer=T,  cex=1.5,  line=0)
	par(op)
	dev.off()

	msg (">>>> Writing QTLs...")
	significativeQTLs = getQTL (data5, snpsAnnFile, gwasModel, ploidy)
	outFile = sprintf ("out/out-gwasp%s-%s-QTLs.tbl", ploidy, gwasModel) 
	write.table (file=outFile, significativeQTLs, quote=F, sep="\t", row.names=F)
}

#-------------------------------------------------------------
# Extracts significant QTL
#-------------------------------------------------------------
getQTL <- function(data,snpsAnnFile, gwasModel, ploidy, traits=NULL,models=NULL) 
{
	stopifnot(inherits(data,"GWASpoly.thresh"))
	sc <<-data@scores

	if (is.null(traits)) traits <- names(data@scores)
	else stopifnot(is.element(traits,names(data@scores)))

	if (is.null(models)) models <- colnames(data@scores[[1]])
	else stopifnot(is.element(models,colnames(data@scores[[1]])))

	if (!is.null (snpsAnnFile)) {
		msg (">>>> Reading associations...")
		snpsAnnotations <- read.csv (file=snpsAnnFile, header=T)
	}

	n.model <- length(models)
	n.trait <- length(traits)
	output <- data.frame(NULL)
	for (j in 1:n.model) {
		ix <- which(data@scores[[traits[1]]][,models[j]] > data@threshold[traits[1],models[j]])
		markers <-  data.frame (SNP=data@map[ix,c("Marker")])
		msg ("         ", models [j],": ", as.character (markers[,1]))

		if (!is.null (snpsAnnFile)) 
			snpAnn  <- merge (markers, snpsAnnotations, by.x="SNP",by.y="SNP_id", sort=F)[,c(2,7)]
		else
			snpAnn = "None"

		scores <- data@scores[[1]][,models[j]]
		datax = calculateInflationFactor (scores)

		n.ix <- length(ix)
		
		df = data.frame(Ploidy=rep (ploidy, n.ix),
						Type=rep (gwasModel, n.ix),
						GC=rep(datax$delta,n.ix), 
						Model=rep(models[j],n.ix),
						Score=round(data@scores[[traits[1]]][ix,models[j]],2),
						Threshold=round(rep(data@threshold[traits[1],models[j]],n.ix),2),
						Effect=round(data@effects[[traits[1]]][ix,models[j]],2),
						data@map[ix,])
						#snpAnn) 
						#stringsAsFactors=F,check.names=F)

		output <- rbind(output, df)
	}
	#out <-cbind (Type=gwasModel, output)
	output = output [order(-output$GC,-output$Score),]
	return(output)
}
#-------------------------------------------------------------
# QQ plot
#-------------------------------------------------------------
qqPlot <- function(data,trait,model,cex=1,filename=NULL) 
{
	stopifnot(inherits(data,"GWASpoly.fitted"))
	sc <<- data@scores
	traits <- names(data@scores)

	stopifnot(is.element(trait,traits))
	models <- colnames(data@scores[[trait]])
	stopifnot(is.element(model,models))
	scores <- data@scores[[trait]][,model]

	datax = calculateInflationFactor (scores)

	n <- length(datax$scores)
	unif.p <- -log10(ppoints(n))
	if (!is.null(filename)) {postscript(file=filename,horizontal=FALSE)}
	par(pty="s")
	plot(unif.p, datax$scores, pch=16,cex=cex,
		 xlab=expression(paste("Expected -log"[10],"(p)",sep="")),
		 ylab=expression(paste("Observed -log"[10],"(p)",sep="")),
		 main=paste(trait," (",model,") ",sep=""))

	mtext (bquote(lambda[GC] == .(datax$delta)), side=3, line=-2, cex=0.7)

	lines(c(0,max(unif.p)),c(0,max(unif.p)),lty=2)
	if (!is.null(filename)) {dev.off()}
	return(datax$delta)
}

#-------------------------------------------------------------
# Calculate the inflation factor from -log10 values
#-------------------------------------------------------------
calculateInflationFactor <- function (scores)
{
	remove <- which(is.na(scores))
	if (length(remove)>0) 
		x <- sort(scores[-remove],decreasing=TRUE)
	else 
		x <- sort(scores,decreasing=TRUE)

	pvalues = 10^-x
	chisq <- na.omit (qchisq(1-pvalues,1))
	delta  = round (median(chisq)/qchisq(0.5,1), 3)

	return (list(delta=delta, scores=x))
}

#-------------------------------------------------------------
# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
#-------------------------------------------------------------
initGWAS <- function (phenotypeFile, genotypeFile, ploidy, format, data1) 
{
	msg();msg ("Initializing GWAS...");msg()
	# When data is previously loaded
	if (!is.null (data)) {msg(">>>> Loading GWAS data..."); return (data)}

	data1 <- read.GWASpoly (ploidy = ploidy, pheno.file = phenotypeFile, 
							geno.file = genotypeFile, format = format, n.traits = 1, delim=",")

	return (data1)
}

#-------------------------------------------------------------
# Get params from config file and define models according to ploidy
#-------------------------------------------------------------
getConfigurationParameters <- function (configFile) 
{
	msg("Reading config file...")

	params = config::get (file=configFile) 

	message ("------------------------------------------------")
	message ("Summary of input parameters:")
	message ("------------------------------------------------")
	msg ("Phenotype filename : ", params$phenotypeFile) 
	msg ("Genotype filename  : ", params$genotypeFile) 
	msg ("GwAS model         : ", params$gwasModel) 
	msg ("Genotype format    : ", params$format) 
	msg ("MIND               : ", params$MIND) 
	msg ("GENO               : ", params$GENO) 
	msg ("MAF                : ", params$MAF) 
	msg ("HWE                : ", params$HWE) 
	message ("------------------------------------------------")

	runCommand (sprintf ("cp %s out/", configFile))

	return (params)
}
#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

#-------------------------------------------------------------
# Call main 
#-------------------------------------------------------------
main (args)


