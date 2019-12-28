#!/usr/bin/Rscript

library (stringr)
library (dplyr)


#------------------------------------------------------------------------
#------------------------------------------------------------------------
createVennDiagramsMarkers <- function (A,B,C,D){
	x <- list()
	x$Gwasp4 = as.character (A)
	x$Gwasp2 = as.character (B)
	x$Plink  = as.character (C)
	x$Tassel = as.character (D)

	require(VennDiagram)
	v0 <<-venn.diagram(x, height=9000, width=9000,
										col = c("red", "blue", "green", "yellow"),
										fill = c("red", "blue", "green", "yellow"), 
										alpha = 0.5, filename = NULL)
	grid.draw(v0)
	overlaps <- calculate.overlap(x)
	overlaps <- rev(overlaps)

	posOverlap = as.numeric (gsub ("a","", (names (overlaps))))
	for (i in 1:length(overlaps)){
		pos = posOverlap [i]
		v0[[pos+8]]$label <- paste(overlaps[[i]], collapse = "\n")
	}

	pdf("venn.pdf")
	grid.draw(v0)
	dev.off()
}


#------------------------------------------------------------------------
#------------------------------------------------------------------------

options (width=300)
#options(scipen=999)
nROWS = 5

model="Naive"

files =  list.files(".", pattern="^(.*(Naive).*(tbl)[^$]*)$")

files = c("out-Gwasp4-Naive-significativeQTLs.tbl", "out-Gwasp2-Naive-significativeQTLs.tbl", 
		  "out-Plink-Naive-assoc.linear.adjusted.tbl", "out-Tassel-Naive-GLM_Stats_geno+pheno.tbl")
#files = c("out-Gwasp4-Naive-significativeQTLs.tbl", 
summTable = data.frame ()
pThreshold = round (1/(578*2068), 10)
for (f in files) {
	data = read.table (file=f, header=T)
	print (f)
	if (str_detect(f, "Gwasp")) {
		if (str_detect(f, "Gwasp4")) tool = "Gwasp4" else tool = "Gwasp2"
		data    = data [1:nROWS,]
		snps    = data$Marker
		pVal	  = round (10^(-data$Score),10)
		chrom   = data$Chrom
		pos	    = data$Position
		signf   = data$Score >= data$Threshold
		pscores = data$Score
		tscores = round (-log10 (pThreshold), 4)
	}else if (str_detect (f, "Plink")) {
		data    = data [1:nROWS,]
		tool    = "Plink"
		snps    = data$SNP
		pVal    = data$UNADJ
		chrom   = data$CHR
		pos	    = NA
		pscores = round (-log10 (pVal), 4)
		tscores = round (-log10 (pThreshold), 4)
		signf   = pscores >= tscores
		
	}else if (str_detect (f, "Tassel")) {
		data = data %>% top_n (-1*nROWS, p)
		tool   = "Tassel"
		snps   = data$Marker
		pVal   = data$p
		chrom  = data$Chr
		pos	= data$Pos
		pscores  = round (-log10 (pVal), 4)
		tscores = round (-log10 (pThreshold), 4)
		signf    = pscores >= tscores
	}

	df = data.frame (TOOL=tool, MODEL=model, CHR=chrom, POS=pos, SNP=snps, P = pVal, pTHR=pThreshold, pSCORE=pscores, tSCORE=tscores, Signf=signf )
	df = df %>% distinct (SNP, .keep_all=T)
	#df = data.frame (TOOL=tool, MODEL=model, SNPs=snps, P = pVal, CHR=chrom, POS=pos,pTHR=pThreshold)
	summTable = rbind (summTable, df)
	#summTableSorted = summTable %>% arrange (SNPs, P)
	summTableSorted = summTable %>% add_count (SNP, sort=T, name="N1") %>% arrange (desc(N1))
	summTableSorted = summTableSorted %>% add_count (SNP, Signf, sort=T, name="Ns") %>% arrange (desc(Ns))
	write.table (file="summary-gwas.tbl", summTable, row.names=F,quote=F, sep="\t")
	write.table (file="summary-gwas-sorted.tbl", summTableSorted, row.names=F,quote=F, sep="\t")

}
# Create Venn diagram of common markers
markersGwasp4 = summTable %>% filter (TOOL %in% "Gwasp4") %>% select (SNP) #%>% as.character
markersGwasp2 = summTable %>% filter (TOOL %in% "Gwasp2") %>% select (SNP) #%>% as.character
markersPlink  = summTable %>% filter (TOOL %in% "Plink") %>% select (SNP)  #%>% as.character
markersTassel = summTable %>% filter (TOOL %in% "Tassel") %>% select (SNP) #%>% as.character

#createVennDiagrams (markersGwasp4$SNP, markersGwasp2$SNP, markersPlink$SNP)
createVennDiagramsMarkers (markersGwasp4$SNP, markersGwasp2$SNP, markersPlink$SNP, markersTassel$SNP)
