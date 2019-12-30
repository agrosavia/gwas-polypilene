#!/usr/bin/Rscript
## Log: r0.8: Create venn diagrams, summary table of first Ns"
##     

library (stringr)
library (dplyr)

#------------------------------------------------------------------------
# Create Venn diagram of common markers using info from summary table
#------------------------------------------------------------------------
createVennDiagramsMarkers <- function (summaryTable){
	require(VennDiagram)
	flog.threshold(ERROR)
	x <- list()
	x$Gwasp4 = summaryTable %>% filter (TOOL %in% "Gwasp4") %>% select (SNP) %>% .$SNP
	x$Gwasp2 = summaryTable %>% filter (TOOL %in% "Gwasp2") %>% select (SNP) %>% .$SNP
	x$Plink  = summaryTable %>% filter (TOOL %in% "Plink")  %>% select (SNP) %>% .$SNP
	x$Tassel = summaryTable %>% filter (TOOL %in% "Tassel") %>% select (SNP) %>% .$SNP

	v0 <<-venn.diagram(x, height=9000, width=9000,
										col = c("red", "blue", "green", "yellow"),
										fill = c("red", "blue", "green", "yellow"), 
										alpha = 0.5, filename = NULL)

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
SIGNIFICANCE_LEVEL = 0.05
nMARKERS = 2068
THRESHOLD = round (-log10 (SIGNIFICANCE_LEVEL/nMARKERS),4)

model="Naive"
files =  list.files(".", pattern="^(.*(Naive).*(tbl)[^$]*)$")
files = c("out-Gwasp4-Naive-significativeQTLs.tbl", "out-Gwasp2-Naive-significativeQTLs.tbl", 
		  "out-Plink-Naive-assoc.linear.adjusted.tbl", "out-Tassel-Naive-GLM_Stats_geno+pheno.tbl")
#files = c("out-Gwasp4-Naive-significativeQTLs.tbl", 
summaryTable = data.frame ()

for (f in files) {
	data = read.table (file=f, header=T)
	print (f)
	if (str_detect(f, "Gwasp")) {
		if (str_detect(f, "Gwasp4")) tool = "Gwasp4" else tool = "Gwasp2"
		data    = data [1:nROWS,]
		snps    = data$Marker
		pVal	= round (10^(-data$Score),10)
		chrom   = data$Chrom
		pos	    = data$Position
		pscores = data$Score
		tscores = data$Threshold
		signf   = data$Score >= data$Threshold
	}else if (str_detect (f, "Plink")) {
		data    = data [1:nROWS,]
		tool    = "Plink"
		snps    = data$SNP
		pVal    = data$BONF
		chrom   = data$CHR
		pos	    = NA
		pscores = round (-log10 (data$BONF), 4)
		tscores = THRESHOLD
		signf   = pscores >= tscores
		
	}else if (str_detect (f, "Tassel")) {
		data = data %>% top_n (-1*nROWS, p)
		
		tool    = "Tassel"
		snps    = data$Marker
		#pVal    = nMARKERS*min (data$p, data$add_p, data$dom_p, na.rm=T)
		pVal    = data %>% rowwise %>% mutate (minP=min(data$p, data$add_p, data$dom_p, na.rm=T)) %>% select (minP)
		chrom   = data$Chr
		pos		= data$Pos
		pscores = round (-log10 (pVal),4)
		tscores = THRESHOLD
		signf   = pscores >= tscores
	}

	df = data.frame (TOOL=tool, MODEL=model, CHR=chrom, POS=pos, SNP=snps, 
					P = pVal, SCORE=pscores, THRESHOLD=tscores, Signf=signf )
	df = df %>% distinct (SNP, .keep_all=T)
	summaryTable = rbind (summaryTable, df)
	summTableSorted = summaryTable %>% add_count (SNP, sort=T, name="N1") %>% arrange (desc(N1))
	summTableSorted = summTableSorted %>% add_count (SNP, Signf, sort=T, name="Ns") %>% arrange (desc(Ns))
	write.table (file="summary-gwas.tbl", summaryTable, row.names=F,quote=F, sep="\t")
	write.table (file="summary-gwas-sorted.tbl", summTableSorted, row.names=F,quote=F, sep="\t")

}

# Create Venn diagram of common markers
createVennDiagramsMarkers (summaryTable)
