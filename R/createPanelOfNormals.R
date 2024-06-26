#' Runs createPanelOfNormals for ichorCNA
#'
#' @param gcWig    GC Wig file for reference genome  
#' @param mapWig    Mappabiliy Wig file for reference genome  
#' @param repTimeWig  type = "character", default=NULL, help ="Path to replication timing WIG file.    
#' @param filelist    List of of wig files.  
#' @param outfile    Output file.  
#' @param centromere    File containing Centromere locations  
#' @param flankLength    Length of region flanking centromere to remove.    
#' @param chrs    Specify chromosomes to analyze.  
#' @param genomeStyle    NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. 
#' @param genomeBuild    Geome build.    
#' @param chrNormalize    Specify chromosomes to normalize GC/mappability biases  
#' @param minMapScore    Include bins with a minimum mappability score of this value. 
#' @param maleChrXLogRThres    ChrX Log ratio threshold to confirm as male gender.  
#' @param fracReadsInChrYForMale    Threshold for fraction of reads in chrY to assign as male.    
#' @param exons.bed    Path to bed file containing exon regions.  
#' @param method    Median or Mean.  
#' @param ylim    Y-limits for plotting of mean/median log ratios  
#' @param plotChrPanels    Plot PoN values.
#' @export
createPanelOfNormals <- function(gcWig, mapWig, repTimeWig = NULL, filelist, outfile, centromere, flankLength = 1e5,
                                 chrs = "c(1:22,\"X\")", genomeStyle = "NCBI", genomeBuild = "hg19", 
                                 chrNormalize = "c(1:22)", minMapScore = 0.0, maleChrXLogRThres = -0.80, 
                                 fracReadsInChrYForMale = 0.001, exons.bed = NULL, method = "median", ylim = "c(-2,2)", 
                                 plotChrPanels = FALSE) {
  require(HMMcopy)
  require(GenomicRanges)
  require(ggplot2)

options(stringsAsFactors=FALSE, scipen=0)
options(bitmapType='cairo')

ylim <- eval(parse(text =ylim))
chrs <- as.character(eval(parse(text =chrs)))
chrNormalize <- as.character(eval(parse(text=chrNormalize))); 
seqlevelsStyle(chrs) <- genomeStyle
seqlevelsStyle(chrNormalize) <- genomeStyle


if (!is.null(centromere)){
	centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
}
seqlevelsStyle(centromere$Chr) <- genomeStyle

files <- read.delim(filelist, header = FALSE, stringsAsFactors=FALSE, sep ="\t")[, 1]



## LOAD GC/MAP/REPTIME WIG FILES ###
message("Reading GC and mappability files")
gc <- wigToGRanges(gcWig)
if (is.null(gc)){
    stop("GC wig file not provided but is required")
}
map <- wigToGRanges(mapWig)
if (is.null(map)){
  message("No mappability wig file input, excluding from correction")
}
repTime <- wigToGRanges(repTimeWig)
if (is.null(repTime)){
  message("No replication timing wig file input, excluding from correction")
}else{
  if (mean(repTime$value, na.rm = TRUE) > 1){
    repTime$value <- repTime$value / 100 ## values in [0,1] - for LNCaP_repTime_10kb_hg38.txt
  }
}

normalGR <- NULL
info <- data.frame(Sample = character(), sex = character(), chrYcov = numeric(),
	chrXMedian = numeric(), normChrX = numeric())
for (i in 1:length(files)){
	### LOAD NORMAL FILES ###
	sid <- gsub(".wig","",basename(files[i]))
	message("Loading normal file:", files[i])
	normal_reads <- wigToGRanges(files[i])
		
	## FILTER BY EXONS IF PROVIDED ##
	## add gc and map to RangedData object ##
	if (!is.null(exons.bed)){
		targetedSequences <- read.delim(exons.bed, header=F, sep="\t", skip=86)
	}else{
		targetedSequences <- NULL
	}
	# normal_counts <- loadReadCountsFromWig(normal_reads, chrs=chrs, gc=gc, map=map, 
	# 				centromere=centromere, targetedSequences=targetedSequences)
	normal_counts <- loadReadCountsFromWig(normal_reads, chrs = chrs, gc = gc, map = map, 
									   repTime = repTime, centromere = centromere, 
									   flankLength = flankLength, targetedSequences = targetedSequences, 
									   chrXMedianForMale = maleChrXLogRThres, 
									   genomeStyle = genomeStyle, fracReadsInChrYForMale = fracReadsInChrYForMale,
                                       chrNormalize = chrNormalize, mapScoreThres = minMapScore)
	gender <- normal_counts$gender
	### CORRECT TUMOUR DATA FOR GC CONTENT AND MAPPABILITY BIASES ###
	message("Correcting ", sid, " sex: ", gender$gender, 
		" chrYcov: ", gender$chrYCovRatio, 
		" chrXMedian: ", gender$chrXMedian)
	if (is.null(normalGR)){
		normalGR <- normal_counts$counts
		values(normalGR) <- values(normalGR)$copy
		colnames(values(normalGR)) <- sid
	}else{
		values(normalGR)[[sid]] <- normal_counts$counts$copy
	}
	
	## Normalize chrX -- if male sex ##
	chrXStr <- grep("X", chrs, value = TRUE)
	chrXInd <- as.character(seqnames(normalGR)) == chrXStr		
	if (gender$gender == "male"){
		chrXMedian <- gender$chrXMedian		
		values(normalGR)[[sid]][chrXInd] <- values(normalGR)[[sid]][chrXInd] - chrXMedian
	}
	
	info <- rbind(info, data.frame(Sample = sid, sex = gender$gender, 
		chrYcov = gender$chrYCovRatio, chrXMedian = gender$chrXMedian,
		normChrX = median(values(normalGR)[[sid]][chrXInd], na.rm=T)))


}
print (info)

mat <- values(normalGR)
if (method == "median"){
   medianVal <- apply(mat, 1, median, na.rm = TRUE)
}else if (method == "mean"){
   medianVal <- apply(mat, 1, mean, na.rm = TRUE)
}else{
  stop("method is not specified as median or mean.")
}
values(normalGR)[["Median"]] <- medianVal

####################################
## Save the output ##
####################################

## save GR object ##
saveRDS(normalGR, file = paste0(outfile, "_", method, ".rds"))

## output to text file ##
write.table(as.data.frame(normalGR[,"Median"]), file=paste0(outfile, "_", method, ".txt"), col.names=TRUE, row.names=F, quote=F, sep="\t")

####################################
## Plots ##
####################################
## plot the median profile of the PoN
normalGR <- filterEmptyChr(normalGR)
chrLvls <- levels(seqnames(normalGR))
chrsToUse <- as.vector(seqnames(normalGR))
starts <- start(ranges(normalGR))
ends <- end(ranges(normalGR))
midPts <- starts + ((ends - starts) / 2) + 1
# get chr info from USCS using GenomeInfoDB #
seqInfo <- Seqinfo(genome = genomeBuild)
seqlevelsStyle(seqInfo) <- genomeStyle

#print(seqInfo)

seqInfo <- keepSeqlevels(seqInfo, chrs)
if (plotChrPanels){
	coord <- NULL
  	coord$posns <- midPts
}else{
	coord <- getGenomeWidePositions(chrsToUse, midPts, seqInfo)
}
mat <- data.frame(Chrs = chrsToUse, Positions = coord$posns, Median = normalGR$Median)
mat$Chrs <- factor(mat$Chrs, levels = chrLvls)

gp <- ggplot(mat, aes(x = Positions, y = Median)) +
	geom_point(alpha = 0.75, size = 1, shape = 19) +# + geom_line() 
	ylim(ylim) +
	theme(panel.background = element_blank(),
				panel.border = element_rect(colour = "black", fill = NA),
				panel.grid.major = element_blank(),#element_line(colour = "grey75", linetype="dashed"), 
				panel.grid.minor = element_blank(),#element_blank(), axis.title=element_text(size = 14, face = "bold"),
				axis.text = element_text(size = 12))


if (plotChrPanels){ # additional chr panel attributes #
	gp <- gp + facet_wrap( ~ Chrs, ncol = 4, scales = "free_x") + 
	  scale_x_continuous(breaks = NULL) + expand_limits(x = 0)
}else{
	numLines <- length(coord$chrBkpt)
	mid <- (coord$chrBkpt[1:(numLines - 1)] + coord$chrBkpt[2:numLines]) / 2
	gp <- gp + scale_x_continuous(name = "Chromosome", labels = chrs, breaks = mid, limits = c(coord$chrBkpt[1], tail(coord$chrBkpt ,1))) +
	  geom_vline(xintercept = coord$chrBkpt, linetype = "dotted")
}
outplotfile <- paste0(outfile, "_", method, ".png")
ggsave(gp, file = outplotfile, width = 16, height = 4)	



}

