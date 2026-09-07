library(optparse)

option_list <- list(
  make_option(c("--WIG"), type = "character", help = "Path to tumor WIG file. Required."),
  make_option(c("--NORMWIG"), type = "character", default=NULL, help = "Path to normal WIG file. Default: [%default]"),
  make_option(c("--gcWig"), type = "character", help = "Path to GC-content WIG file; Required"),
  make_option(c("--mapWig"), type = "character", default=NULL, help = "Path to mappability score WIG file. Default: [%default]"),
  make_option(c("--repTimeWig"), type = "character", default=NULL, help ="Path to replication timing WIG file. Default: [%default]"),
  make_option(c("--normalPanel"), type="character", default=NULL, help="Median corrected depth from panel of normals. Default: [%default]"),
  make_option(c("--sex"), type = "character", default = NULL, help = "User specified gender: male or female [Default: %default]"),
  make_option(c("--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions. Default: [%default]"),
  make_option(c("--id"), type = "character", default="test", help = "Patient ID. Default: [%default]"),
  make_option(c("--centromere"), type="character", default=NULL, help = "File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package. Default: [%default]"),
  make_option(c("--minMapScore"), type = "numeric", default=0.9, help="Include bins with a minimum mappability score of this value. Default: [%default]."),
  make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove. Default: [%default]"),
  make_option(c("--normal"), type="character", default="0.5", help = "Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [%default]"),
  make_option(c("--normal.init"), type="character", default="c(0.5, 0.5)", help = "Specific initialization of normal contamination for multiple samples. Default: [%default]"),
  make_option(c("--scStates"), type="character", default="NULL", help = "Subclonal states to consider. Default: [%default]"),
  make_option(c("--scPenalty"), type="numeric", default=0.1, help = "Penalty for subclonal state transitions. 0.1 penalizes subclonal states by ~10%. Default: [%default]"),
  make_option(c("--normal2IgnoreSC"), type="numeric", default=1.0, help="Ignore subclonal analysis when normal proportion is >= this value. Default: [%default]"),
  make_option(c("--coverage"), type="numeric", default=NULL, help = "PICARD sequencing coverage. Default: [%default]"),
  make_option(c("--likModel"), type="character", default="t", help="Likelihood model to use. \"t\" or \"gaussian\". Use \"gaussian\" for faster runtimes. Default: [%default]"),
  make_option(c("--lambda"), type="character", default="NULL", help="Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data. Default: [%default]"),
  make_option(c("--lambdaScaleHyperParam"), type="numeric", default=3, help="Hyperparameter (scale) for Gamma prior on Student's-t precision. Default: [%default]"),
  make_option(c("--ploidy"), type="character", default="2", help = "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [%default]"),
  make_option(c("--maxCN"), type="numeric", default=7, help = "Total clonal CN states. Default: [%default]"),
  make_option(c("--estimateNormal"), type="logical", default=TRUE, help = "Estimate normal. Default: [%default]"),
  make_option(c("--estimateScPrevalence"), type="logical", default=TRUE, help = "Estimate subclonal prevalence. Default: [%default]"),
  make_option(c("--estimatePloidy"), type="logical", default=TRUE, help = "Estimate tumour ploidy. Default: [%default]"),
  make_option(c("--maxFracCNASubclone"), type="numeric", default=0.7, help="Exclude solutions with fraction of subclonal events greater than this value. Default: [%default]"),
  make_option(c("--maxFracGenomeSubclone"), type="numeric", default=0.5, help="Exclude solutions with subclonal genome fraction greater than this value. Default: [%default]"),
  make_option(c("--minSegmentBins"), type="numeric", default=50, help="Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."),
  make_option(c("--altFracThreshold"), type="numeric", default=0.05, help="Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [%default]"),
  make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases. Default: [%default]"),
  make_option(c("--chrTrain"), type="character", default="c(1:22)", help = "Specify chromosomes to estimate params. Default: [%default]"),
  make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze. Default: [%default]"),
  make_option(c("--genomeBuild"), type="character", default="hg19", help="Geome build. Default: [%default]"),
  make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
  make_option(c("--normalizeMaleX"), type="logical", default=TRUE, help = "If male, then normalize chrX by median. Default: [%default]"),
  make_option(c("--fracReadsInChrYForMale"), type="numeric", default=0.001, help = "Threshold for fraction of reads in chrY to assign as male. Default: [%default]"),
  make_option(c("--includeHOMD"), type="logical", default=FALSE, help="If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: [%default]"),
  make_option(c("--txnE"), type="numeric", default=0.9999999, help = "Self-transition probability. Increase to decrease number of segments. Default: [%default]"),
  make_option(c("--txnStrength"), type="numeric", default=1e7, help = "Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [%default]"),
  make_option(c("--multSampleTxnStrength"), type="numeric", default=1, help="Strength of same state transition between multiple samples. Default: [%default]"),
  make_option(c("--plotFileType"), type="character", default="pdf", help = "File format for output plots. Default: [%default]"),
  make_option(c("--plotYLim"), type="character", default="c(-2,2)", help = "ylim to use for chromosome plots. Default: [%default]"),
  make_option(c("--outDir"), type="character", default="./", help = "Output Directory. Default: [%default]"),
  make_option(c("--libdir"), type = "character", default=NULL, help = "Script library path. Usually exclude this argument unless custom modifications have been made to the ichorCNA R package code and the user would like to source those R files. Default: [%default]"),
  make_option(c("--cores"), type="numeric", default = 1, help = "Number of cores to use for EM. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)

if (!is.null(opt$libdir) && opt$libdir != "None"){
	source(paste0(opt$libdir,"/R/utils.R"))
	source(paste0(opt$libdir,"/R/segmentation.R"))
	source(paste0(opt$libdir,"/R/EM.R"))
	source(paste0(opt$libdir,"/R/output.R"))
	source(paste0(opt$libdir,"/R/plotting.R"))
	source(paste0(opt$libdir,"/R/runIchorCNA.R"))
} else {
	library(ichorCNA)
}

run_ichorCNA(tumor_wig = opt$WIG, normal_wig = opt$NORMWIG, gcWig = opt$gcWig, mapWig = opt$mapWig,
             repTimeWig = opt$repTimeWig, normal_panel = opt$normalPanel, sex = opt$sex,
             exons.bed = opt$exons.bed, id = opt$id, centromere = opt$centromere,
             minMapScore = opt$minMapScore, flankLength = opt$rmCentromereFlankLength,
             normal = opt$normal, normal.init = opt$normal.init, scStates = opt$scStates,
             scPenalty = opt$scPenalty, normal2IgnoreSC = opt$normal2IgnoreSC, coverage = opt$coverage,
             likModel = opt$likModel, lambda = opt$lambda, lambdaScaleHyperParam = opt$lambdaScaleHyperParam,
             ploidy = opt$ploidy, maxCN = opt$maxCN, estimateNormal = opt$estimateNormal,
             estimateScPrevalence = opt$estimateScPrevalence, estimatePloidy = opt$estimatePloidy,
             maxFracCNASubclone = opt$maxFracCNASubclone, maxFracGenomeSubclone = opt$maxFracGenomeSubclone,
             minSegmentBins = opt$minSegmentBins, altFracThreshold = opt$altFracThreshold,
             chrNormalize = opt$chrNormalize, chrTrain = opt$chrTrain, chrs = opt$chrs,
             genomeBuild = opt$genomeBuild, genomeStyle = opt$genomeStyle, normalizeMaleX = opt$normalizeMaleX,
             fracReadsInChrYForMale = opt$fracReadsInChrYForMale, includeHOMD = opt$includeHOMD,
             txnE = opt$txnE, txnStrength = opt$txnStrength, multSampleTxnStrength = opt$multSampleTxnStrength,
             plotFileType = opt$plotFileType, plotYLim = opt$plotYLim, outDir = opt$outDir, cores = opt$cores)
