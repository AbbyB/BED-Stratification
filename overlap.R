
# Program to analyze overlap in BED groups in False Positives, False Negatives, and Not Assessed
# Input: All of the filename extension '.csv' files in the current directory, 
#         which should be csv charts which include a column 'VCF' with the following values:
#          True Positives, False Positives, False Negatives, and Not Assessed
#         and the columns of the chart should include:
#          gc15, gc15to20, gc20to25, gc25to30, gc30to55, gc55to60, gc60to65, gc65to70, gc70to75, gc75to80, gc80to85, gc85, 
#          lt7.lt51bp, lt7.51to200bp, lt7.gt200bp, gt6.lt51bp, gt6.51to200bp, gt6.gt200bp, 
#          compositional, cpg.islands,
#          SRHomopolymer.3to5, SRHomopolymer.6to10, SRHomopolymer.gt10, SRDiTR.11to50, SRDiTR.51to200, SRDiTR.gt200, 
#          SRTriTR.11to50, SRTriTR.51to200, SRTriTR.gt200, SRQuadTR.11to50, SRQuadTR.51to200, SRQuadTR.gt200
# 
# Output: 1 CSV file of the overlap, 
#          moving left to right: leftmost column ("Name") has the chart file name without the extension next to 
#          the descriptor of the overlap (ainb where overlap is the counts in both sand notOverlap is those in b not in a); 
#          the counts for TP, FP, FN, and NA; the total counts; rate of overlap for FP, FN, and NA; 
#          the counts for TP, FP, FN, and NA not overlapping; the total counts for not overlapping; rates of not overlapping for FP, FN, and NA 
#          ratios of overlap to not overlap for TP, FP, FN, and NA; total ratio of overlap to not overlap; odds ratios by FP, FN, and NA
# 
# Run the commented 'install.packages' lines the first time this program is run only

#install.packages("data.table")   # only run this the first time to install
library(data.table)

# Input
# Finds all the files in the current directory with the file extension '.csv'
files <- list.files(pattern = "\\.csv$")
# **If you want to remove files from this list use:
#    files <- files[grep("(word to search for to remove file name)", files, invert = TRUE)]
#   Or if you want to select certain files use:
#    files <- files[grep("(keyword of file names)", files)]

# Reads all files found in the directory
charts <- list()
for (z in 1:length(files)) {
  charts[[z]] <- fread(files[[z]], header=TRUE)
}

# Labels the chart files using the file name without the '.csv' extension
chartName <- gsub('.csv','',files)
# **If you want different labels, change to:
#    chartName <- c("(label 1)", "(label 2)", "(etc. for the total number of files)")

# List of the VCF names
VCFNames <- unique(charts[[1]][,VCF])
# Set data table to sort by VCF
lapply(1:length(charts),function(z){setkey(charts[[z]],VCF)})

# Sum groups of columns by row
# Save the column names to a variable
SR <- names(charts[[1]])[grepl("SR", names(charts[[1]]))]
LC <- c("lt7.lt51bp","lt7.51to200bp","lt7.gt200bp","gt6.lt51bp","gt6.51to200bp","gt6.gt200bp")
# Loop through the charts and row sum the relevant columns
lapply(1:length(charts),function(z){
  # Low Complexity
  charts[[z]][,LC := rowSums(charts[[z]][,LC,with = FALSE])]
  # Simple Repeats
  charts[[z]][,SR := rowSums(charts[[z]][,SR,with = FALSE])]
  # Low Complexity and Simple Repeats
  charts[[z]][,LCSR := LC+ SR]
})

# Initialize the variables
overlap <- list()
names <- list()
bases <- list()
counts <- list()
# Loop through each of the charts to count each of the relevant overlaps
# Create a data table with row 1-4 FN, FP, NA, and TP respectively
for (z in 1:length(charts)){
  #GC content overlap with both sets of repeats
  counts[[z]] <- cbind(charts[[z]][LCSR > 0,list(sum(gc15),
                                                 sum(gc15to20),
                                                 sum(gc20to25),
                                                 sum(gc25to30),
                                                 sum(gc30to55),
                                                 sum(gc55to60),
                                                 sum(gc60to65),
                                                 sum(gc65to70),
                                                 sum(gc70to75),
                                                 sum(gc75to80),
                                                 sum(gc80to85),
                                                 sum(gc85)),by = VCF][,VCF := NULL],
                       
                       #Low complexity overlap with simple repeats
                       charts[[z]][SR > 0,list(sum(lt7.lt51bp),
                                               sum(lt7.51to200bp),
                                               sum(lt7.gt200bp),
                                               sum(gt6.lt51bp),
                                               sum(gt6.51to200bp),
                                               sum(gt6.gt200bp)),by = VCF][,VCF := NULL],
                       
                       #Compositional overlap with GC content, repeated twice
                       rep(charts[[z]][compositional == 1,list(sum(gc15),
                                                               sum(gc15to20),
                                                               sum(gc20to25),
                                                               sum(gc25to30),
                                                               sum(gc30to55),
                                                               sum(gc55to60),
                                                               sum(gc60to65),
                                                               sum(gc65to70),
                                                               sum(gc70to75),
                                                               sum(gc75to80),
                                                               sum(gc80to85),
                                                               sum(gc85)),by = VCF][,VCF := NULL],2),
                       
                       #Compositional overlap with low complexity
                       charts[[z]][compositional == 1,list(sum(lt7.lt51bp),
                                                           sum(lt7.51to200bp),
                                                           sum(lt7.gt200bp),
                                                           sum(gt6.lt51bp),
                                                           sum(gt6.51to200bp),
                                                           sum(gt6.gt200bp)),by = VCF][,VCF := NULL],
                       
                       #Compositional overlap with simple repeats
                       charts[[z]][compositional == 1,list(sum(SRHomopolymer.3to5),
                                                           sum(SRHomopolymer.6to10),
                                                           sum(SRHomopolymer.gt10),
                                                           sum(SRDiTR.11to50),
                                                           sum(SRDiTR.51to200),
                                                           sum(SRDiTR.gt200),
                                                           sum(SRTriTR.11to50),
                                                           sum(SRTriTR.51to200),
                                                           sum(SRTriTR.gt200),
                                                           sum(SRQuadTR.11to50),
                                                           sum(SRQuadTR.51to200),
                                                           sum(SRQuadTR.gt200)),by = VCF][,VCF := NULL],
                       
                       #Compositional overlap with both sets of repeats
                       charts[[z]][LCSR > 0,sum(compositional),by = VCF][,VCF := NULL],
                       
                       #GC content overlap with CPG Islands
                       charts[[z]][cpg.islands == 1,list(sum(gc15),
                                                         sum(gc15to20),
                                                         sum(gc20to25),
                                                         sum(gc25to30),
                                                         sum(gc30to55),
                                                         sum(gc55to60),
                                                         sum(gc60to65),
                                                         sum(gc65to70),
                                                         sum(gc70to75),
                                                         sum(gc75to80),
                                                         sum(gc80to85),
                                                         sum(gc85)),by = VCF][,VCF := NULL])
  
  # Names of overlaps
  names[[z]] <- c(paste(chartName[z],"LCSRin15"),paste(chartName[z],"LCSRin15to20"),paste(chartName[z],"LCSRin20to25"),paste(chartName[z],"LCSRin25to30"),paste(chartName[z],"LCSRin30to55"),paste(chartName[z],"LCSRin55to60"),paste(chartName[z],"LCSRin60to65"),paste(chartName[z],"LCSRin65to70"),paste(chartName[z],"LCSRin70to75"),paste(chartName[z],"LCSRin75to80"),paste(chartName[z],"LCSRin80to85"),paste(chartName[z],"LCSRin85"),
                  paste(chartName[z],"SRinlt7.lt51bp"),paste(chartName[z],"SRinlt7.51to200bp"),paste(chartName[z],"SRinlt7.gt200bp"),paste(chartName[z],"SRingt6.lt51bp"),paste(chartName[z],"SRingt6.51to200bp"),paste(chartName[z],"SRingt6.gt200bp"),
                  paste(chartName[z],"compin15"),paste(chartName[z],"compin15to20"),paste(chartName[z],"compin20to25"),paste(chartName[z],"compin25to30"),paste(chartName[z],"compin30to55"),paste(chartName[z],"compin55to60"),paste(chartName[z],"compin60to65"),paste(chartName[z],"compin65to70"),paste(chartName[z],"compin70to75"),paste(chartName[z],"compin75to80"),paste(chartName[z],"compin80to85"),paste(chartName[z],"compin85"),
                  paste(chartName[z],"15incomp"),paste(chartName[z],"15to20incomp"),paste(chartName[z],"20to25incomp"),paste(chartName[z],"25to30incomp"),paste(chartName[z],"30to55incomp"),paste(chartName[z],"55to60incomp"),paste(chartName[z],"60to65incomp"),paste(chartName[z],"65to70incomp"),paste(chartName[z],"70to75incomp"),paste(chartName[z],"75to80incomp"),paste(chartName[z],"80to85incomp"),paste(chartName[z],"85incomp"),
                  paste(chartName[z],"compinlt7.lt51bp"),paste(chartName[z],"compinlt7.51to200bp"),paste(chartName[z],"compinlt7.gt200bp"),paste(chartName[z],"compingt6.lt51bp"),paste(chartName[z],"compingt6.51to200bp"),paste(chartName[z],"compingt6.gt200bp"),
                  paste(chartName[z],"compinSRH3to5"),paste(chartName[z],"compinSRH6to10"),paste(chartName[z],"compinSRHgt10"),paste(chartName[z],"compinSRD11to50"),paste(chartName[z],"compinSRD51to200"),paste(chartName[z],"compinSRDgt200"),paste(chartName[z],"compinSRT11to50"),paste(chartName[z],"compinSRT51to200"),paste(chartName[z],"compinSRTgt200"),paste(chartName[z],"compinSRQ11to50"),paste(chartName[z],"compinSRQ51to200"),paste(chartName[z],"compinSRQgt200"),
                  paste(chartName[z],"LCSRincomp"),
                  paste(chartName[z],"islandsin15"),paste(chartName[z],"islandsin15to20"),paste(chartName[z],"islandsin20to25"),paste(chartName[z],"islandsin25to30"),paste(chartName[z],"islandsin30to55"),paste(chartName[z],"islandsin55to60"),paste(chartName[z],"islandsin60to65"),paste(chartName[z],"islandsin65to70"),paste(chartName[z],"islandsin70to75"),paste(chartName[z],"islandsin75to80"),paste(chartName[z],"islandsin80to85"),paste(chartName[z],"islandsin85"))
  
  # Total bases that potentially could have overlapped
  gcbases <- charts[[z]][,list(sum(gc15),
                               sum(gc15to20),
                               sum(gc20to25),
                               sum(gc25to30),
                               sum(gc30to55),
                               sum(gc55to60),
                               sum(gc60to65),
                               sum(gc65to70),
                               sum(gc70to75),
                               sum(gc75to80),
                               sum(gc80to85),
                               sum(gc85)),by = VCF][,VCF := NULL]
  irbases <- charts[[z]][,list(sum(lt7.lt51bp),
                               sum(lt7.51to200bp),
                               sum(lt7.gt200bp),
                               sum(gt6.lt51bp),
                               sum(gt6.51to200bp),
                               sum(gt6.gt200bp)),by = VCF][,VCF := NULL]
  srbases <- charts[[z]][,list(sum(SRHomopolymer.3to5),
                               sum(SRHomopolymer.6to10),
                               sum(SRHomopolymer.gt10),
                               sum(SRDiTR.11to50),
                               sum(SRDiTR.51to200),
                               sum(SRDiTR.gt200),
                               sum(SRTriTR.11to50),
                               sum(SRTriTR.51to200),
                               sum(SRTriTR.gt200),
                               sum(SRQuadTR.11to50),
                               sum(SRQuadTR.51to200),
                               sum(SRQuadTR.gt200)),by = VCF][,VCF := NULL]
  compbases <- charts[[z]][,sum(compositional),by = VCF][,VCF := NULL]
  # Order the totals to match the counts order
  bases[[z]] <- cbind(gcbases,irbases,gcbases,rep(compbases,12),irbases,srbases,compbases,gcbases)
  
  # Make the data frame to display results
  # Add the names and counts for TP, FP, FN, and NA
  overlap[[z]] <- cbind(Name = names[[z]],as.data.table(cbind(overlapTP = as.integer(counts[[z]][4,]), overlapFP = as.integer(counts[[z]][2,]), overlapFN = as.integer(counts[[z]][1,]), overlapNA = as.integer(counts[[z]][3,])))[, lapply(.SD, as.numeric)])
  # Add a total counts column
  overlap[[z]][, totalOverlap := (overlapTP + overlapFP + overlapFN + overlapNA)]
  # Add FP, FN, and NA rates of overlap
  overlap[[z]][, inFPTP := overlapFP/overlapTP][, inFNTP := overlapFN/overlapTP][, inNATP := overlapNA/overlapTP]
  # Add the counts for TP, FP, FN, and NA not overlapping
  overlap[[z]][, notOverlapTP := as.integer(bases[[z]][4,])-overlapTP][, notOverlapFP := as.integer(bases[[z]][2,])-overlapFP][, notOverlapFN := as.integer(bases[[z]][1,])-overlapFN][, notOverlapNA := as.integer(bases[[z]][3,])-overlapNA]
  # Add a total counts column for not overlapping
  overlap[[z]][, totalNotOverlap := (notOverlapTP + notOverlapFP + notOverlapFN + notOverlapNA)]
  # Add FP, FN, and NA rates of not overlapping
  overlap[[z]][, notInFPTP := notOverlapFP/notOverlapTP][, notInFNTP := notOverlapFN/notOverlapTP][, notInNATP := notOverlapNA/notOverlapTP]
  # Add ratios of overlap to not overlap by FP, FN, and NA
  overlap[[z]][, ratioTP := overlapTP/notOverlapTP][, ratioFP := overlapFP/notOverlapFP][, ratioFN := overlapFN/notOverlapFN][, ratioNA := overlapNA/notOverlapNA]
  # Add total ratio of overlap to not overlap
  overlap[[z]][, totalRatio := totalOverlap/totalNotOverlap]
  # Add odds ratios by FP, FN, and NA
  overlap[[z]][, oddsRatioFPTP := inFPTP/notInFPTP][, oddsRatioFNTP := inFNTP/notInFNTP][, oddsRatioNATP := inNATP/notInNATP]
}

# Save the file to a csv file
write.csv(file="overlap.csv", x=rbindlist(overlap))
