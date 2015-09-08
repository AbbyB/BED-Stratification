
# Program to create a chart of the confidence intervals for False Positives, False Negatives, and Not Assessed
#  and then graph the 3 rates for each chart
# Input: All of the filename extension '.csv' files in the current directory, 
#         which should be csv charts which include a column 'VCF' whose values when sorted alphabetically
#         the first four must correspond to the following categories:
#          False Negatives, False Positives, Not Assessed, and True Positives
#         This program treats the first column that only contains 0s and 1s as the beginning of the BED files
# 
# Output: 1 CSV file of the confidence intervals, 
#          moving left to right: leftmost column ("VCF") has the chart file name without the extension next to the VCF name, 
#          then BED file, estimate rate, lower bound, upper bound
# 
# Run the comment-ed 'install.packages' lines the first time this program is run only
# When running this program with different BED or VCF files, take note of comments with '**'

# install.packages("data.table")  # only run this the first time to install
# install.packages("plyr")
# install.packages("Hmisc")
# install.packages("ggplot2")
library(data.table)
library(plyr)
library(Hmisc)
library(ggplot2)

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

# List of the VCF file names
VCFNamesSorted <- array(dim=c(length(files),4))
for (z in 1:length(files)) {
  VCFNamesSorted[z,] <- sort(unique(charts[[z]][,VCF]))
}
# 1. FN, 2. FP, 3. NA, 4. TP

# **If you want to only process selected BED files use:
# **Select which columns
# bedFiles <- c('VCF','gc15','gc15to20','gc20to25','gc25to30','gc30to55','gc55to60','gc60to65','gc65to70','gc70to75','gc75to80','gc80to85','gc85',
#               'compositional','hg19SelfChainSplit','cpg.islands',
#               'SRHomopolymer.3to5','SRHomopolymer.6to10','SRDiTR.11to50','SRDiTR.51to200',
#               'SRTriTR.11to50','SRTriTR.51to200','SRQuadTR.11to50','SRQuadTR.51to200',
#               'lt7.lt51bp','lt7.51to200bp','lt7.gt200bp','gt6.lt51bp','gt6.51to200bp','gt6.gt200bp', 'normRegion')
# 
# for (z in 1:length(charts)) {
#   charts[[z]] <- charts[[z]][,bedFiles,with=FALSE]
# }

# The number of columns
chartEnd <- ncol(charts[[1]])
# The column number of the first column of the chart
#  of BED files and not variant info,
#  these values should be the same for all charts
# Finds the first column that has only binary values
x<-1
while(!all(unique(charts[[1]][[x]]) %in% c(0,1))) {
  if(x==chartEnd){stop("Please set the value for firstBED")}
  x<-x+1}
firstBED <- x
# **Set firstBED to first BED file column if the error message "Please set the value for firstBED" appears:
#    firstBED <- 

# Compute
# Create the first row of the binomial confidence interval
#  creates a new row that is:
#  Chart name, BED name, estimate point, lower bound, upper bound
confInt <- lapply(1:length(charts),function(z){
  rbindlist(list(as.data.table(t(as.data.table(lapply(firstBED:chartEnd, function(x){cbind(paste(chartName[z],"False Negatives"),names(charts[[z]])[x],binconf(sum(charts[[z]][VCF==VCFNamesSorted[z,1],x,with=FALSE]),sum(charts[[z]][VCF==VCFNamesSorted[z,1] | VCF==VCFNamesSorted[z,4],x,with=FALSE])))})))),
                 as.data.table(t(as.data.table(lapply(firstBED:chartEnd, function(x){cbind(paste(chartName[z],"False Positives"),names(charts[[z]])[x],binconf(sum(charts[[z]][VCF==VCFNamesSorted[z,2],x,with=FALSE]),sum(charts[[z]][VCF==VCFNamesSorted[z,2] | VCF==VCFNamesSorted[z,4],x,with=FALSE])))})))),
                 as.data.table(t(as.data.table(lapply(firstBED:chartEnd, function(x){cbind(paste(chartName[z],"Not Assessed"),names(charts[[z]])[x],binconf(sum(charts[[z]][VCF==VCFNamesSorted[z,3],x,with=FALSE]),sum(charts[[z]][VCF==VCFNamesSorted[z,3] | VCF==VCFNamesSorted[z,4],x,with=FALSE])))}))))))
})
confInt <- rbindlist(confInt)
# Rename the first 2 columns
setnames(confInt,c("VCF","BED.File","PointEst","Lower","Upper"))

# Create a data frame of the total variants from each VCF
totals <- ldply(1:length(charts), .fun=function(z){
  rbind(nrow(charts[[z]][VCF==VCFNamesSorted[z,1],]),
        nrow(charts[[z]][VCF==VCFNamesSorted[z,2],]),
        nrow(charts[[z]][VCF==VCFNamesSorted[z,3],]),
        nrow(charts[[z]][VCF==VCFNamesSorted[z,4],]))
})

# Save the chart of confidence intervals
write.csv(file="Confidence_Intervals.csv", x=confInt)


#-----------------------------------------------------------------------------------


# Graphing

# Set the cell type to numeric
confInt[, PointEst := as.numeric(confInt[,PointEst])]
confInt[, Lower := as.numeric(confInt[,Lower])]
confInt[, Upper := as.numeric(confInt[,Upper])]

# Sets indefinite cells to 0
#  the graph will show no bar or error bars
confInt[is.nan(confInt[,PointEst]) | is.na(confInt[,PointEst]), PointEst := 0]

# List of the VCF file names
VCFNames <- unique(confInt[,VCF])
# List of the bed file names
bedFiles <- confInt[VCF == VCFNames[1],BED.File]

# Create a vector of the average rate in the whole genome
VCFNum <- c(1:(length(charts)*4))
lines <- lapply(VCFNum[! VCFNum %in% seq(0,length(charts)*4,4)], function(z){
  totals[z, ]/ sum(totals[z,],totals[ceil(z/4)*4,])
})

# Each bar is grouped by the BED File type, legend present with the BED file types
# Create a vector that assigns each BED file to a group
# **Change these values
bedGroups <- c()
bedGroups[1:12] <- "GC Content"
bedGroups[13:18] <- "Unimask"
bedGroups[19:21] <- "Miscellaneous"
bedGroups[22:33] <- "Simple Repeats"
bedGroups[34:41] <- "Low Complexity"
bedGroups[42] <- "Ideal Region"

# Save each VCF graph to a variable
graphs <- lapply(1:length(VCFNames),function(y){
  # Creates simple bar graph for 1 VCF, assigns the data to x, y, and fill
  ggplot(data=confInt[VCF == VCFNames[y],], aes(x=factor(bedFiles, levels=bedFiles), y=PointEst, fill=bedGroups)) + 
    # Create a continous y axis
    geom_bar(stat="identity") + 
    # Add legend for BED file groupings
    scale_fill_discrete(name="BED File Groupings", breaks=unique(bedGroups)) + 
    # Set all text color to black
    theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
          # Make the x axis horizontal
          axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1, size=12),
          # Adjust text size and face
          axis.text.y  = element_text(size=12),
          legend.text = element_text(size=12),
          axis.title.x = element_text(size=15, face="bold"),
          axis.title.y = element_text(size=15, face="bold"),
          legend.title = element_text(size=15), 
          plot.title = element_text(size=18, face="bold")) + 
    # Add error bars
    geom_errorbar(aes(ymax = Upper, ymin = Lower), position=position_dodge(width=0.9), width=0.25) + 
    # Scale the y axis
    ylim(0, 1) + 
    # Add a title
    ggtitle(paste(VCFNames[y],"Rate")) + 
    # Add horizontal line for typical genome regions
    geom_hline(yintercept = lines[[y]], colour="#555555", linetype="dashed") + 
    # Label the axes
    labs(x = "BED File", y = "Rate")
})

# Print all 6 graphs separately
graphs[[1]]
graphs[[2]]
graphs[[3]]
graphs[[4]]
graphs[[5]]
graphs[[6]]

#----------------------------

# Show all graphs in one frame

# Reformat the confidence interval data frame
IS <- c(rep("Indels",nrow(confInt)/2),rep("SNPs",nrow(confInt)/2))
# Create Indels or SNPs column
confInt2 <- cbind(VCFp1 = IS, confInt)
setnames(confInt2,"VCF","VCFp2")
# Create false negatives, false positives, or not assessed column
VCFp2Descriptor <- c("False Negatives","False Positives","Not Assessed","False Negatives","False Positives","Not Assessed")
for (y in 1:length(VCFNames)) {
  confInt2[VCFp2 == VCFNames[y],VCFp2 := VCFp2Descriptor[y]]
}
# Reformat the average rate data frame
linesChart <- data.frame(VCFp1 = IS, VCFp2 = VCFp2Descriptor, V1 = as.data.frame(t(as.data.frame(lines))))

# Assign x, y, and fill
allGraph <- ggplot(confInt2, aes(x=BED.File, y=PointEst, fill=c(rep(bedGroups,6)))) + 
  # Create a continous y axis
  geom_bar(stat="identity") + 
  # Keep the x-axis in the correct order
  scale_x_discrete(limits=bedFiles) +
  # Split into six graphs
  facet_grid(VCFp1 ~ VCFp2) +
  # Add legend for BED file groupings
  scale_fill_discrete(name="BED File Groupings", breaks=unique(bedGroups)) + 
  # Make the x axis horizontal
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1, size=12, colour = "black"),
        # Adjust text size, color, and face
        axis.text.y  = element_text(size=12, colour = "black"),
        legend.text = element_text(size=12, colour = "black"),
        strip.text.x = element_text(size=15, face="bold"),
        strip.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=13, face="bold"),
        axis.title.y = element_text(size=13, face="bold"),
        legend.title = element_text(size=13, colour = "black")) + 
  # Add error bars
  geom_errorbar(aes(ymax = Upper, ymin = Lower), position=position_dodge(width=0.9), width=0.25) + 
  # Scale the y axis
  ylim(0, 1) + 
  # Add horizontal line for typical genome regions
  geom_hline(data = linesChart, aes(yintercept = V1), colour="#555555", linetype="dashed") + 
  # Label the axes
  labs(x = "BED File", y = "Rate")

# Plot all graphs in one frame
allGraph
