
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

# install.packages("Hmisc")  # only run this the first time to install
# install.packages("data.table")
# install.packages("ggplot2")
library(Hmisc)
library(data.table)
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
  rbindlist(list(as.data.table(t(as.data.table(lapply(firstBED:chartEnd, function(x){cbind(paste(chartName[z],"False Negatives"),names(charts[[z]])[x],binconf(sum(charts[[z]][VCF==VCFNamesSorted[1],x,with=FALSE]),sum(charts[[z]][VCF==VCFNamesSorted[1] | VCF==VCFNamesSorted[4],x,with=FALSE])))})))),
                 as.data.table(t(as.data.table(lapply(firstBED:chartEnd, function(x){cbind(paste(chartName[z],"False Positives"),names(charts[[z]])[x],binconf(sum(charts[[z]][VCF==VCFNamesSorted[2],x,with=FALSE]),sum(charts[[z]][VCF==VCFNamesSorted[2] | VCF==VCFNamesSorted[4],x,with=FALSE])))})))),
                 as.data.table(t(as.data.table(lapply(firstBED:chartEnd, function(x){cbind(paste(chartName[z],"Not Assessed"),names(charts[[z]])[x],binconf(sum(charts[[z]][VCF==VCFNamesSorted[3],x,with=FALSE]),sum(charts[[z]][VCF==VCFNamesSorted[3] | VCF==VCFNamesSorted[4],x,with=FALSE])))}))))))
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

# Run first - function to allow multiplot to work
# Source: R cookbook
#         http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Set the cell type to numeric
confInt[, PointEst := as.numeric(confInt[,PointEst])]
confInt[, Lower := as.numeric(confInt[,Lower])]
confInt[, Upper := as.numeric(confInt[,Upper])]

# Sets indefinite cells to 0
#  the graph will show no bar or error bars
confInt[is.nan(confInt[,PointEst]), PointEst := 0]
# **If reading in the confidence intervals file,
#    use the line below to remove blank cells:
#     confInt[is.na(confInt[,PointEst]), PointEst := 0]

# List of the VCF file names
VCFNames <- unique(confInt[,VCF])
# List of the bed file names
bedFiles <- confInt[VCF == VCFNames[1],BED.File]

# Create a vector of the average rate in the whole genome
sequence <- c(1:(length(charts)*4))
lines <- lapply(sequence[! sequence %in% seq(0,length(charts)*4,4)], function(z){
  totals[z, ]/ sum(totals[z,],totals[ceil(z/4)*4,])
})

#GRAPHING WAY 1 -------------

# Each bar is colored differently, no legend
# Save each VCF graph to a variable
graphs <- lapply(1:length(VCFNames),function(y){
  # Creates simple bar graph for 1 VCF, assigns the data to x and y
  ggplot(data=confInt[VCF == VCFNames[y],], aes(x=factor(bedFiles, levels=bedFiles), y=PointEst, fill=factor(bedFiles, levels=bedFiles))) + 
    # Create a continous y axis
    geom_bar(stat="identity") + 
    # Remove legend
    guides(fill=FALSE) + 
    # Make the x axis horizontal
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1)) + 
    # Add error bars
    geom_errorbar(aes(ymax = Upper, ymin = Lower), position=position_dodge(width=0.9), width=0.25) + 
    # Scale the y axis
    ylim(0, 1) + 
    # Add a title
    ggtitle(paste(VCFNames[y],"Rate")) + 
    # Add horizontal line for typical genome regions
    geom_hline(yintercept = lines[y], colour="#555555", linetype="dashed") + 
    # Label the axes
    labs(x = "BED File", y = "Rate")
})

#----------------------------

#GRAPHING WAY 2 -------------

# Each bar is grouped by the BED File type, legend present
# Create a vector that assigns each BED file to a group
# **Change these values
bedGroups <- c()
bedGroups[1:12] <- "GC Content"
bedGroups[13:20] <- "Imperfect Repeats"
bedGroups[21:26] <- "Unimask"
bedGroups[27:29] <- "Miscellaneous"
bedGroups[30:41] <- "Simple Repeats"
bedGroups[42] <- "Ideal Region"
# Create a vector for each group once in the order they appear along the x axis
totalGroups <- unique(bedGroups)

# Save each VCF graph to a variable
graphs <- lapply(1:length(VCFNames),function(y){
  # Creates simple bar graph for 1 VCF, assigns the data to x and y
  ggplot(data=confInt[VCF == VCFNames[y],], aes(x=factor(bedFiles, levels=bedFiles), y=PointEst, fill=bedGroups)) + 
    # Create a continous y axis
    geom_bar(stat="identity") + 
    # Add legend for BED file groupings
    scale_fill_discrete(name="BED File Groupings", breaks=totalGroups) + 
    # Make the x axis horizontal
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1)) + 
    # Add error bars
    geom_errorbar(aes(ymax = Upper, ymin = Lower), position=position_dodge(width=0.9), width=0.25) + 
    # Scale the y axis
    ylim(0, 1) + 
    # Add a title
    ggtitle(paste(VCFNames[y],"Rate")) + 
    # Add horizontal line for typical genome regions
    geom_hline(yintercept = lines[y], colour="#555555", linetype="dashed") + 
    # Label the axes
    labs(x = "BED File", y = "Rate")
})

#----------------------------

# Print all 6 graphs separately
graphs[[1]]
graphs[[2]]
graphs[[3]]
graphs[[4]]
graphs[[5]]
graphs[[6]]

# Plot all graphs on one page (click zoom to see)
multiplot(plotlist = graphs, cols=round(sqrt(length(graphs))))
