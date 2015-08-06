
# Program runs binomical conf interval tests and logistical regression for FN, FP, and NA rates for stratification bed files and graphs them
# Input: All of the filename extension '.csv' files in the current directory, 
#         which should be csv charts which include a column 'VCF' whose values when sorted alphabetically
#         the first four must correspond to the following categories:
#          False Negatives, False Positives, Not Assessed, and True Positives
# 
# Output: 1 CSV file that includes all the VCFs, 
#          moving left to right: leftmost column ("VCF") has the file name without the extension next to FN, FP, or NA, 
#          then estimate , lower 2.5% bound, upper 97.5% bound
# 
# Run the commentted 'install.packages' lines the first time this program is run only

# install.packages("Hmisc") #only run this the first time to install
# install.packages("data.table")
# install.packages("ggplot2")
library(Hmisc)
library(data.table)
library(ggplot2)

# Input
# Finds all the files in the current directory with the file extension '.csv'
files <- list.files(pattern = "\\.csv$")
# If you want to remove files from this list use:
#  files <- files[grep("(word to search for to remove file name)", files,invert = TRUE)]

# Reads all files found in the directory
charts <- list()
for (z in 1:length(files)) {
  charts[[z]] <- fread(files[[z]], header=TRUE)
}

# Labels the chart files using the file name without the '.csv' extension
chartName <- gsub('.csv','',files)
# **If you want different labels, change to:
#    chartName <- c((label 1), (label 2), (etc. for the total number of files))

# Number of regression models ran per chart
numReg <- 5

# List of the VCF file names
VCFNamesSorted <- sort(unique(charts[[1]][,VCF]))[1:4]
# 1. FN, 2. FP, 3. NA, 4. TP

# Initialize the variables
chartsFPTP <- list()
chartsFNTP <- list()
chartsNATP <- list()
logitReg <- list()
logitRegResults <- list()
# Loop through each chart
for (z in 1:length(charts)) {
  # Create columns for FN, FP, NA, and TP
  charts[[z]][,FN := as.numeric(charts[[z]][,VCF==VCFNamesSorted[1]])][,FP := as.numeric(charts[[z]][,VCF==VCFNamesSorted[2]])]
  charts[[z]][,NotA := as.numeric(charts[[z]][,VCF==VCFNamesSorted[3]])][,TP := as.numeric(charts[[z]][,VCF==VCFNamesSorted[4]])]
  # Create new data tables for each logistic regression
  chartsFPTP[[z]] <- charts[[z]][FP==1 | TP==1,]
  chartsFNTP[[z]] <- charts[[z]][FN==1 | TP==1,]
  chartsNATP[[z]] <- charts[[z]][NotA==1 | TP==1,]

  # Initialize y to be a counter for the results  
  y <- numReg*(z-1)+1
  
  #test  1
  logReg[[y]] <- glm(FP ~ gc15 + gc15to20 + gc20to25 + gc25to30 + gc30to55 + gc55to60 + gc60to65 + gc65to70 + gc70to75 + gc75to80 + gc80to85 + gc85 , data=chartsFPTP[[z]], family = binomial(logit))
  summary(logitReg[[y]])
  logRegResults[[y]] <- as.data.table(cbind(File = paste(chartName, "False Positives"),exp(cbind(Estimate = logReg[[y]]$coef, confint(logReg[[y]]))), keep.rownames=TRUE))
  setnames(logRegResults[[y]],c("BED","Estimate","Lower2.5","Upper97.5"))
  y <- y+1
  
  #test 2
  logReg[[y]] <- glm(FP ~ gc15 + gc15to20 + gc20to25 + gc25to30 + gc65to70 + gc70to75 + gc75to80 + gc80to85 + gc85 , data=chartsFPTP[[z]], family = binomial(logit))
  summary(logReg[[y]])
  logRegResults[[y]] <- as.data.table(cbind(File = paste(chartName, "False Positives"),exp(cbind(Estimate = logReg[[y]]$coef, confint(logReg[[y]]))), keep.rownames=TRUE))
  setnames(logRegResults[[y]],c("BED","Estimate","Lower2.5","Upper97.5"))
  y <- y+1
  
  # Logistical regression on FPs
  logReg[[y]] <- glm(FP ~ gc15 + gc15to20 + gc20to25 + gc25to30 + gc55to60 + gc60to65 + gc65to70 + gc70to75 + gc75to80 + gc80to85 + gc85 + lt7.lt51bp + lt7.51to200bp + gt6.lt51bp + gt6.51to200bp + gt6.gt200bp + hg19SelfChainSplit	+ refseqUnionCds + hs37.d5.mask75_50 + hs37.d5.mask35_50 + structural2 + compositional1 + cpg.islands + SRHomopolymer.3to5 + SRHomopolymer.6to10 + SRDiTR.11to50 + SRDiTR.51to200 + SRTriTR.11to50 + SRTriTR.51to200 + SRQuadTR.11to50 + SRQuadTR.51to200, data=chartsFPTP[[z]], family = binomial(logit))
  summary(logReg[[y]])
  #this generates the ratios and confidence intervals
  logRegResults[[y]] <- as.data.table(cbind(File = paste(chartName, "False Positives"),exp(cbind(Estimate = logReg[[y]]$coef, confint(logReg[[y]]))), keep.rownames=TRUE))
  setnames(logRegResults[[y]],c("BED","Estimate","Lower2.5","Upper97.5"))
  y <- y+1
  
  # Logistical regression on FNs
  logReg[[y]] <- glm(FN ~ gc15 + gc15to20 + gc20to25 + gc25to30 + gc55to60 + gc60to65 + gc65to70 + gc70to75 + gc75to80 + gc80to85 + gc85 + lt7.lt51bp + lt7.51to200bp + gt6.lt51bp + gt6.51to200bp + gt6.gt200bp + hg19SelfChainSplit	+ refseqUnionCds + hs37.d5.mask75_50 + hs37.d5.mask35_50 + structural2 + compositional1 + cpg.islands + SRHomopolymer.3to5 + SRHomopolymer.6to10 + SRDiTR.11to50 + SRDiTR.51to200 + SRTriTR.11to50 + SRTriTR.51to200 + SRQuadTR.11to50 + SRQuadTR.51to200, data=chartsFNTP[[z]], family = binomial(logit))
  summary(logReg[[y]])
  logRegResults[[y]] <- as.data.table(cbind(File = paste(chartName, "False Negatives"),exp(cbind(Estimate = logReg[[y]]$coef, confint(logReg[[y]]))), keep.rownames=TRUE))
  setnames(logRegResults[[y]],c("BED","Estimate","Lower2.5","Upper97.5"))
  y <- y+1
  
  # Logistical regression on NAs
  logReg[[y]] <- glm(NotA ~ gc15 + gc15to20 + gc20to25 + gc25to30 + gc55to60 + gc60to65 + gc65to70 + gc70to75 + gc75to80 + gc80to85 + gc85 + lt7.lt51bp + lt7.51to200bp + gt6.lt51bp + gt6.51to200bp + gt6.gt200bp + hg19SelfChainSplit	+ refseqUnionCds + hs37.d5.mask75_50 + hs37.d5.mask35_50 + structural2 + compositional1 + cpg.islands + SRHomopolymer.3to5 + SRHomopolymer.6to10 + SRDiTR.11to50 + SRDiTR.51to200 + SRTriTR.11to50 + SRTriTR.51to200 + SRQuadTR.11to50 + SRQuadTR.51to200, data=chartsNATP[[z]], family = binomial(logit))
  summary(logReg[[y]])
  logRegResults[[y]] <- as.data.table(cbind(File = paste(chartName, "Not Assessed"),exp(cbind(Estimate = logReg[[y]]$coef, confint(logReg[[y]]))), keep.rownames=TRUE))
  setnames(logRegResults[[y]],c("BED","Estimate","Lower2.5","Upper97.5"))
}

# Save the chart of logistic regression
write.csv(file="Logistic_Regression.csv", x=rbindlist(logRegResults))


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

# List of the bed file names
bedFiles = logRegResults[[1]][,BED]

#GRAPHING WAY 1 -------------

# Each bar is colored differently, no legend
# Save each VCF graph to a variable
graphs <- lapply(1:length(VCFNames),function(y){
  # Creates simple bar graph for 1 VCF, assigns the data to x and y
  ggplot(data=logRegResults[[y]], aes(x=factor(bedFiles, levels=bedFiles), y=Estimate, fill=factor(bedFiles, levels=bedFiles))) + 
    # Create a continous y axis
    geom_bar(stat="identity") + 
    # Remove legend
    guides(fill=FALSE) + 
    # Make the x axis horizontal
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1)) + 
    # Add error bars
    geom_errorbar(aes(ymax = Upper97.5, ymin = Lower2.5), position=position_dodge(width=0.9), width=0.25) + 
    # Add a title
    ggtitle(paste(unique(logRegResults[[y]][,VCF]),"Odds")) + 
    # Label the axes
    labs(x = "BED File", y = "Odds")
})

#----------------------------

#GRAPHING WAY 2 -------------

# Each bar is grouped by the BED File type, legend present
# Create a vector that assigns each BED file to a group
# **Change these values
bedGroups <- c()
bedGroups[1] <- "Intercept"
bedGroups[2:12] <- "GC Content"
bedGroups[13:17] <- "Imperfect Repeats"
bedGroups[20:23] <- "Unimask"
bedGroups[c(18,19,24)] <- "Miscellaneous"
bedGroups[25:32] <- "Simple Repeats"
# Create a vector for each group once in the order they appear along the x axis
totalGroups <- unique(bedGroups)

# Save each VCF graph to a variable
graphs <- lapply(1:length(VCFNames),function(y){
  # Creates simple bar graph for 1 VCF, assigns the data to x and y
  ggplot(data=logRegResults[[y]], aes(x=factor(bedFiles, levels=bedFiles), y=Estimate, fill=bedGroups)) + 
    # Create a continous y axis
    geom_bar(stat="identity") + 
    # Add legend for BED file groupings
    scale_fill_discrete(name="BED File Groupings", breaks=totalGroups) + 
    # Make the x axis horizontal
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1)) + 
    # Add error bars
    geom_errorbar(aes(ymax = Upper97.5, ymin = Lower2.5), position=position_dodge(width=0.9), width=0.25) + 
    # Add a title
    ggtitle(paste(unique(logRegResults[[y]][,VCF]),"Odds")) + 
    # Label the axes
    labs(x = "BED File", y = "Odds")
})

#----------------------------

# Print all 3 graphs separately
graphs[[1]]
graphs[[2]]
graphs[[3]]

# Plot all graphs on one page (click zoom to see)
multiplot(plotlist = graphs, cols=floor(sqrt(length(graphs))))
