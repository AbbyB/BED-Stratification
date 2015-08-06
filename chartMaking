
# Program to create a chart with the BED files along the top and the SNPs variants along the y axis
#  with binary output (1 or 0 in each cell) of whether the variant is in each BED file
# Input: All of the filename extension '.tabular' files in the current directory, 
#         which should be tab delimited VCF files that have been annotated for the following BED files:
#          gc15, gc15to20, gc20to25, gc25to30, gc30to55, gc55to60, gc60to65, gc65to70, gc70to75, gc75to80, gc80to85, gc85, 
#          lt7.lt51bp, lt7.lt101bp, lt7.51to200bp, lt7.gt200bp, gt6.lt51bp, gt6.lt101bp, gt6.51to200bp, gt6.gt200bp, 
#          um75.hs37.d5, um35.hs37.d5, hs37.d5.mask75_50, hs37.d5.mask35_50, structural, compositional, 
#          hg19SelfChainSplit, refseqUnionCds, cpg.islands, 
#          SRHomopolymer.3to5, SRHomopolymer.6to10, SRHomopolymer.gt10, SRDiTR.11to50, SRDiTR.51to200, SRDiTR.gt200, 
#          SRTriTR.11to50, SRTriTR.51to200, SRTriTR.gt200, SRQuadTR.11to50, SRQuadTR.51to200, SRQuadTR.gt200
# 
# Output: 1 CSV file that includes all the VCFs, 
#          moving left to right: leftmost column ("VCF") has the VCF file name without the extension, 
#          then columns with variant information, then chart columns
# 
# Run the commentted 'install.packages' lines the first time this program is run only
# When running this program with different BED or VCF files, take note of comments with '**'

# install.packages("data.table")  # only run this the first time to install
library(data.table)

# Input
# Finds all the files in the current directory with the file extension '.tabular'
files <- list.files(pattern = "\\.tabular$")

# Reads all files found in the directory
VCFs <- list()
for (y in 1:length(files)) {
  VCFs[[y]] <- fread(files[[y]], header=TRUE, sep="\t", colClasses="character")
}

# Labels the VCF files using the file name without the '.tabular' extension
VCFName <- gsub('.tabular','',files)
# **If you want different labels, change to:
#    VCFName <- c((label 1), (label 2), (etc. for the total number of files))

# Only keep columns that are important
# Variant information
# **Choose which variant info to keep 
variantInfo <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','GT')

# The names of the annotations in Galaxy, starting with GC annotations
annotations <- c('gc_content_15','gc_content_15to20','gc_content_20to25','gc_content_25to30','gc_content_30to55','gc_content_55to60',
                 'gc_content_60to65','gc_content_65to70','gc_content_70to75','gc_content_75to80','gc_content_80to85','gc_content_85',
                 # Unimask annotations
                 'um75_hs37d5','um35_hs37d5','hs37d5.mask75_50','hs37d5.mask35_50','structural2','compositional1',
                 # Miscellaneous annotations
                 'hg19_self_chain_split','refseq_union_cds','model_based_cpg',
                 # Simple repeats annotations
                 'SR_homopolymer_3to5','SR_homopolymer_6to10','SR_homopolymer_gt10','SR_diTR_11to50','SR_diTR_51to200','SR_diTR_gt200',
                 'SR_triTR_11to50','SR_triTR_51to200','SR_triTR_gt200','SR_quadTR_11to50','SR_quadTR_51to200','SR_quadTR_gt200',
                 # Imperfect and perfect repeats annotations
                 'lt7_lt51bp','lt7_lt101bp','lt7_51to200bp','lt7_gt200bp','gt6_lt51bp','gt6_lt101bp','gt6_51to200bp','gt6_gt200bp')
# **To add another BED file
#    1. Add the name of the annotation to the previous vector
#    Note: don't start names with a number or use dashes (-)

# Combine the names of all the columns to keep
keepColumns <- c(variantInfo,annotations)

# Create columns for the chart
addColumns <- c('gc15','gc15to20','gc20to25','gc25to30','gc30to55','gc55to60','gc60to65','gc65to70','gc70to75','gc75to80','gc80to85','gc85',
                'lt7.lt51bp','lt7.lt101bp','lt7.51to200bp','lt7.gt200bp','gt6.lt51bp','gt6.lt101bp','gt6.51to200bp','gt6.gt200bp',
                'um75.hs37.d5','um35.hs37.d5','hs37.d5.mask75_50','hs37.d5.mask35_50','structural','compositional',
                'hg19SelfChainSplit','refseqUnionCds','cpg.islands',
                'SRHomopolymer.3to5','SRHomopolymer.6to10','SRHomopolymer.gt10','SRDiTR.11to50','SRDiTR.51to200','SRDiTR.gt200',
                'SRTriTR.11to50','SRTriTR.51to200','SRTriTR.gt200','SRQuadTR.11to50','SRQuadTR.51to200','SRQuadTR.gt200')
# **To add another BED file
#    2. Add the name of the column of the chart to the previous vector
#    Note: don't start names with a number, don't use dashes (-), don't use the same name as the annotation

# Compute
# Loop through all the files
for (y in 1:length(VCFs)) {
  
  # Select the columns to keep
  VCFs[[y]] <- VCFs[[y]][,keepColumns,with=FALSE]
  # Add the chart columns
  VCFs[[y]][,addColumns := 0,with=FALSE]
  # Label the VCF file
  VCFs[[y]] <- cbind(VCF = VCFName[y], VCFs[[y]])
  
  # Create GC part of the chart
  VCFs[[y]][,gc15 := as.numeric(VCFs[[y]][,gc_content_15] != "")]
  VCFs[[y]][,gc15to20 := as.numeric(VCFs[[y]][,gc_content_15to20] != "")]
  VCFs[[y]][,gc20to25 := as.numeric(VCFs[[y]][,gc_content_20to25] != "")]
  VCFs[[y]][,gc25to30 := as.numeric(VCFs[[y]][,gc_content_25to30] != "")]
  VCFs[[y]][,gc30to55 := as.numeric(VCFs[[y]][,gc_content_30to55] != "")]
  VCFs[[y]][,gc55to60 := as.numeric(VCFs[[y]][,gc_content_55to60] != "")]
  VCFs[[y]][,gc60to65 := as.numeric(VCFs[[y]][,gc_content_60to65] != "")]
  VCFs[[y]][,gc65to70 := as.numeric(VCFs[[y]][,gc_content_65to70] != "")]
  VCFs[[y]][,gc70to75 := as.numeric(VCFs[[y]][,gc_content_70to75] != "")]
  VCFs[[y]][,gc75to80 := as.numeric(VCFs[[y]][,gc_content_75to80] != "")]
  VCFs[[y]][,gc80to85 := as.numeric(VCFs[[y]][,gc_content_80to85] != "")]
  VCFs[[y]][,gc85 := as.numeric(VCFs[[y]][,gc_content_85] != "")]
  
  # Create unimask part of the chart
  VCFs[[y]][,hs37.d5.mask75_50 := as.numeric(VCFs[[y]][,hs37d5.mask75_50] != "")]
  VCFs[[y]][,hs37.d5.mask35_50 := as.numeric(VCFs[[y]][,hs37d5.mask35_50] != "")]
  VCFs[[y]][,um75.hs37.d5 := as.numeric(VCFs[[y]][,um75_hs37d5] != "")]
  VCFs[[y]][,um35.hs37.d5 := as.numeric(VCFs[[y]][,um35_hs37d5] != "")]
  VCFs[[y]][,structural := as.numeric(VCFs[[y]][,structural2] != "")]
  VCFs[[y]][,compositional := as.numeric(VCFs[[y]][,compositional1] != "")]
  
  # Create miscellaneous part of the chart
  VCFs[[y]][,hg19SelfChainSplit := as.numeric(VCFs[[y]][,hg19_self_chain_split] != "")]
  VCFs[[y]][,refseqUnionCds := as.numeric(VCFs[[y]][,refseq_union_cds] != "")]
  VCFs[[y]][,cpg.islands := as.numeric(VCFs[[y]][,model_based_cpg] != "")]
  
  # Create simple repeats part of the chart
  VCFs[[y]][,SRHomopolymer.3to5 := as.numeric(VCFs[[y]][,SR_homopolymer_3to5] != "")]
  VCFs[[y]][,SRHomopolymer.6to10 := as.numeric(VCFs[[y]][,SR_homopolymer_6to10] != "")]
  VCFs[[y]][,SRHomopolymer.gt10 := as.numeric(VCFs[[y]][,SR_homopolymer_gt10] != "")]
  VCFs[[y]][,SRDiTR.11to50 := as.numeric(VCFs[[y]][,SR_diTR_11to50] != "")]
  VCFs[[y]][,SRDiTR.51to200 := as.numeric(VCFs[[y]][,SR_diTR_51to200] != "")]
  VCFs[[y]][,SRDiTR.gt200 := as.numeric(VCFs[[y]][,SR_diTR_gt200] != "")]
  VCFs[[y]][,SRTriTR.11to50 := as.numeric(VCFs[[y]][,SR_triTR_11to50] != "")]
  VCFs[[y]][,SRTriTR.51to200 := as.numeric(VCFs[[y]][,SR_triTR_51to200] != "")]
  VCFs[[y]][,SRTriTR.gt200 := as.numeric(VCFs[[y]][,SR_triTR_gt200] != "")]
  VCFs[[y]][,SRQuadTR.11to50 := as.numeric(VCFs[[y]][,SR_quadTR_11to50] != "")]
  VCFs[[y]][,SRQuadTR.51to200 := as.numeric(VCFs[[y]][,SR_quadTR_51to200] != "")]
  VCFs[[y]][,SRQuadTR.gt200 := as.numeric(VCFs[[y]][,SR_quadTR_gt200] != "")]
  
  # Create imperfect & perfect repeats part of the chart
  VCFs[[y]][,lt7.lt51bp := as.numeric(VCFs[[y]][,lt7_lt51bp] != "")]
  VCFs[[y]][,lt7.lt101bp := as.numeric(VCFs[[y]][,lt7_lt101bp] != "")]
  VCFs[[y]][,lt7.51to200bp := as.numeric(VCFs[[y]][,lt7_51to200bp] != "")]
  VCFs[[y]][,lt7.gt200bp := as.numeric(VCFs[[y]][,lt7_gt200bp] != "")]
  VCFs[[y]][,gt6.lt51bp := as.numeric(VCFs[[y]][,gt6_lt51bp] != "")]
  VCFs[[y]][,gt6.lt101bp := as.numeric(VCFs[[y]][,gt6_lt101bp] != "")]
  VCFs[[y]][,gt6.51to200bp := as.numeric(VCFs[[y]][,gt6_51to200bp] != "")]
  VCFs[[y]][,gt6.gt200bp := as.numeric(VCFs[[y]][,gt6_gt200bp] != "")]
  
  # **To add another BED file
  #    3. Add this line
  #        VCFs[[y]][,(name of the column in the chart (step 2)) := as.numeric(VCFs[[y]][,(name of the annotation in the VCF file (step 1))] != "")] - 
  
  # Removes the annotation columns
  VCFs[[y]][, ((length(variantInfo)+2):(length(keepColumns)+1)) := NULL, with = FALSE]
  
  # Sums each column in the chart 
  #  for (n in (length(variantInfo)+2):(length(variantInfo)+length(addColumns)+1)) {
  #    print(sum(VCFs[[y]][, n, with = FALSE]))
  #  }
  
}

# -------------------------------------------------

# Creates one large data table of all the VCFs
#  data <- rbindlist(VCFs)

# Creates one large data table of all VCFs and writes it to a CSV file
write.csv(file="VCF_and_BED_Chart.csv", x=rbindlist(VCFs))
