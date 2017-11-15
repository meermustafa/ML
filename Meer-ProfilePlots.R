# Meer Mustafa
# Oct 11, 2017

# set working directory

 # Figure 1a)
# 
# 1) check if the values in the matrixes and vectors are extractable
# if not:
#     2) update BED and BAM files
# 





bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  library("GenomicRanges")
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}



# # sort by different columns in the bed file and check if their outputs are the same
# - pick the most enriched enhancers by sorting by the -log10pvalue, -log10qvalue
# 
# - function to normalize bed files to an argument defined size e.g. 1000bp. use the summit position, to get the model V2+V10 then expand out n bp down and upstream



#### function
getBEDWithoutPromoters

setwd("~/Dropbox (Sanjana Lab)/SLab histone/ProfilePlots")
refGeneTSS = bed_to_granges('refGene_hg19_TSS.bed')
head(refGeneTSS)

# flank the TSS by 2.5 Kb on both sides
flanked2.5Kb_refGeneTSS = flank(refGeneTSS, width = 2500, both = T)

write.table(x = flanked2.5Kb_refGeneTSS, file = 'flanked2.5Kb_refGeneTSS.bed', sep = '\t',
            row.names = F, col.names = F, quote = F)



#### function
SortAndSummitOfBEDRange


### function
ChooseTopNGRanges









# - function to calculate the bam reads at this genomic interval =
#   RSamTools

library(Rsamtools)

# get reads from BAM file in a certain region
bf <- BamFile(file = 'A549.LUNG_ADENOCARCINOMA.ENCFF122XDP.H3K27ac.bam')
# set region of genome in BAM file to search
param <- ScanBamParam(which = GRanges("chr8",IRanges(start=126000000, end=133000000)) )
param = ScanBamParam(which = A549_MYC_Bed_rG_minus2500bpFlankTSS[2,])

# quick summary
quickBamFlagSummary(bf, param=param)

testPileup = pileup(bf, scanBamParam=param)
max(testPileup$count)

class(AllBAMNamesWithABGA[1])
deparse(substitute(testPileup))
class(deparse(substitute(testPileup)))

# create a new R object name using paste to paste the BAM filename
newname = paste0(AllBAMNamesWithABGA[1], deparse(substitute(testPileup)), sep = '__')
newname
paste(AllBAMNamesWithABGA[1], deparse(substitute(testPileup)), sep = '__') = pileup(bf, scanBamParam=param) # pileupParam=PileupParam(max_depth=8000, min_mapq=0, min_base_quality=0)
head(testPileup)
tail(testPileup)
hist(testPileup$count,breaks = c(100))
AllBAMNamesWithABGA


# resave pileup column
paste(AllBAMNamesWithABGA[1], deparse(substitute(testPileup)), sep = '__')$count
paste("samtools view", AllBAMNamesWithABGA[1], " | wc -l",sep=" ")




# assign data frame a new name in a for loop
for (i in 1:5) {
  assign( paste0("b", i), data.frame(oneBAMSMatrix))
}
b1

for (i in 1:5) {
  assign(paste0("readPileup", 
                sapply(strsplit(AllBAMNamesWithABGA[1], split='.', fixed=TRUE), function(x) (x[1]))# split the bam file name to keep only cell line
  ), 
  data.frame(testPileup)) # what to assign to the named data frame
}



# FUNCTION TO GET READ PILEUP FROM BAM WITHIN THE REGION OF A GRANGE
getBAMReadsInAGRange = function (inputBAMFileName, inputGRangeObject, bpSizeToFlankGRangeBy, binSizebp) {
  # print out passed-in args
  print(inputBAMFileName)
  print(inputGRangeObject)
  print(deparse(substitute(inputBAMFileName)))
  
  # read in BED file as GRange object using genomation's readGenereic function
  readGeneric()
  
  # set the bam file name i.e. character
  bf <- BamFile(file = inputBAMFileName)
  
  
  # iterate through the GRange ranges and get the read pileup within the range
  for (rowIndex in 1:length(inputGRangeObject)) {
    
    # set the GRange in which to take the reads from
    param = ScanBamParam(which = inputGRangeObject[rowIndex, ])
    
    # use pileup from RSamTools to compute reads in the BED bounds
    # this returns multiple read counts for positions that have reads on both DNA strands 
    rowPileup = pileup(bf, scanBamParam=param)
    
    
    consolidate / aggregate to pos and take the max read count (keep 1 read count per DNA position)
    save as readPileup
    
    
    
    should have a single row of scored bins for each Grange Enhancer
    rbind to a master readPileup
    

  }
  
  # normalize the score by calling samtools from command line
  divide by library size and save in a new column inside the existing df
  # calculate library size and divide by 1 million 1e6
  CurrentBAMLibrarySizeInMillionsOfReads =  ( as.numeric(system(paste("samtools view", 
                                                                      inputBAMFileName, 
                                                                      " | wc -l",sep=" "), 
                                                                intern = T)) 
                                              / 1000000 )
  print(CurrentBAMLibrarySizeInMillionsOfReads)
  
  
  
  
  # SAVE THE DF USING THE NAME OF THE SPLIT NAME OF THE BAM FILE THAT IS LOADED IN
  assign(paste0("readPileup", 
                sapply(strsplit(AllBAMNamesWithABGA[1], split='.', fixed=TRUE), function(x) (x[1]))# split the bam file name to keep only cell line
                ), 
         data.frame(Pileup)) # what to assign to the named data frame
  
  
}





### modular testing!!!


# give individual vectors a name (e paste0 e.g. e1, e2, e3)
# rbind to larger matrix that contains all enhancers for one BAM (name this using the BAM file input name)
# output should be a collection of matrices named after BAM file that was used for pileup counts



list.files(pattern = '.bam$')



# iterate through the GRange ranges and get the read pileup within the range
for (rowIndex in 1:length(inputGRangeObject)) {
  
  # set the GRange in which to take the reads from
  param = ScanBamParam(which = inputGRangeObject[rowIndex, ])
  
  # use pileup from RSamTools to compute reads in the BED bounds
  # this returns multiple read counts for positions that have reads on both DNA strands 
  rowPileup = pileup(bf, scanBamParam=param)
  
  
  consolidate / aggregate to pos and take the max read count (keep 1 read count per DNA position)
  save as readPileup
  
  
  
  should have a single row of scored bins for each Grange Enhancer
  rbind to a master readPileup
  
  
}

# normalize the score by calling samtools from command line
divide by library size and save in a new column inside the existing df
# calculate library size and divide by 1 million 1e6
CurrentBAMLibrarySizeInMillionsOfReads =  ( as.numeric(system(paste("samtools view", 
                                                                    inputBAMFileName, 
                                                                    " | wc -l",sep=" "), 
                                                              intern = T)) 
                                            / 1000000 )
print(CurrentBAMLibrarySizeInMillionsOfReads)





# singular code blocks!!!

setwd("~/Dropbox (Sanjana Lab)/SLab histone/ProfilePlots")


require(genomation)
if(require(randomForest1) == T) {
  cat('successfully loaded genomation library \n')
}

# put a counter on the bams
for (i in 1:length(list.files(pattern = '.bam$'))) {
  # set all the BAM file names as a vector to query from
  AllBamNamesInDirectory = list.files(pattern = '.bam$')
  
  cat ('Working on the' ,i, "bam. It's name is", AllBamNamesInDirectory[i] ,'\n')
  
}




AllBamNamesInDirectory = list.files(pattern = '.bam$')
AllBamNamesInDirectory

library(genomation)
A549_MYC_Bed_rG_minus2500bpFlankTSS = 
              readGeneric('lung_adenocarcinoma.cancer.bed.sorted.enrichmentRanked.minus2500bpFlankingTSS',
              meta.cols = list(score = 5),
              header = F, keep.all.metadata = T)


# make an artificial GRange object
testGrange = GRanges( seqnames = 'chrX', ranges = IRanges(start = 9995288, end = 9996288))
testGrange
testGrange1 = GRanges( seqnames = 'chrX', ranges = IRanges(start = 73466361, end = 73467361))
testGrange1
testGrange1 = GRanges( seqnames = 'chrX', ranges = IRanges(start = 71839346, end = 71840346))
testGrange1

# combine GRange
combinedtestrange= c(testGrange, testGrange1)
combinedtestrange

# make artificial grange
gr <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(start = 101:110, end = 111:120, 
                   names = paste0('e', 1:10)) # name the enhancers
  )
gr

head(A549_MYC_Bed_rG_minus2500bpFlankTSS)











# load in the BED file - 
require(genomation)

# read in the cancer-wide compiled enhancer BED
readGeneric()

list.files(pattern = 'cancer.bed.sorted$')

read in all cancer beds
concatenate them all using GRange c

# grab only the granges
granges(gr)

# grabonly the mcols
mcols(gr)


### concatenate all H3K27ac cancer BEDS
# 1 use 11 from ProfilePlots
list.files(pattern = 'cancer.bed.sorted$')

# 2 use several 100s GEO-wide
AllGEOCistromH3K27acBEDs = list.files(path = '../Pancancer_H3K27ac/All BED files/', pattern = '.bed')
# pick only the cancer samples
db = read.table('../Pancancer_H3K27ac/All BED files/all_BED_peak_lengths_and_counts_w_cellline_treatment.txt',header = T)
head(db)
which(db$cancer_status == 'yes')
AllGEOCistromH3K27acBEDs[which(db$cancer_status == 'yes')]

# loop through all the files, read them in as GRange Objects, concatenate all the 


# use all H3K27ac BEDs
AllCancerBEDs = system(paste("ls *cancer.bed.sorted"), intern = T)
AllCancerBEDs

system(paste("ls", 
             inputBAMFileName, 
             " | wc -l",sep=" "), 
       intern = T))



lung_adenocarcinoma.cancer.bed.sorted.enrichmentRanked.minus2500bpFlankingTSS


getBAMReadsWithinBinnedEnhancers = function (inputBEDfile, topNEnhancersToSelect, enhancerBinSizebp) {
  
  # load in the BED file - 
  require(genomation)
  
  # read in the cancer-wide compiled enhancer BED
  fullBEDGRange = readGeneric(inputBEDfile,
                              meta.cols = list(score = 5),
                              header = F, keep.all.metadata = T)
  
  n = topNEnhancersToSelect
  top_fullBEDGRange = fullBEDGRange[1:n]
  
  
  # loop through the bam files. make sure each .bam has a .bai index
  for (BAM in list.files(pattern = '.bam$')) {
    
    
    # set names for all BAMs in current directory
    AllBAMNamesInDirectory = list.files(pattern = '.bam$')
    
    # assign bam file name
    bf = BamFile(file = BAM) 
    
    # initialise matrix to hold all enhancers
    currentCancerMatrix = matrix(
      nrow = length(top_fullBEDGRange), 
      ncol = length(seq(start(top_fullBEDGRange[1,]), end(top_fullBEDGRange[1,]) , enhancerBinSizebp)))
    
    cat('The dimensions of the current cancers df is:',dim(currentCancerMatrix))
    
    # loop through the enhancers in the BED file 
    for (BEDRowIndex in 1:length(top_fullBEDGRange)) {
      
      # set chr
      binPositionChr = seqnames(top_fullBEDGRange[BEDRowIndex, ])
      
      # set the edge positions of the bins in the enhancer
      binPositionEdges = seq(start(top_fullBEDGRange[BEDRowIndex, ]), end(top_fullBEDGRange[BEDRowIndex, ]) , enhancerBinSizebp)
      
      # set the singular enhancer region as a variable to be used later in the ScanBamParam
      singularEnhancer = top_fullBEDGRange[BEDRowIndex, ]
      cat('This is the singular enhancer that is about to be analysed:') ; print(singularEnhancer)
      
      
      # initialise an empty vector to hold the MAX read pileup
      #### NAME THIS E1, E2, ETC... IN PARALLEL WITH THE ROW NUMBER BEING PROCESSED
      
      # assign( paste0('e', BEDRowIndex) # add the number of the enhancer
      #         , value = c() ) # should be empty vector
      
      binPileupForCurrentEnhancer = c()
      
      # segment enhancers and compute pileups in bins of the enhancers
      # use length of loop to stop at length(bins - 1)
      for (binIndex in 1: (length( seq(start(singularEnhancer), end(singularEnhancer) , enhancerBinSizebp)) - 1) ) {
        
        print(binIndex)
        
        # use the i and i + 1 bin. Starting with bin 1 and 2. Stop at length(bins - 1)
        # set this GRange as the ScanBamParam for RSamTools to find the pileups per bp in this genomic range
        param = ScanBamParam(which = (seqnames = binPositionChr ,start = binPositionEdges[i], end = binPositionEdges[i + 1]))
        cat('This is the singular enhancer that is about to be analysed:') ; print(singularEnhancer)
        
        # compute pileup in the current GRange
        binPileup = pileup(bf, scanBamParam = param)
        
        # Take the max pileup in the binned genomic region
        maxCountInBin = max (binPileup$count)
        
        # append this bin's pileup to the running vector of max read count pileups for this enhancer
        # add a value of the next binned pilup to the vector
        
        # assign( paste0('e', BEDRowIndex) # add the number of the enhancer
        #         , append( paste0('e', BEDRowIndex), maxCountInBin) ) # adding to the vector
        
        binPileupForCurrentEnhancer = append(binPileupForCurrentEnhancer, maxCountInBin)
        
        
        
        
      } # end of binIndex looping
      
      # append the most recent enhancer's binned pileups () to the full cancer matrix
      currentCancerMatrix[BEDRowIndex, ] = binPileupForCurrentEnhancer
      cat('The dimensions of the current cancers matrix after adding the',BEDRowIndex,'enhancer is:',dim(currentCancerMatrix))
      
      
    } # end of enhancer looping 
    
    # normlize the read count to the library size
    CurrentBAMLibrarySizeInMillionsOfReads =  ( as.numeric(system(paste("samtools view", 
                                                                        BAM, 
                                                                        " | wc -l",sep=" "), 
                                                                  intern = T)) 
                                                / 1000000 ) # divide by 1e6 to get Reads Per Million
    cat('The RPM value for',BAM,'is',print(CurrentBAMLibrarySizeInMillionsOfReads))
    
    
    # divide the matrix by this RPM value
    currentCancerMatrix = currentCancerMatrix / CurrentBAMLibrarySizeInMillionsOfReads
    cat('The dimensions of the current cancers matrix after normalizing to library size is:',dim(currentCancerMatrix))
    
    # normalize values from 0 to 1
    #normalizedBedGraphs = apply(TransposedUniformBGAWindowedValuesOnly, 2, function(x) ((x - min(x)) / (max(x) - min(x))) )
    
    
    # after looping through all the enhancer GRanges change the name of the currentCancerMatrix 
    # to a name that has the BAM cancer cell line name
    assign( paste0( sapply( strsplit( BAM, split='.', fixed=TRUE ), # pull name from BAM file name currently loaded in
                            function(x) (x[1])) # split the bam file name to keep only cell line
                    , 'binnedEnhancerPileups') # add to the name that this is the binnedEnhancerPileups
            , value = currentCancerMatrix ) # should be empty df on 1st iteration, should have
    
    
    # loops to the next BAM file now
    
    
  } # end of BAM looping
  
  
} # end of function



#turn the below for loop into a function to change the bin sizes
getBAMReadsWithinBinnedEnhancers = function (inputBEDfile, topNEnhancersToSelect, enhancerBinSizebp) {
  




AFTER RUNNING FIRST PILEUP, CHECK THE HEATMAP AND REORDER THE ROWS (ENHANCERS) 
TO GO FROM MOST ENRICHED (HIGHEST ROWSUMS OR ROWMEANS) to lowest
save this grange again and use this going forward
---- USE WHILE LOOP TO DO THIS. WHILE REORDERBEDCOUNTER < 2 so it only runs once

orrrrr just do this before looping through all BAMs



# check to see if all the pileup labels come from the same GRange
length(which(duplicated(rowPileup$which_label))) == nrow(rowPileup)-1
which(duplicated(rowPileup$pos))
rowPileup[185:190,]






# DON'T NEED TO COMBINE SEPARATE 
##### UNLESS WANT TO CBIND ONE BAM ONTO THE OTHERS AND THEN IN PHEATMAP, USE A LINE BREAK TO DELINEATE DIFFERENT BAM SAMPLES
# concatenating all of the matrices/dfs from the list into one data frame
docalledMatrix = do.call(rbind, list1)
docalledMatrix = rbindlist(oneBAMSMatrix, oneBAMSMatrixv2)
docalledMatrix






#strsplit(AllBAMNamesWithABGA[1], split='.', fixed=TRUE, FUN = function(x) (x[1])
sapply(strsplit(AllBAMNamesWithABGA[1], split='.', fixed=TRUE), function(x) (x[1]))

bamNameStrings[1]

getBAMReadsInAGRange(AllBAMNamesWithABGA[1], A549_MYC_Bed_rG_minus2500bpFlankTSS)

assign( paste0( sapply(strsplit(AllBAMNamesInDirectory[2], split='.', fixed=TRUE), # pull name from BAM file name
                       function(x) (x[1])) # split the bam file name to keep only cell line)
  , 'binnedEnhancerPileups') # specify that this is the name of binnedEnhancerPileups
  , matrix())

DND41binnedEnhancerPileups
DND41binnedEnhancerPileups = rbind(DND41binnedEnhancerPileups, oneBAMSMatrix)
DND41binnedEnhancerPileups

assign(paste0("b", 
              sapply(strsplit(AllBAMNamesWithABGA[2], split='.', fixed=TRUE), # pull name from BAM file name
                     function(x) (x[1])) # split the bam file name to keep only cell line)
              , matrix(oneBAMSMatrix)) # what to assign to the named data frame
       )

       
for loop through the BAM file names using list.files
call the function that was just created


# get the names
bamNameStrings = sapply(strsplit(list.files(path = "~/Dropbox (Sanjana Lab)/SLab histone/ProfilePlots/",
                                            pattern = '.bam$'),
                                 split='.', fixed=TRUE), function(x) (x[1]))
bamNameStrings


for (BAMname in (list.files(path = "~/Dropbox (Sanjana Lab)/SLab histone/ProfilePlots/",
                                      pattern = '.bam$'))) {
  
  print(BAMname)
  
  
  
}










# DON'T NEED TO INCLUDE BGA NAMES

AllbgaNameStrings = list.files(path = "~/Dropbox (Sanjana Lab)/SLab histone/Pancancer_H3K27ac/ENCODE_CancerVsNormalPairs/BAMs/cancer and normal",
                               pattern = '*bga')
AllbgaNameStrings
bgaNameStrings = sapply(strsplit( list.files(path = "~/Dropbox (Sanjana Lab)/SLab histone/Pancancer_H3K27ac/ENCODE_CancerVsNormalPairs/BAMs/cancer and normal",
                                             pattern = '*bga'),
                                  split='.', fixed=TRUE), function(x) (x[1]))
bgaNameStrings



bamNameStrings = sapply(strsplit(list.files(path = "~/Dropbox (Sanjana Lab)/SLab histone/ProfilePlots/",
                                            pattern = '.bam$'),
                                 split='.', fixed=TRUE), function(x) (x[1]))
bamNameStrings

commonBGAandBAMnames = bamNameStrings %in% bgaNameStrings
commonBGAandBAMnames

# use the system command to call samtools view <bam> | wc -l to normalize the number of readsto library size
AllBAMNames = list.files(path = "~/Dropbox (Sanjana Lab)/SLab histone/ProfilePlots/",
                               pattern = '.bam$')
AllBAMNames

AllBAMNamesWithABGA=list.files(path = "~/Dropbox (Sanjana Lab)/SLab histone/ProfilePlots/",
                               pattern = '.bam$') [commonBGAandBAMnames]
AllBAMNamesWithABGA


# singular system commands
# system('samtools view A549.LUNG_ADENOCARCINOMA.ENCFF122XDP.H3K27ac.bam | wc -l')
test = system(paste("samtools view", AllBAMNamesWithABGA[1], " | wc -l",sep=" "), intern = T)
test

# for loop to get all library sizes
BAMLibrarySizeInMillionsOfReads = c()
for (BAMNamesWithABGA in AllBAMNamesWithABGA) {
  print(BAMNamesWithABGA)
  print( class(BAMNamesWithABGA))
  
  # calculate library size and divide by 1 million 1e6
  CurrentBAMLibrarySizeInMillionsOfReads =  ( as.numeric(system(paste("samtools view", 
                                                                      BAMNamesWithABGA, 
                                                                      " | wc -l",sep=" "), 
                                                                intern = T)) 
                                              / 1000000 )
  print(CurrentBAMLibrarySizeInMillionsOfReads)
  
  # append to vector
  BAMLibrarySizeInMillionsOfReads = append(CurrentBAMLibrarySizeInMillionsOfReads, BAMLibrarySizeInMillionsOfReads)
  print(BAMLibrarySizeInMillionsOfReads)
  
}
BAMLibrarySizeInMillionsOfReads

# appended in to front so reverse the order
BAMLibrarySizeInMillionsOfReads = rev(BAMLibrarySizeInMillionsOfReads)

BAMLibrarySizeInMillionsOfReads # number in millions of reads in order of listed files in the directory
AllBAMNamesWithABGA # all BAM names that have the library size compute








library(Rsamtools)

# initialise list to hold all normalized data frames



# loop through the BAM names in order to get their per bp
# set region of genome in BAM file to search
####### ONCE THE BED REGIONS ARE DEFINED (E.G. ENHANCER SITES ARE DEFINED) 
# PASS IN LINE-BY-LINE THE GRANGE OBJECT (rows of the readGeneric object) 
# and get the read counts for that row, then save that vector and 

param <- ScanBamParam(which = GRanges("chr8", IRanges(start = 126000000, end = 133000000)))



for (i in 1:length(AllBAMNamesWithABGA)) {
  print(i)
  
  # get reads from BAM file in a certain region
  bf <- BamFile(file = AllBAMNameWithABGA)
  
  # quick summary
  quickBamFlagSummary(bf)
  
  # compute per bp pileup
  testPileup = pileup(bf ) # to get a certain genome region >scanBamParam = param  # pileupParam=PileupParam(max_depth=8000, min_mapq=0, min_base_quality=0)
  head(testPileup)
  tail(testPileup)
  hist(testPileup$count,breaks = c(100))
  
  # resave pileup column
  
}






#### final ------

# load in the BED file - 
require(genomation)

# read in the cancer-wide compiled enhancer BED
fullBEDGRange = readGeneric('lung_adenocarcinoma.cancer.bed.sorted.enrichmentRanked.minus2500bpFlankingTSS',
                            meta.cols = list(score = 5),
                            header = F, keep.all.metadata = T)

n = 5000
top_fullBEDGRange = fullBEDGRange[1:n]

# make all the widths the same
width(top_fullBEDGRange) = 1000
width(top_fullBEDGRange) == 100
head(top_fullBEDGRange)


# take V10 summit position, add it to start of range, 
# use this value and set it as both the start and the end, 
# flank it by 1000 or 2000 bp (use mean of median of width of all the ranges)




combinedtestrange




# create function to:
# read in a BED file (likely the one that will be the concantenated pancancer H3K27ac)
# 






# loop through the bam files. make sure each .bam has a .bai index
for (BAM in list.files(pattern = '.bam$')) {
  
  # assign bam file name
  bf = BamFile(file = BAM) 
  
  # initialise matrix to hold all enhancers
  currentCancerMatrix = matrix(
    nrow = length(combinedtestrange), 
    ncol = length(seq(start(top_fullBEDGRange[1,]), end(top_fullBEDGRange[1,]) , 50)) - 1 ) # 
  
  cat('The dimensions of the current cancers df is:',dim(currentCancerMatrix))
  
  # loop through the enhancers in the BED file 
  for (BEDRowIndex in 1:length(top_fullBEDGRange)) {
    
    # set chr
    binPositionChr = as.character(seqnames(top_fullBEDGRange[BEDRowIndex, ]))
    
    # set the edge positions of the bins in the enhancer
    binPositionEdges = seq(start(top_fullBEDGRange[BEDRowIndex, ]), end(top_fullBEDGRange[BEDRowIndex, ]) , 50)
    
    # set the singular enhancer region as a variable to be used later in the ScanBamParam
    singularEnhancer = top_fullBEDGRange[BEDRowIndex, ]
    #cat('This is the singular enhancer that is about to be analysed:') ; print(singularEnhancer)
    
    # initialise an empty vector to hold the MAX read pileup
    #### NAME THIS E1, E2, ETC... IN PARALLEL WITH THE ROW NUMBER BEING PROCESSED
    binPileupForCurrentEnhancer = c()
    
    # segment enhancers and compute pileups in bins of the enhancers
    # use length of loop to stop at length(bins - 1)
    for (binIndex in 1: (length( seq(start(singularEnhancer), end(singularEnhancer) , 50)) - 1))  {
      
      cat('This is the bin index that pileup will be computed for:',binIndex,'\n')
      
      # use the i and i + 1 bin. Starting with bin 1 and 2. Stop at length(bins - 1)
      # set this GRange as the ScanBamParam for RSamTools to find the pileups per bp in this genomic range
      param = ScanBamParam(which = GRanges(seqnames = binPositionChr, IRanges(start = binPositionEdges[binIndex], end = binPositionEdges[binIndex + 1])))
      
      # compute pileup in the current GRange
      binPileup = pileup(bf, scanBamParam = param)
      
      # Take the max pileup in the binned genomic region
      maxCountInBin = round(max (binPileup$count), digits = 1)
      cat('the max read count in the bin is:',maxCountInBin,'\n')
      
      # if the maxCountInBin is a -Inf or a non-numeric value then convert it to 0
      if (maxCountInBin == -Inf) {
        cat('Found a wrong value, replacing it with a 0 read count!\n')
        maxCountInBin = 0
      }
      
      cat('the NEW max read count in the bin is:',maxCountInBin,'\n')
      
      # append this bin's pileup to the running vector of max read count pileups for this enhancer
      # add a value of the next binned pilup to the vector
      binPileupForCurrentEnhancer = append(binPileupForCurrentEnhancer, maxCountInBin)
      
    } # end of bins
    
    cat('about the append this enhancers bins to the full matrix \n')
    cat('the size of the binned enhancer vector is',length(binPileupForCurrentEnhancer),'and the # of columns of the matrix is', ncol(currentCancerMatrix),'\n')
    
    # check to make sure that the length of bin values should == # of cols of the matrix
    if (length(binPileupForCurrentEnhancer == ncol(currentCancerMatrix))) {
      
      # append the most recent enhancer's binned pileups () to the full cancer matrix
      currentCancerMatrix[BEDRowIndex, ] = binPileupForCurrentEnhancer
      cat('The new dimensions of the cancers matrix after adding the',BEDRowIndex,'enhancer is:',dim(currentCancerMatrix),'\n')
    }
    
  } # end of all enhancers
  
  cat('About to compute library size of',BAM)
  # normlize the read count to the library size
  CurrentBAMLibrarySizeInMillionsOfReads =  ( as.numeric(system(paste("samtools view", 
                                                                      BAM, 
                                                                      " | wc -l",sep=" "), 
                                                                intern = T)) 
                                              / 1000000 ) # divide by 1e6 to get Reads Per Million
  cat('The RPM value for',BAM,'is',print(CurrentBAMLibrarySizeInMillionsOfReads))
  
  # divide the matrix by this RPM value
  currentCancerMatrix = currentCancerMatrix / CurrentBAMLibrarySizeInMillionsOfReads
  cat('The dimensions of the current cancers matrix after normalizing to library size is:',dim(currentCancerMatrix),'\n')
  
  # normalize values from 0 to 1
  cat('About to normalize values in matrix from 0-1 range!')
  # take values rowwise and normalize, transpose to keep in row format
  currentCancerMatrix = t(apply(currentCancerMatrix, 1, function(x) ((x - min(x)) / (max(x) - min(x))) ))
  print(currentCancerMatrix)
  
  # after looping through all the enhancer GRanges change the name of the currentCancerMatrix 
  # to a name that has the BAM cancer cell line name
  assign( paste0( sapply( strsplit( BAM, split='.', fixed=TRUE ), # pull name from BAM file name currently loaded in
                          function(x) (x[1])) # split the bam file name to keep only cell line
                  , 'binnedEnhancerPileups') # add to the name that this is the binnedEnhancerPileups
          , value = currentCancerMatrix ) # should be empty df on 1st iteration, should have
  
  # print statement to notify completion of the loop for each BAM
  cat('Successfully computed binned BED range read pileup for',
      paste0( sapply( strsplit( BAM, split='.', fixed=TRUE ), # pull name from BAM file name currently loaded in
                      function(x) (x[1]))),'BAM file. \n')
  
}


paste0( sapply( strsplit( AllBamNamesInDirectory[1], split='.', fixed=TRUE ), # pull name from BAM file name currently loaded in
                function(x) (x[1])) # split the bam file name to keep only cell line
        , 'binnedEnhancerPileups')




# plotting -----

# COMBINE SEPARATE MATRICES
##### UNLESS WANT TO CBIND ONE BAM ONTO THE OTHERS AND THEN IN PHEATMAP, USE A LINE BREAK TO DELINEATE DIFFERENT BAM SAMPLES
# concatenating all of the matrices/dfs from the list into one data frame
dfList = list(df1, df2)
for (binnedEnhancer in ls(pattern = '*binnedEnhancer*')) {
  print(binnedEnhancer)
  #dfList[[i]] = binnedEnhancer
  
}
list1
docalledMatrix = do.call(rbind, list1)
docalledMatrix = rbindlist(oneBAMSMatrix, oneBAMSMatrixv2)
docalledMatrix




# cbind all of the cancers and plot on heatmap
ls(pattern = '*binnedEnhancer*')
allCancersMatrix = data.frame()
for (singleEnhancer in ls(pattern = '*binnedEnhancer*')) {
  
  # convert the vector name to be able to call that name in the workspace
  singleEnhancerMatrix = eval(expr = as.name(singleEnhancer))
  print(singleEnhancerMatrix)
  # binnedEnhancerMatrix = as.name(binnedEnhancer)
  # print(class(binnedEnhancerMatrix))
  allCancersMatrix = cbind(allCancersMatrix, singleEnhancerMatrix)
  # print(allCancersMatrix)
  
}

allCancersMatrix

for(dfname in ls(pattern = '*binnedEnhancer*')){ 
  eval(substitute({ 
    print(dfname)
  }, 
  list( 
    df=as.name(dfname) 
  ) 
  ) 
  ) 
  
} 


allCancersMatrix



allCancersMatrix = rbind(
  A549binnedEnhancerPileups,
  K562binnedEnhancerPileups,
  MCF7binnedEnhancerPileups,
  MDA_MB_231binnedEnhancerPileups
)

allCancersMatrix

library(pheatmap)
pheatmap(allCancersMatrix,
         show_rownames = T,
         cluster_rows = F,
         cluster_cols = F)











