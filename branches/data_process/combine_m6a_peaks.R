### combine_m6a_peaks.R ####################################################################
# Description
# Generates a matrix of combined m6a peaks 
### HISTORY #######################################################################################
# Version       Date            Developer           Comments
# 0.01          2018-12-15      Schen       
### NOTES #########################################################################################
# named.bed
# column 1 : chr
# column 2 : start
# column 3 : end
# column 4 : gene_name
# column 5 : score
# column 6 : strand
# column 7 : sample_name
# ...
### PREAMBLE ######################################################################################
library(getopt);
params <- matrix(
	c(
		'input_dir', 'i', '0', 'character',
		'pattern', 'p', '0', 'character',
		'stranded', 's', '1', 'logical', 
		'output_name', 'o', '0', 'character'
		),
	ncol = 4,
	byrow = TRUE
	);
opt <- getopt(params);
### set default value for stranded to FALSE
if (is.null(opt$stranded)){opt$stranded = FALSE};
if (is.null(opt$input_dir)){opt$input_dir = '/mnt/work1/users/hansengroup/Yong/m6a/peak/metpeak/peak_30M/'};
if (is.null(opt$pattern)){opt$pattern = '_peaks_sorted_named.bed'};
if (is.null(opt$output_name)){opt$output_name = 'test'};
### validate options
if(is.null(opt$input_dir)) {
	stop('No input directory provided')
} else if (is.null(opt$input_dir)) {
	stop('No input file pattern provided')
	} else if (is.null(opt$output_name)) {
	stop('No output name provided')
};
###
flist <- list.files(path = opt$input_dir, pattern=opt$pattern, full.names = TRUE);
fname <- gsub('.*/', '', flist);
fname <- gsub(opt$pattern, '', fname);
i <- 1;
tmp <- read.table(flist[i], as.is = TRUE)[, 1:7];
tmp <- unique(tmp[, c(1:4, 6)]);
rownames(tmp) <- apply(tmp[, c(1:5)], 1, paste0, collapse = '_');
mydata <- data.frame(matrix(nrow = nrow(tmp), ncol = length(fname)));
rownames(mydata) <- rownames(tmp);
colnames(mydata) <- fname;
mydata[rownames(tmp), fname[i]] <- 1;
for(i in seq(2, length(fname))){
	tmp <- read.table(flist[i], as.is = TRUE)[, 1:7];
	tmp <- unique(tmp[, c(1:4, 6)]);
	rownames(tmp) <- apply(tmp[, c(1:5)], 1, paste0, collapse = '_');
	mydata[rownames(tmp), fname[i]] <- 1;
};
mydata[is.na(mydata)] <- 0;
save(mydata, file = paste0(Sys.Date(), '_peaks_uniq.rda'));
#### create bed file, merge overlapping regions
myinfo <- strsplit(rownames(mydata), '_');
mybed <- data.frame(chr = gsub('_.*', '', rownames(mydata)), 
	start = sapply(myinfo, function(x) x[2]),
	end = sapply(myinfo, function(x) x[3]),
	name = sapply(myinfo, function(x) x[4]),
	score = 0,
	strand = sapply(myinfo, function(x) x[5])
	);
mybed$start <- as.numeric(as.vector(mybed$start));
mybed <- mybed[order(mybed$chr, mybed$start), ];
mybed$start <- as.vector(mybed$start);
mybed$end <- as.vector(as.numeric(as.vector(mybed$end)));
mybed <- cbind(mybed, mydata);
write.table(mybed, paste0(Sys.Date(), '_peaks_uniq.bed'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t');
cmd <- paste0('/mnt/work1/software/bedtools/2.27.1/bin/bedtools merge -s -c 4,5,6,', paste0(seq(7, 69), collapse=','), ' -o distinct,distinct,distinct,', paste0(rep('max', 63), collapse = ','), ' -i ', paste0(Sys.Date(), '_peaks_uniq.bed'), '> ', paste0(Sys.Date(), '_peaks_merge.bed'));
system2('/bin/bash', args = c('-c', shQuote(cmd)));
cmd1 <- paste0('cut -f 1,2,3,4,5,6 ', paste0(Sys.Date(), '_peaks_merge.bed'), '> ', paste0(Sys.Date(), '_peaks_merge_anno.bed'))
system2('/bin/bash', args = c('-c', shQuote(cmd1)));
#### read in the new file
mdat <- read.table(paste0(Sys.Date(), '_peaks_merge.bed'), as.is = TRUE);
colnames(mdat)[7:ncol(mdat)] <- fname;
save(mdat, file = paste0(Sys.Date(), '_peaks_merge.rda'));
