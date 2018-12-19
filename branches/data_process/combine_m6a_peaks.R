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
