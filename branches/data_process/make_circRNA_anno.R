### make_circRNA_anno.R ####################################################################
# Description
# Generates bed annotation file for circRNA from the ref all file
### HISTORY #######################################################################################
# Version       Date            Developer           Comments
# 0.01          2018-12-19      Schen       
### NOTES #########################################################################################
#
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general);
library(getopt);
conf <- read.config.file('~/landscape/rnaseq_landscape/master_config_rnaseq.R');
params <- matrix(
        c(
                'input_rpkm', 'i', '0', 'character',
                'output_name', 'o', '0', 'character'),
        ncol = 4,
        byrow = TRUE
        );
opt <- getopt(params);
ref <- read.csv(conf$ref_all, row.names = 1, as.is = TRUE);
rpkm.circ <- read.csv(opt$input_rpkm, row.names = 1, as.is = TRUE);
ref <- ref[rownames(rpkm.circ), ];
write.csv(ref, paste0(Sys.Date(), '_circ_annotation_', opt$output_name, '.csv'));
ref.annot <- data.frame(chr = ref$chr, start = ref$start, end = ref$end, name = rownames(ref), score = 0, strand = ref$strand);
write.table(ref.annot, paste0(Sys.Date(), '_coordinate_circRNA_', opt$output_name, '.bed'), quote = FALSE,
	row.names = FALSE, col.names = FALSE, sep = '\t');
