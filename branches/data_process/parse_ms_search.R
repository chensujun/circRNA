### parse_ms_search.R ####################################################################
# Description
# Parse the MS/MS search results, save the xls file to rda
# initially characterization of circRNAs with peptide detection
# compare with circRNA associated with m6a
### HISTORY #######################################################################################
# Version       Date            Developer           Comments
# 0.01          2018-12-19      Schen       
### NOTES #########################################################################################
#
### PREAMBLE ######################################################################################
#### start working
library(BoutrosLab.plotting.general);
library(VennDiagram);
library(gdata);
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/m6a/translation/data");
source('~/landscape/rnaseq_landscape/analysis/my_functions/plot_functions.R');
source('~/landscape/rnaseq_landscape/analysis/my_functions/create.boxplot.R');
conf <- read.config.file('~/landscape/rnaseq_landscape/master_config_rnaseq.R');
mysearch <- read.xls('Circular miRNA search_PDX TMT data_Aug22_2018.xlsx');
save(mysearch, file = paste0(Sys.Date(), '_peptide_search.rda'));
mycirc <- gsub('pre=.*,post=.*|\\(|\\)', '', mysearch$Protein);
mycirc <- unlist(strsplit(mycirc, ';'));
mycirc <- unique(gsub('_[0-9]+$', '', mycirc));
##
circ.pp <- gsub('pre=.*,post=.*|\\(|\\)', '', mysearch$ProteinName);
circ.pp <- gsub('_[0-9]+$', '', circ.pp);
circ.pp <- data.frame(table(circ.pp));

### characteristics of circRNAs with peptide detection
#### length distribution
ref <- read.csv(conf$ref_m6a, as.is = TRUE, row.names = 1);
anno <- read.table('2018-12-19_coordinate_circRNA_m6a.bed', as.is = TRUE);
anno$len <- apply(ref[anno$V4, ], 1, function(x) sum(as.numeric(as.vector(unlist(strsplit(x[6], ','))))));
len.all <- density(log10(anno$len), bw=0.25);
len.pp <- density(log10(anno[anno$V4%in%mycirc, ]$len), bw = 0.25)
col_b <- default.colours(2)
col_f <- unlist(lapply(col_b, function(x) t_col(x, 50)))
pdf(generate.filename('density_length', 'circRNA_peptide', 'pdf'));
par(xaxs='i',yaxs='i', ps = 15, font.axis = 1, font.lab= 1)
plot(len.all, main = '', xlim = c(0, 6), ylim = c(0, 1), col = col_b[1], 
	xlab = expression('Length'), axes = FALSE);
axis(1, pos=0, at = seq(0, 6), labels = format(10^seq(0, 6), scientific=FALSE));
axis(2, pos=0, las = 2);
polygon(len.all, col = col_f[1], border = col_b[1]);
polygon(len.pp, col = col_f[2], border = col_b[2]);
legend(3.8, 0.95, border = col_b, fill = col_f, bty = 'n', legend = c('All', 'Peptide detection'));
dev.off();
#### mean abundance
rpkm.m6a <- read.csv('../../2018-06-12_circ_RPKM_m6a.csv', row.names = 1, as.is = TRUE);
anno$mean <- rowMeans(rpkm.m6a[rownames(anno), ]);
anno$ip <- rowMeans(rpkm.m6a[anno$V4, grep('IP', colnames(rpkm.m6a))]);
anno$ipt <- rowMeans(rpkm.m6a[anno$V4, grep('INPUT', colnames(rpkm.m6a))]);
wilcox.test(anno[anno$V4%in%mycirc, ]$ip, anno$ip);
wilcox.test(anno[anno$V4%in%mycirc, ]$ipt, anno$ipt);
#### overlap between circRNAs with peptide detection and m6a association 
m6a.circ <- read.table('2018-12-19_overlap_circRNA_m6a.bed', as.is = TRUE);
venn.diagram(main='',
        list('peptide circ' = mycirc,
        'm6a circ' = unique(m6a.circ$V4)
        ),
        fill=c("dodgerblue","goldenrod1"),
        filename = generate.filename('overlap_circRNA', 'peptide_m6a', 'tiff'),
        alpha=c(0.5,0.5),
        col=c("dodgerblue","goldenrod1"),
        main.cex=2,
        imagetype = "tiff",
        cex = 2,
        cat.cex = 2,
        cat.pos = c(-15, 15)
);
phyper(length(intersect(m6a.circ$V4, mycirc)), length(unique(m6a.circ$V4)), 
	nrow(rpkm.m6a)-length(unique(m6a.circ$V4)), length(mycirc), lower.tail = FALSE);
#### m6a peaks associated with peptide circRNA
load('2018-12-19_peaks_merge.rda');
anno.reg <- read.table('2018-12-17_peaks_merge_anno_hg19_name.bed');
rownames(mdat) <- paste0('peak_', seq(nrow(mdat)));
tmp <- data.frame(name = paste(mdat$V1, mdat$V2, mdat$V3, mdat$V4, mdat$V6, sep = '_'), id = rownames(mdat));
tmp <- tmp[tmp$name%in%anno.reg$V7), ];
anno.reg$id <- tmp[match(anno.reg$V7, tmp$name), ]$id;
mdat <- mdat[, 7:69];
anno.reg$nsamp <- rowSums(mdat[anno.reg$id, ]>0);
anno.reg$normal <- rowSums(mdat[anno.reg$id, grep('normal', colnames(mdat))]>0);
anno.reg$tumor <- rowSums(mdat[anno.reg$id, grep('tumor', colnames(mdat))]>0);
anno.reg$id_hg19 <- paste(anno.reg$V1, anno.reg$V2, anno.reg$V3, anno.reg$V4, anno.reg$V6, sep = '_');
m6a.circ$id_hg19 <- paste(m6a.circ$V7, m6a.circ$V8, m6a.circ$V9, m6a.circ$V10, m6a.circ$V12, sep = '_');
anno.reg$pp_circ <- factor(ifelse(anno.reg$id_hg19 %in% m6a.circ[m6a.circ$V4 %in% mycirc, ]$id_hg19, 1, 0));
circ.peaks <- data.frame(table(m6a.circ$V4));
colnames(circ.peaks) <- c('circ', 'npeak');
circ.peaks$pp_circ <- factor(ifelse(circ.peaks$circ%in%mycirc, 1, 0));
#### number of samples with m6a: all vs peptide circRNAs
pval <- scientific.notation(wilcox.test(anno.reg$nsamp~anno.reg$pp_circ)$p.value)
create.boxplot(
        formula = nsamp ~ pp_circ,
        data =  anno.reg,
        add.stripplot = TRUE,
        points.cex = .3,
        border.col = default.colours(2),
        lwd = 5,
        umb.lwd = 5,
        xlab.label = 'Peptide circRNAs',
        ylab.label = '#samples with m6A',
        xaxis.lab = c('No', 'Yes'),
        filename = generate.filename('peptide_circRNA','m6a_nsamples', 'tiff'),
        box.ratio = 0.6,
        ylab.cex = 1.5,
        xlab.cex = 1.5,
        width = 3,
        height = 4,
        style = 'Nature',
        key = list(
                        text = list(
                        lab = pval,
                        cex = 1,
                        col = 'black'
                        ),
                x = 0.25,
                y = 0.95
                ),
);
pval <- scientific.notation(wilcox.test(anno.reg$normal~anno.reg$pp_circ)$p.value)
create.boxplot(
        formula = normal ~ pp_circ,
        data =  anno.reg,
        add.stripplot = TRUE,
        points.cex = .3,
        border.col = default.colours(2),
        lwd = 5,
        umb.lwd = 5,
        xlab.label = 'Peptide circRNAs',
        ylab.label = '#normal samples with m6A',
        xaxis.lab = c('No', 'Yes'),
        filename = generate.filename('peptide_circRNA','m6a_nnormal', 'tiff'),
        box.ratio = 0.6,
        ylab.cex = 1.5,
        xlab.cex = 1.5,
        width = 3,
        height = 4,
        style = 'Nature',
        key = list(
                        text = list(
                        lab = pval,
                        cex = 1,
                        col = 'black'
                        ),
                x = 0.25,
                y = 0.95
                ),
);
###
pval <- scientific.notation(wilcox.test(anno.reg$tumor~anno.reg$pp_circ)$p.value);
create.boxplot(
        formula = tumor ~ pp_circ,
        data =  anno.reg,
        add.stripplot = TRUE,
        points.cex = .3,
        border.col = default.colours(2),
        lwd = 5,
        umb.lwd = 5,
        xlab.label = 'Peptide circRNAs',
        ylab.label = '#tumor samples with m6A',
        xaxis.lab = c('No', 'Yes'),
        filename = generate.filename('peptide_circRNA','m6a_ntumor', 'tiff'),
        box.ratio = 0.6,
        ylab.cex = 1.5,
        xlab.cex = 1.5,
        width = 3,
        height = 4,
        style = 'Nature',
        key = list(
                        text = list(
                        lab = pval,
                        cex = 1,
                        col = 'black'
                        ),
                x = 0.25,
                y = 0.95
                ),
);
#### number of m6a associated with peptide circRNAs
to.plot <- circ.peaks;
to.plot[to.plot$npeak>2, ]$npeak <- 2;
npeak.non <- density(to.plot[to.plot$pp_circ==0, ]$npeak, bw = 0.25);
npeak.pp <- density(to.plot[to.plot$pp_circ==1, ]$npeak, bw = 0.25);
col_b <- default.colours(2)
col_f <- unlist(lapply(col_b, function(x) t_col(x, 50)))
pdf(generate.filename('density_npeak', 'circRNA_peptide', 'pdf'));
par(xaxs='i',yaxs='i', ps = 15, font.axis = 1, font.lab= 1)
plot(npeak.non, main = '', xlim = c(0, 3), ylim = c(0, 1.3), col = col_b[1], 
        xlab = expression('#peaks'), axes = FALSE);
axis(1, pos=0, at = seq(0, 3), labels = c(0, 1, 2, 3));
axis(2, pos=0, las = 2);
polygon(npeak.non, col = col_f[1], border = col_b[1]);
polygon(npeak.pp, col = col_f[2], border = col_b[2]);
legend(2.3, 1.2, border = col_b, fill = col_f, bty = 'n', legend = c('No', 'Yes'));
dev.off();
to.plot <- data.frame(table(to.plot$npeak, to.plot$pp_circ));
create.barplot(
	formula = Freq~Var1,
	data = to.plot,
	groups = Var2,
	col = default.colours(2),
	border.col = default.colours(2),
	ylab.label = '# circRNAs',
	xlab.label = '# m6A peaks',
	xaxis.lab = c(1, '>=2'),
	yat = seq(0, 2000, 500),
	ylim = c(0, 2300),
	filename = generate.filename('npeak', 'circRNA_peptide', 'pdf'),
	style = 'Nature',
	key = list(
        points = list(
        	col = default.colours(2),
        	pch = 22,
        	cex = 1.5,
        	fill = default.colours(2)
        	),
        text = list(
        	lab = c('No', 'Yes')
        	),
    x = 0.75,
    y = 0.99
	),
	width = 4
	);
####