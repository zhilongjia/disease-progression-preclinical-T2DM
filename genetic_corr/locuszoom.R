library(locuszoomr)
rm(list=ls())
library(EnsDb.Hsapiens.v75)


gwas_dm <- read.csv('data/gwas_summary/t2dm/res_subtype1_pc10.txt.gz', sep = ' ')


gwas_dm <- subset(gwas_dm, select = c('CHR', 'SNP', 'BP', 'P', 'BETA', 'SE'))

colnames(gwas_dm) <- c('chrom', 'rsid', 'pos', 'p', 'beta', 'se')

distance_ld <- 250000


if (require(EnsDb.Hsapiens.v75)) {
  loc <- locus(gwas_dm, gene = 'TRANK1',
               flank = c(distance_ld, distance_ld), 
               ens_db = "EnsDb.Hsapiens.v75")
  loc <- link_LD(loc, token = "19b0a7d7cbe2", pop = "EUR")
}

gwas_dis <- read.csv('data/gwas_summary/disease/BD.txt.gz', sep = ' ')
gwas_dis <- subset(gwas_dis, select = c('CHR', 'SNP', 'BP', 'P', 'BETA', 'SE'))
colnames(gwas_dis) <- c('chrom', 'rsid', 'pos', 'p', 'beta', 'se')

if (require(EnsDb.Hsapiens.v75)) {
  loc2 <- locus(gwas_dis, gene = 'TRANK1',
               flank = c(distance_ld, distance_ld), 
               ens_db = "EnsDb.Hsapiens.v75")
  loc2 <- link_LD(loc2, token = "19b0a7d7cbe2", pop = "EUR")
}
source('scatter_plot.R')
# save to jpeg
jpeg("results/figs/coloc_s1_BD2.jpg", width = 6, height = 8, units = "in", res = 600)
oldpar <- set_layers(2)

# region plot for Subtype 1
scatter_plot(loc, xticks = FALSE, labels = c('rs9834970'), 
             label_y = c(25), col = "grey", 
             legend_pos = "topright",lwd = 0.8,
             # font size
             cex = 1.6, cex.axis = 1.8, cex.lab = 1.8, cex.legend = 1.5)

# region plot for BD
scatter_plot(loc2, col = "orange", xticks = FALSE, labels = c('rs9834970'), 
             legend_pos = "topright", lwd = 0.5,
             cex = 1.6, cex.axis = 1.8, cex.lab = 1.8, cex.legend = 1.5)

# gene track
genetracks(loc, maxrows = 3, filter_gene_biotype = 'protein_coding',
           gene_col = 'grey', exon_col = 'red', exon_border = 'darkgrey',
           cex.axis = 1.8, cex.text = 1.6, cex.lab = 1.8
           )

par(oldpar)  # revert par() settings
dev.off()
