library(tidyverse)
library(limma)
library(edgeR)
library(RColorBrewer)
library(ggrepel)

options(ggrepel.max.overlaps = Inf)
extrafont::loadfonts()

# colorblind palette
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Reading data ------------------------------------------------------------

sample_no <- 'Valentine_CD8Tfc_2021' # change accordingly

# directories
workingdir <- paste0(getwd(), '/',sample_no)
plotdir <- 'plots'
datadir <- 'data'
resultsdir <- 'results'

# check for sample directory
ifelse(!dir.exists(file.path(workingdir)), 
       dir.create(file.path(workingdir)), FALSE)

# check for plot directory
ifelse(!dir.exists(file.path(workingdir, plotdir)), 
       dir.create(file.path(workingdir, plotdir)), FALSE)

# check for data directory
ifelse(!dir.exists(file.path(workingdir, datadir)), 
       dir.create(file.path(workingdir, datadir)), FALSE)

# check for results directory
ifelse(!dir.exists(file.path(workingdir, resultsdir)), 
       dir.create(file.path(workingdir, resultsdir)), FALSE)


# read in GSE177911 raw counts
tfc_countdata <- read.csv(file = paste0(workingdir,'/data/','GSE177911_featureCounts_rawcounts_matrix.csv'),
                          row.names = 1,
                          header = TRUE)

tfc_countdata <- tfc_countdata[,c(3,4,1,2)] # reorder for downstream analysis


# read in GSE112540 reanalysis counts
bulkcountdata <- read.csv(file = paste0(workingdir,'/data/','GSE112540_reanalysis_counts.csv'),
                          row.names = 1,
                          header = TRUE)


cols2 <- c("CD8IL2KO_1", "CD8WT_1", "CD8IL2KO_2", "CD8WT_2", 
           "CD8WT_3", "CD8IL2KO_3", "CD8WT_4", "CD8IL2KO_4") 

colnames(bulkcountdata) <- cols2

bulkcountdata <- bulkcountdata[,c(1,3,6,8)] # subset to only contain IL-2-KO samples

countdata <- cbind(bulkcountdata, tfc_countdata)


# QC and Filtering --------------------------------------------------------

# Obtain CPMs
myCPM <- cpm(countdata)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 1

table(rowSums(thresh))

keep <- rowSums(thresh) >= 3

# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep)

plot(myCPM[,1],countdata[,1])

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
abline(v=1)


# Setting up the grouping and factoring for downstream analysis
group <- gsub("_.*", "", colnames(counts.keep)) # clean colnames
group <- factor(group, levels=c("CD8IL2KO","CD8Tfc","CD4Tfh"))

mod.matrix <- model.matrix(~0 + group)
colnames(mod.matrix) <- c("CD8IL2KO","CD8Tfc","CD4Tfh")
mod.matrix

y <- DGEList(counts.keep)

# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)

# barplot for library sizes
y$samples %>%
  ggplot(aes(x = rownames(.), y=lib.size)) +
  geom_bar(stat = "identity") + 
  xlab('') +
  ylab('Library Size') +
  coord_flip() +
  ggtitle("Barplot of library sizes") + 
    theme_bw(base_size = 16)
ggsave(filename = paste0(workingdir,'/plots/',sample_no,"_librarysizes.pdf"), 
       units="in", width=8, height=6, dpi=200)


# check distributions of samples using boxplots
logcounts %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  pivot_longer(!gene, names_to = 'sample', values_to = 'logcpm') %>%
  ggplot(aes(x = sample, y = logcpm)) + 
  geom_boxplot() + 
  xlab('') +
  ylab("logCPM (Unnormalized)") +
  geom_hline(aes(yintercept = median(logcpm)), 
             color="dodgerblue", 
             linetype = 'dashed', 
             size = 2) + 
  theme_bw(base_size = 16)
ggsave(filename = paste0(workingdir,'/plots/',sample_no,"_logcpmboxplots.pdf"), 
       units="in", width=12, height=8, dpi=200)


# plots for filtering lowly expressed genes
samplenames <- colnames(y)
samplenames

L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6

rawlcpm <- cpm(countdata, log = TRUE)

lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(y)
col <- brewer.pal(nsamples, "Paired")

png(paste0(workingdir,'/plots/',sample_no,"_filteringlowgenes.png"), units="in", height=6, width=12, res=300)
par(mfrow=c(1,2))
plot(density(rawlcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(rawlcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
plot(density(logcounts[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(logcounts[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
dev.off()


# Exploratory Analysis ----------------------------------------------------

# Dimension Reduction - MDS
png(paste0(workingdir,'/plots/',sample_no,"_mdsplot.png"), units="in", width=6, height=6, res=300)
plotMDS(y)
dev.off()


# most variable genes heatmap
var_genes <- apply(logcounts, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
highly_variable_lcpm <- logcounts[select_var,]

dim(highly_variable_lcpm)
head(highly_variable_lcpm)

CellType <- as.data.frame(group)
rownames(CellType) <- colnames(highly_variable_lcpm)
colnames(CellType) <- 'CellType'


## Get some nicer colors
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c(cbbPalette)[group]
cn <- group
CellType_2 = list(CellType = c(CD8IL2KO = "#999999", CD8Tfc="#E69F00", CD4Tfh="#56B4E9"))


# Plot the heatmap of top 500 most variable genes 
png(paste0(workingdir,'/plots/',sample_no,"_top500vargenes_heatmap.png"), units="in", width=6, height=6, res=300)
gplots::heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),
                  trace="none", 
                  main="Top 500 most variable genes\n across samples",
                  ColSideColors=col.cell,
                  scale="row", 
                  margins = c(10,5), 
                  colCol = col.cell, 
                  labCol = cn, 
                  labRow = F, 
                  Colv = FALSE,
                  dendrogram = "none")
dev.off()



# Normalization and Transformation ----------------------------------------
y.norm <- calcNormFactors(y, method="TMM") # TMM normalization

par(mfrow=c(1,1))
# voom transformation (with weighting)
png(paste0(workingdir,'/plots/',sample_no,"_voom_transformation.png"), units="in", height=6, width=8, res=300)
vwts <- voomWithQualityWeights(y.norm, mod.matrix, normalize.method = "none", plot = TRUE)
vwts
dev.off()


png(paste0(workingdir,'/plots/',sample_no,"_voom_MDS.png"), units="in", height=6, width=6, res=300)
plotMDS(vwts, col=col.cell)
dev.off()

png(paste0(workingdir,'/plots/',sample_no,"_logcpm_boxplots.png"), units="in", height=6, width=12, res=300)
par(mfrow=c(1,2))
# plot 1 unnormalized
boxplot(logcounts, xlab="", ylab="logCPM",las=2,main="Unnormalised logCPM", cex.axis=0.8)
abline(h=median(logcounts),col="dodgerblue", lwd=3, lty=2)

# plot 2 voom trans logCPM
boxplot(vwts$E, xlab="", ylab="logCPM",las=2,main="Voom transformed logCPM", cex.axis=0.8)
abline(h=median(vwts$E),col="dodgerblue", lwd=3, lty=2)
dev.off()



 # Differential Expression -------------------------------------------------

fit <- lmFit(vwts, mod.matrix)
head(coef(fit))

cont.matrix <- makeContrasts(Diffs = CD8Tfc - CD4Tfh,
                             levels=mod.matrix)
fit.cont <- contrasts.fit(fit, contrasts=cont.matrix)
efit <- eBayes(fit.cont, robust = TRUE)

# Gene Annotations --------------------------------------------------------

gene_annotations <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
                                     keys=rownames(efit),
                                     columns=c("ENTREZID","SYMBOL","GENENAME"), 
                                     keytype = "SYMBOL")

table(gene_annotations$SYMBOL==rownames(efit))

ann_efit <- efit
ann_efit$genes <- gene_annotations

topTable(ann_efit,coef=1,sort.by="p")

cd8tfcvscd4tfh <- topTable(ann_efit,coef = 1,sort.by = "p", n="Inf")

# save the differential expression data
write_csv(cd8tfcvscd4tfh, 
          file = paste0(workingdir,'/results/',sample_no,"_differential_expression.csv"), 
          col_names = TRUE)



expCountsggplot <- function(x){
  # function for quickly looking at expression differences based on a gene
  
  log_gene_counts <- function(x) {
    vwts$E[which(row.names(vwts$E) == x), ]
  }
  
  counts_df <- data.frame(Samples = as.factor(
    gsub('_\\d','', names(log_gene_counts(x)))),
    LogCPM = log_gene_counts(x))
  
  avg <- mean(counts_df$LogCPM)
  
  df_means <- counts_df %>% 
    group_by(Samples) %>% 
    summarise(means = mean(LogCPM))
  
  ggplot(counts_df, aes(Samples, LogCPM, color = Samples)) +
    geom_jitter(size = 3) + 
    geom_point(data = df_means, 
               aes(Samples, means),
               color = "black", 
               shape = 4,
               size = 6,
               alpha = 1) +
    scale_color_brewer(palette = "Dark2") +
    geom_hline(yintercept=avg, 
               color = "dodgerblue", 
               linetype = 4, 
               size = 2) + 
    xlab('') + 
    ylab('logCPM') +
    ggtitle(paste0("Gene:", x)) + 
    theme_bw(base_size = 16)
}


expCountsggplot("Tox2")

expCountsggplot("Cxcr5")


# Heatmaps ----------------------------------------------------------------

# expression counts for plotting in heatmap based on gene names
log_gene_counts <- function(x) {
  vwts$E[which(row.names(vwts$E) == x), ]
}

markergenes <- read.csv(file = paste0(workingdir,'/data/','Valentine_CD8Tfc_2021_heatmap_markergenes.csv'),
                        row.names = 1,
                        header = TRUE)

# marker gene lists
follicular_genes <- markergenes %>%
  filter(Gene.Class == 'Follicular Markers') %>%
  pull(gene)

effector_genes <- markergenes %>%
  filter(Gene.Class == 'Effector Markers') %>%
  pull(gene)


memory_genes <- markergenes %>%
  filter(Gene.Class == 'Memory Markers') %>%
  pull(gene)

exhaustion_markers <- markergenes %>%
  filter(Gene.Class == 'Exhaustion Markers') %>%
  pull(gene)



totalheatmap_genes <-
  length(follicular_genes) + length(effector_genes) +
  length(exhaustion_markers) + length(memory_genes)

ann_row_breaks <-
  c(
    length(follicular_genes),
    length(effector_genes),
    length(exhaustion_markers),
    length(memory_genes)
  )

all_genes <- c(follicular_genes,
               effector_genes,
               exhaustion_markers,
               memory_genes)

all_genes_1 <- all_genes[c(1:21)]
all_genes_2 <- all_genes[c(22:36)]


ann_row = data.frame(`Gene Class` = factor(rep(
  c(
    "Follicular Markers",
    "Effector Markers",
    "Exhaustion Markers",
    "Memory Markers"
  ),
  ann_row_breaks
)))

ann_row_sub1 = data.frame(`Gene Class` = factor(rep(
  c("Follicular Markers",
    "Effector Markers"),
  ann_row_breaks[1:2]
)))

ann_row_sub2 = data.frame(`Gene Class` = factor(rep(
  c("Exhaustion Markers",
    "Memory Markers"),
  ann_row_breaks[3:4]
)))




CellType_3 = list(
  CellType = c(
    CD8IL2KO = "#999999",
    CD8Tfc = "#E69F00",
    CD4Tfh = "#56B4E9"
  ),
  Gene.Class = c(
    `Follicular Markers` = "#009E73",
    `Effector Markers` = "#F0E442",
    `Exhaustion Markers` = "#0072B2",
    `Memory Markers` = "#D55E00"
  )
)

CellType_3_sub1 = list(
  CellType = c(
    CD8IL2KO = "#999999",
    CD8Tfc = "#E69F00",
    CD4Tfh = "#56B4E9"
  ),
  Gene.Class = c(
    `Follicular Markers` = "#009E73",
    `Effector Markers` = "#F0E442"
  )
)

CellType_3_sub2 = list(
  CellType = c(
    CD8IL2KO = "#999999",
    CD8Tfc = "#E69F00",
    CD4Tfh = "#56B4E9"
  ),
  Gene.Class = c(`Exhaustion Markers` = "#0072B2",
                 `Memory Markers` = "#D55E00")
)


gene_types_1 <- data.frame(Gene.Class = ann_row$Gene.Class[1:21])
rownames(gene_types_1) <- make.names(all_genes[1:21], unique = TRUE)
gene_types_2 <- data.frame(Gene.Class = ann_row$Gene.Class[22:36])
rownames(gene_types_2) <- make.names(all_genes[22:36], unique = TRUE)



# generate the heatmaps
group_pheatmap1 <- function(z) {
  
  clust <- t(sapply(z, log_gene_counts))
  pheatmap::pheatmap(
    t(t(scale(t(clust)))),
    annotation_col = gene_types_1,
    annotation_colors = CellType_3_sub1,
    annotation_row = CellType,
    clustering_distance_rows = "correlation",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend = TRUE,
    annotation_legend = TRUE,
    treeheight_row = 0,
    fontsize_col = 14,
    angle_col = "45",
    treeheight_col = 0,
    border_color = "black",
    gaps_row = c(4, 6),
    gaps_col = c(13),
    annotation_names_col = FALSE,
    fontsize = 12,
    cellwidth = 16,
    cellheight = 16
  )
}

group_pheatmap2 <- function(z) {
  
  clust <- t(sapply(z, log_gene_counts))
  pheatmap::pheatmap(
    t(t(scale(t(clust)))),
    annotation_col = gene_types_2,
    annotation_colors = CellType_3_sub2,
    annotation_row = CellType,
    clustering_distance_rows = "correlation",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    legend = TRUE,
    annotation_legend = TRUE,
    treeheight_row = 0,
    fontsize_col = 14,
    angle_col = "45",
    treeheight_col = 0,
    border_color = "black",
    gaps_row = c(4, 6),
    gaps_col = c(9),
    annotation_names_col = FALSE,
    fontsize = 12,
    cellwidth = 16,
    cellheight = 16
  )
}



png(paste0(workingdir,'/plots/',sample_no,"_figure4_heatmap1.png"), units="in", height=6, width=10, res=200)
p_heatmap1 <- group_pheatmap1(all_genes_1)
dev.off()

png(paste0(workingdir,'/plots/',sample_no,"_figure4_heatmap2.png"), units="in", height=6, width=10, res=200)
p_heatmap2 <- group_pheatmap2(all_genes_2)
dev.off()


# Volcano Plots -----------------------------------------------------------

cutoff_up <- sort(cd8tfcvscd4tfh$adj.P.Val)[107]
cutoff_down <- sort(cd8tfcvscd4tfh$adj.P.Val)[24]

deg_dat <- cd8tfcvscd4tfh %>% 
  mutate(TopGeneLabel_up=ifelse(logFC >= 2 & adj.P.Val <= cutoff_up, SYMBOL, "")) %>%
  mutate(TopGeneLabel_down=ifelse(logFC <= -2 & adj.P.Val <= cutoff_down, SYMBOL, ""))

# Here we see how many genes labeled
length(unique(deg_dat$TopGeneLabel_up))
length(unique(deg_dat$TopGeneLabel_down))

deg_dat <- deg_dat %>%
  mutate(`-log10(adj.P.Val)` = -log10(adj.P.Val))


cbbPalette <-
  c(
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  )



gg_volc_pub <- function(z){
  ggplot(z) +
    geom_point(mapping = aes(x = logFC, y=`-log10(adj.P.Val)`, 
                             color="LFC<2, p>0.05"),
               size = 1) +
    geom_point(subset(z, adj.P.Val < .05), 
               mapping = aes(x = logFC, y=`-log10(adj.P.Val)`, 
                             color = "LFC<2, p<0.05"),
               size = 1) +
    geom_point(subset(z, abs(logFC)>2), 
               mapping = aes(x = logFC, y=`-log10(adj.P.Val)`, 
                             color = "LFC>2, p>0.05"),
               size = 1) +
    geom_point(subset(z, adj.P.Val<.05 & abs(logFC)>2), 
               mapping = aes(x = logFC, y=`-log10(adj.P.Val)`, 
                             color = "LFC>2, p<0.05"),
               size = 1) +
    geom_text_repel(aes(x = logFC, y=`-log10(adj.P.Val)`), 
                    label=z$TopGeneLabel_up, 
                    force = 3, size = 3, 
                    segment.size = 1, 
                    fontface = "bold", 
                    alpha = 3/4,
                    family = "Arial",
                    max.time = 10,
                    max.iter = Inf,
                    direction = "both",
                    verbose = TRUE,
                    seed = 2021) +
    geom_text_repel(aes(x = logFC, y=`-log10(adj.P.Val)`), 
                    label=z$TopGeneLabel_down, 
                    force = 4, size = 3,
                    segment.size = 1, 
                    fontface="bold", 
                    alpha = 3/4,
                    family = "Arial",
                    max.time = 10,
                    max.iter = Inf,
                    direction = "both",
                    verbose = TRUE,
                    seed = 2021) +
    theme_classic(base_size = 16, base_family = "Arial")
}

p1 <- gg_volc_pub(deg_dat)
p1.done <- p1 + scale_color_manual(name = "", values = cbbPalette) +
  labs(title = "",
       subtitle = "") +
  xlab(expression("log"[2] * "FC")) +
  ylab(expression("-log"[10] * "(adj.P.Val)")) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    legend.background = element_rect(
      color = "black",
      size = 0.5,
      linetype = "solid"
    ),
    text = element_text(family = "Arial", color = "black"),
    axis.text = element_text(family = "Arial", color = "black")
  )

p1.done

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"_figure4_volcanoplot.png"), 
       units="in", height=8, width=9, dpi=200)



# Gene Ontology -----------------------------------------------------------
library(clusterProfiler)

cd8tfcvscd4tfh <- cd8tfcvscd4tfh %>%
  filter(!is.na(adj.P.Val)) %>% 
  mutate(`-log10(adj.P.Val)` = -log10(adj.P.Val))

geneList.diffs <- cd8tfcvscd4tfh[,4]
names(geneList.diffs) <- as.character(cd8tfcvscd4tfh[,2])
geneList.diffs <- sort(geneList.diffs, decreasing = TRUE)


universe.diffs <- cd8tfcvscd4tfh %>% 
  pull(ENTREZID)
sigGenes.diffs <- cd8tfcvscd4tfh %>%
  filter(adj.P.Val < 0.05,!is.na(ENTREZID)) %>% 
  pull(ENTREZID)


# running enrichGO
enrich_go.diffs <- enrichGO(
  gene = sigGenes.diffs,
  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  universe = names(geneList.diffs),
  qvalueCutoff = 0.05,
  readable = TRUE
)

# generate goplot
goplot <- dotplot(enrich_go.diffs, showCategory = 10) + 
  scale_size(range = c(7,14)) +
  theme_bw(base_size = 16, base_family = "Arial") 

goplot + theme(text = element_text(family = "Arial", colour = "black"),
            axis.text = element_text(family = "Arial", color = "black"))


ggsave(filename = paste0(workingdir,'/plots/',sample_no, "_figure4_goplot.png"), 
       dpi = 200, 
       scale = 1.55)



