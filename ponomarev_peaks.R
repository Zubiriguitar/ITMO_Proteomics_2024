library(readxl)
library(limma)
library(ape)
library(dendextend)
library(RColorBrewer)
library(pvclust)
library(gplots)
library(NMF)
library(vegan)
library(Biobase)
library(DT)
library(ggplot2)
library(impute)
library(ggrepel)

# Data uploading and clearing.


peaks <- read.csv('peaks_data.csv', sep=',')
peaks <- peaks[-c(849,1322,3431),]
peaks <- peaks[-(which(peaks$Gene_id %in% '')),]
peaks <- peaks[-(which(peaks$Gene_id %in% 'nan')),]
rownames(peaks) <- peaks$Gene_id
peaks <- peaks[,-1]
peaks <- peaks[,-1]
peaks <- peaks[,-1]

## EDA

#Remove genes with half and more missing values

genes_with_NA <- names(which(rowSums(is.na(peaks)) > ncol(peaks)/2))
peaks <- peaks[!rownames(peaks) %in% genes_with_NA,]

peaks <- as.matrix(peaks)
peaks_trans <- t(peaks)
knn_peaks <- impute.knn(peaks_trans, k = 5)
knn_peaks <- knn_peaks$data
knn_peaks <- as.data.frame(knn_peaks)
knn_peaks <- t(as.matrix(knn_peaks))


peaks_experiment <- as.factor(c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2"))

pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[peaks_experiment]
boxplot(knn_peaks, outline = FALSE, main = "Initial data", col = cols)
legend("topright", levels(peaks_experiment), fill = pal, bty = "n", xpd = T)

# Normalization by log and quantile
peaks_log <- log2(knn_peaks)

peaks_log <- peaks_log[is.finite(rowSums(peaks_log)),]



peaks_norm <- normalizeQuantiles(as.matrix(peaks_log))
boxplot(peaks_norm, outline = FALSE, main = "normalized data", col = cols)
legend("topright", levels(peaks_experiment), fill = pal, bty = "n", xpd = T)

#PCA
peaks_pca <- t(peaks_norm)
terapod_pca <- rda(peaks_pca, scale = TRUE)

rownames(peaks_pca) <- c("BT","BT","BT","BT","BT","BT","BT","BT","BT","BT","BT","BT","BT","BT","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK","CJK")

df_scores <- data.frame(peaks_pca,
                        scores(terapod_pca, display = "sites", choices = c(1, 2, 3), scaling = "sites"))

p_scores <- ggplot(df_scores, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = rownames(peaks_pca)), alpha = 0.5) +
  coord_equal(xlim = c(-3, 3), ylim = c(-3, 3)) + ggtitle(label = "Ordination") + theme_bw()
p_scores

#MA-plot 


maplot <- function(X1, X2, pch = 21, main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2), ...){
  # Coordinates
  X <- (rowMeans(X2) + rowMeans(X1)) / 2
  Y <- rowMeans(X2) - rowMeans(X1)
  # Graph
  scatter.smooth(x = X, y = Y, main = main, pch = pch, xlab = xlab, ylab = ylab, lpars = lpars, ...)
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
}



# BT vs CJK


maplot(peaks_norm[,1:14], peaks_norm[,15:33])

# Differential expression 
expr_data <- as.matrix(peaks_norm)

# Data 
pheno_data <- data.frame(peaks_experiment)
rownames(pheno_data) <- colnames(peaks_norm)
pheno_metadata <- data.frame(
  labelDescription = c("Experimental condition"),
  row.names=c("Condition"))
pheno_data <- new("AnnotatedDataFrame",
                  data = pheno_data,
                  varMetadata = pheno_metadata)

# Protein data
feature_data <- data.frame(Prot = rownames(expr_data))
rownames(feature_data) <- rownames(expr_data)
feature_metadata <- data.frame(
  labelDescription = c("Protain name"),
  row.names = c("Protain"))
f_data <- new("AnnotatedDataFrame",
              data = feature_data,
              varMetadata = feature_metadata)


exp_set <-
  ExpressionSet(assayData = expr_data,
                phenoData = pheno_data,
                featureData = f_data)


X <- model.matrix(~ peaks_experiment, pData(exp_set))
fit <- lmFit(exp_set, design = X, method = "robust", maxit = 1000)
efit <- eBayes(fit)

MA_limma <- function(efit, coef, n = 10, signif = TRUE, fdr = 0.05, lfc = 0, text = TRUE, cex.text = 0.8, col.text = "grey20", main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", pch = 19, pch.signif = 21, col = "darkgreen", alpha = 0.3, cex = 0.3, ...){

  R <- efit$coefficients[, coef]
  I <- efit$Amean

  col_btransp <- adjustcolor(col, alpha.f = alpha)
  # Graph
  plot(I, R, cex = cex, main = main, pch = pch, xlab = xlab, ylab = ylab, col = col_btransp, ...)
  abline(h = 0)
  # labeling differentially expressed genes
  if(signif){
    sign <- p.adjust(efit$p.value[, coef], method = "BH") <= fdr
    large <- abs(efit$coefficients[, coef]) >= lfc
    points(I[sign & large], R[sign & large], cex = cex*2, col = "orange2", pch = pch.signif)
  }
  # labeling proteins with higher expression
  if(text){
    ord <- order(efit$lods[, coef], decreasing = TRUE)
    top_n <- ord[1:n]
    text(I[top_n], R[top_n], labels = efit$genes[top_n, ], pos = 4, cex = cex.text, col = col.text)
  }
}


MA_limma(efit, coef = 2, n = 30)

# Exploring DE proteins 
# Dirst 20 DE proteins
my_list <- topTable(efit, coef = 2, n = 100)
# Filtering ExpressionSet
dif_exp_set <- exp_set[fData(exp_set)$Prot %in% my_list$Prot, ]

dat <- as.matrix(exprs(dif_exp_set))

#Heatmap
pal_blue_red <- colorpanel(75, low = "steelblue", mid = "black", high = "red")
heatmap.2(dat, col = pal_blue_red, scale = "row", key = TRUE, symkey = FALSE, density.info = "none", trace = "none", cexRow = 0.9, cexCol = 1, margins = c(4, 3), keysize = 0.8, key.par = list(mar = c(3, 0.1, 3, 0.1)))


topTable(efit, coef = 2)
numGenes <- nrow(exprs(exp_set))

full_list <- topTable(efit, number = numGenes)
full_list <- full_list[full_list$adj.P.Val <= 0.05,]
write.csv(full_list, 'DE_full')

MA_limma(efit, coef = 2, n = 4)

# MA-plot of first 20 DE proteins with FC > 2
MA_limma(efit, coef = 2, n = 80, text = F, lfc = 1)

my_list <- full_list
dif_exp_set <- exp_set[fData(exp_set)$Prot %in% my_list$Prot, ]



## Tree of DE proteins

diff_prot <- rownames(full_list)
diff_expressed <- as.data.frame(peaks_norm)[diff_prot,]
t_diff_expressed <- t(diff_expressed)
#rownames(t_diff_expressed) <-  as.factor(gsub("_[^_]*", replacement = "", rownames(t_diff_expressed)))
#rownames(t_diff_expressed) <- make.unique(as.character(pheno$yeast_experiment))

d <- dist(x = t_diff_expressed, method = "canberra")

mouse_hc_avg <- hclust(d, method = "average")
mouse_ph_avg <- as.phylo(mouse_hc_avg)
mouse_den_avg <- as.dendrogram(mouse_hc_avg)

get_colours <- function(dend, n_chars, palette = "Dark2"){ #nchars = первые нескольок симовлов которые будут использоваться для разделения фактора на группы
  labs <- get_leaves_attr(dend, "label")
  group <- substr(labs, start = 0, stop = n_chars)
  group <- factor(group)
  cols <- brewer.pal(length(levels(group)), name = palette)[group]
  return(cols)
}

cols <- get_colours(dend = mouse_den_avg, n_chars = 6)
den_avg_c <- color_labels(dend = mouse_den_avg, col = cols)
plot(den_avg_c, horiz = TRUE)

# BT vs CJK


maplot(peaks_norm[,1:14], peaks_norm[,15:33])


# MA-plot первых 20 дифференциально экспрессируемых белков, но таких, чтобы уровень экспрессии различался в 2 раза
MA_limma(efit, coef = 2, n = 20, text = F, lfc = 1)




# Volcano plot 

#First 20 DE proteins
my_list <- full_list
# Фильтруем ExpressionSet
dif_exp_set <- exp_set[fData(exp_set)$Prot %in% my_list$Prot, ]


volcano_list <- full_list

volcano1 <- ggplot(data = volcano_list, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point()

volcano2 <- volcano1 + geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

volcano_list$diffexpressed <- "NO"       


volcano_list$diffexpressed[volcano_list$logFC > 1 & volcano_list$adj.P.Val < 0.05] <- "UP"

volcano_list$diffexpressed[volcano_list$logFC < -1 & volcano_list$adj.P.Val< 0.05] <- "DOWN"

volcanodif1 <- ggplot(data = volcano_list, aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed)) + geom_point() + theme_minimal()


volcanodif2 <- volcanodif1 + geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

volcano_list$Prot <- as.character(volcano_list$Prot)
volcano_list$delabel <- NA
volcano_list$delabel[volcano_list$diffexpressed != "NO"] <- volcano_list$Prot[volcano_list$diffexpressed != "NO"]
volcano_list[volcano_list$adj.P.Val<=0.6,]$delabel <- volcano_list[volcano_list$adj.P.Val<=0.6,]$Prot

plot_proteins <- ggplot(data=volcano_list, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size=3, colour = 'black', max.overlaps = 30)

plot_final <- plot_proteins + geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(0.05), col="red")
plot_final

# GO enrichment analysis



# keep only the significant proteins results
sig <- subset(full_list, adj.P.Val < 0.05)
# get the significant up-regulated proteins
up <- subset(sig, logFC > 0)
# get the significant down-regulated proteins
down <- subset(sig, logFC < 0)

library(gprofiler2)

# needed to convert to enrichResult object
up_names <- gconvert(row.names(up))
down_names <- gconvert(row.names(down))


## Up-regulated proteins

# enrichment analysis using proteins names
multi_gp_up_reg <- gost(list("up-regulated" = up_names$name), multi_query = FALSE, evcodes =TRUE)
# modify the g:Profiler data frame
gp_mod_up = multi_gp_up_reg$result[, c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
gp_mod_up <- gp_mod_up[order(gp_mod_up$p_value), ]
gp_mod_up_BP <- gp_mod_up[gp_mod_up$source == "GO:BP", ]
gp_mod_up_BP$GeneRatio <- paste0(gp_mod_up_BP$intersection_size,  "/", gp_mod_up_BP$query_size)
gp_mod_up_BP$BgRatio <- paste0(gp_mod_up_BP$term_size, "/", gp_mod_up_BP$effective_domain_size)
names(gp_mod_up_BP) <- c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
gp_mod_up_BP$geneID <- gsub(",", "/", gp_mod_up_BP$geneID)
row.names(gp_mod_up_BP) <- gp_mod_up_BP$ID
gp_mod_enrich_up_BP <- new("enrichResult", result = gp_mod_up_BP)


#Draw enrichment plot:

enrichplot::dotplot(gp_mod_enrich_up_BP, showCategory = 10) + ggplot2::labs(title = "up-regulated") + ggplot2::scale_color_gradient(low = "lightseagreen", high = "darkorange1")


## Down-regulated proteins

# enrichment analysis using gene names
multi_gp_down_reg <- gost(list("down-regulated" = down_names$name), multi_query = FALSE, evcodes =TRUE)
# modify the g:Profiler data frame
gp_mod_down = multi_gp_down_reg$result[, c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
gp_mod_down <- gp_mod_down[order(gp_mod_down$p_value), ]
# BP
gp_mod_down_BP <- gp_mod_down[gp_mod_down$source == "GO:BP", ]
gp_mod_down_BP$GeneRatio <- paste0(gp_mod_down_BP$intersection_size,  "/", gp_mod_down_BP$query_size)
gp_mod_down_BP$BgRatio <-  paste0(gp_mod_down_BP$term_size, "/", gp_mod_down_BP$effective_domain_size)
names(gp_mod_down_BP) <- c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
gp_mod_down_BP$geneID <- gsub(",", "/", gp_mod_down_BP$geneID)
gp_mod_enrich_down_BP <- new("enrichResult", result = gp_mod_down_BP)


#Draw enrichment plot:

enrichplot::dotplot(gp_mod_enrich_down_BP, showCategory = 10) + ggplot2::labs(title = "down-regulated") + ggplot2::scale_color_gradient(low = "lightseagreen", high = "darkorange1")


#The enrichment of these GO terms in both up and down regulated  suggests that 
#there are significant alterations in cellular processes related to stress response,
#cellular organization, cytoskeleton dynamics, and biomolecular localization within valve tissues undergoing calcification.
#These changes may contribute to the pathogenesis of valve calcification and could represent potential targets
#for therapeutic intervention or diagnostic biomarkers.