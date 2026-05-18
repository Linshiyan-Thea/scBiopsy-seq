library(Rtsne)
library(mclust)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(pheatmap)
library(ggplotify)

mode <- "count_clean"
data <- read.table(paste0("DGE_", mode, ".txt"))
data <- data + 1
sample <- colnames(data)
cell <- do.call(rbind, strsplit(sample, split="\\."))[,1]
cell <- unique(cell)

for (c in cell){
  c0 <- paste0(c, ".0")
  c1 <- paste0(c, ".1")
  c4 <- paste0(c, ".4")
  check <- (c0 %in% sample) + (c1 %in% sample) + (c4 %in% sample)
  if (check == 3){
    e0 <- data[, c0]
    e1 <- data[, c1]
    e4 <- data[, c4]
    fc01 <- log2(e1 / e0)
    fc04 <- log2(e4 / e0)
    fc14 <- log2(e4 / e1)
    fc01 <- data.frame(fc01)
    fc04 <- data.frame(fc04)
    fc14 <- data.frame(fc14)
    colnames(fc01) <- c
    colnames(fc04) <- c
    colnames(fc14) <- c
    if (c == cell[1]){
      rownames(fc01) <- paste0('t01.', rownames(data))
      rownames(fc04) <- paste0('t04.', rownames(data))
      rownames(fc14) <- paste0('t14.', rownames(data))
      FC01 <- data.frame(fc01)
      FC04 <- data.frame(fc04)
      FC14 <- data.frame(fc14)
    } else {
      FC01 <- cbind(FC01, fc01)
      FC04 <- cbind(FC04, fc04)
      FC14 <- cbind(FC14, fc14)
    }
  }
}

tsne_plot <- function(FC, time='0-1-4', perplexity=10){
  FC <- data.frame(t(na.omit(t(FC))))
  FC <- FC[which(rowSums(FC) > 0), ]
  FC_ <- FC
  FC_[FC_ != 0] <- 1
  FC <- FC[which(rowSums(FC_) > dim(FC_)[2] * 0.9), ]
  
  pca <- prcomp(t(FC), scale. = FALSE)
  rd1 <- pca$x[,1:2]
  rd2 <- Rtsne(t(FC), perplexity = perplexity)$Y
  rd2 <- data.frame(rd2)
  colnames(rd2) <- c('tSNE1', 'tSNE2')
  rd2$Sample <- colnames(FC)
  cl2 <- kmeans(rd1, centers = 2)$cluster
  rd2$k2 <- as.character(cl2)
  rd2$time <- time
  
  p <- ggplot(rd2, aes(x = tSNE1, y = tSNE2, fill = time)) + geom_point() + theme_classic()
  ps <- ggplot(rd2, aes(x = tSNE1, y = tSNE2, fill = time)) +
          geom_point() +
          geom_text_repel(aes(label = Sample), size = 2, box.padding = 0.1, min.segment.length = 0.15, show.legend = FALSE) +
          theme_classic()
  pk2 <- ggplot(rd2, aes(x = tSNE1, y = tSNE2, color = k2)) +
           geom_point() +
           geom_text_repel(aes(label = Sample), size = 2, box.padding = 0.1, min.segment.length = 0.15, show.legend = FALSE) +
           theme_classic()
  
  pdf(paste0('logFC.', mode, '.', time, '.pdf'), width = 14, height = 3)
  print(plot_grid(p, ps, pk2, ncol = 3))
  dev.off()
  
  return(rd2)
}

FC <- rbind(FC01, FC14)
rd2 <- tsne_plot(FC, time = '0-1-4', perplexity = 4)

FC <- data.frame(t(na.omit(t(FC))))
FC <- FC[which(rowSums(FC) > 0), ]
FC_ <- FC
FC_[FC_ != 0] <- 1
FC <- FC[which(rowSums(FC_) > dim(FC_)[2] * 0.9), ]
write.table(FC, paste0('logFC.', mode, '.0-1-4.txt'), sep = '\t')

data <- FC
anno <- data.frame(cell = rd2$Sample, clus = paste0('clu', rd2$k2))

d0 <- data[, anno$cell[which(anno$clus == 'clu1')]]
d1 <- data[, anno$cell[which(anno$clus == 'clu2')]]
pv <- c()
for (i in 1:dim(data)[1]){
  a <- wilcox.test(as.numeric(d0[i,]), as.numeric(d1[i,]))
  pv <- c(pv, a$p.value)
}

de_info <- data.frame(row.names = rownames(data), pvalue = pv)
de_info$is_DE <- de_info$pvalue < 0.05
de_info <- de_info[order(de_info$pvalue), ]
de <- rownames(de_info)[which(de_info$is_DE == TRUE)]
write.table(de_info, paste0('pvalue.', mode, '.0-1-4.txt'), sep = '\t')
write(de, paste0('DE.', mode, '.0-1-4.txt'), sep = '\t')
write.table(data[de,], paste0('logFC_DE.', mode, '.0-1-4.txt'), sep = '\t')
write.table(de_info[de,], paste0('pvalue_DE.', mode, '.0-1-4.txt'), sep = '\t')
