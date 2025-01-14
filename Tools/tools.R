# 加载必要库
library(tidyverse)
library(shiny)
library(ggplot2)
library(bslib)
library(glue)
library(colourpicker)
library(patchwork)
library(gplots)
library(pheatmap)
library(ggbeeswarm)
library(DT)
library(fgsea)
library(GSEABase)

# ===========================
# Metadata Analysis Functions
# ===========================

# 绘制元数据表摘要
draw_table_Summary <- function(metadata) {
  col1 <- c('Age', 'Diagnosis', 'Postmortem interval (hours)', 'Onset Age', 'Disease Duration (years)')
  col2 <- c('double', 'factor', 'double', 'double', 'double')
  
  row_age <- paste(round(mean(metadata$age_of_death, na.rm = TRUE), 2), 
                   "(+/-", round(sd(metadata$age_of_death, na.rm = TRUE), 2), ")")
  row_diag <- "Normal, Huntington's Disease"
  row_pmi <- paste(round(mean(metadata$pmi, na.rm = TRUE), 2), 
                   "(+/-", round(sd(metadata$pmi, na.rm = TRUE), 2), ")")
  row_onset <- paste(round(mean(metadata$age_of_onset, na.rm = TRUE), 2), 
                     "(+/-", round(sd(metadata$age_of_onset, na.rm = TRUE), 2), ")")
  row_dd <- paste(round(mean(metadata$Duration, na.rm = TRUE), 2), 
                  "(+/-", round(sd(metadata$Duration, na.rm = TRUE), 2), ")")
  
  col3 <- c(row_age, row_diag, row_pmi, row_onset, row_dd)
  
  tibble(
    'Column Name' = col1,
    'Type' = col2,
    'Mean (SD) or Distinct Values' = col3
  )
}

# 绘制元数据表
draw_table_metadata <- function(metadata) {
  data <- metadata %>%
    select(Run, age_of_death, Diagnosis, mrna.seq_reads, pmi, rin)
  
  DT::datatable(
    data, 
    extensions = 'Buttons', 
    class = "display",
    options = list(
      paging = TRUE, searching = TRUE, fixedColumns = TRUE,
      autoWidth = TRUE, ordering = TRUE, dom = 'Bfrtip',
      buttons = c('copy', 'csv')
    )
  )
}

# 绘制小提琴图
draw_violin_metadata <- function(metadata, y_axis) {
  ggplot(metadata, aes(x = Diagnosis, y = !!sym(y_axis), fill = Diagnosis)) +
    geom_violin(na.rm = TRUE) +
    ylab(y_axis)
}

# ===========================
# Counts Matrix Functions
# ===========================

# 判断是否通过过滤
ifpassed <- function(p1, p2) {
  p1 && p2
}

# 绘制Counts矩阵的摘要表
draw_summary_matrix <- function(matrix) {
  gene_num <- nrow(matrix)
  sample_num <- ncol(matrix) - 7
  pass_num <- sum(matrix$passed == TRUE)
  pass_per <- round(pass_num / gene_num * 100, 2)
  nonpass_num <- gene_num - pass_num
  nonpass_per <- round(nonpass_num / gene_num * 100, 2)
  
  tibble(
    Category = c("Number of Genes", "Number of Samples", 
                 "Genes Passing Filter", "Genes Not Passing Filter"),
    Number_Percentage = c(
      gene_num, sample_num, 
      paste0(pass_num, " (", pass_per, "%)"), 
      paste0(nonpass_num, " (", nonpass_per, "%)")
    )
  )
}

# 绘制Counts散点图
draw_scatter_matrix <- function(matrix) {
  plot1 <- ggplot(matrix, aes(x = mid, y = var, col = passed)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10(oob = scales::squish_infinite)
  
  plot2 <- ggplot(matrix, aes(x = mid, y = non_zero, col = passed)) +
    geom_point() +
    scale_x_log10(oob = scales::squish_infinite)
  
  plot1 + plot2
}

# 绘制Counts热图
draw_heatmap_matrix <- function(matrix, metadata) {
  data <- matrix %>% filter(passed == TRUE)
  data <- data[, 2:(ncol(data) - 6)] %>% as.matrix()
  
  color <- if_else(metadata$Diagnosis == "Huntington's Disease", 'red', 'blue')
  
  heatmap.2(
    data, ColSideColors = color, scale = 'row',
    symkey = FALSE, density.info = "none", trace = "none",
    ylab = "Genes", xlab = "Samples", main = "Counts Heatmap",
    labRow = FALSE
  )
}

# PCA绘图
draw_PCA_matrix <- function(matrix, meta, PC1, PC2, topPC) {
  data <- matrix %>% filter(passed == TRUE)
  data_PCA <- data[, 2:(ncol(data) - 6)] %>% as.matrix()
  rownames(data_PCA) <- data$X
  
  pca <- prcomp(t(data_PCA), center = FALSE, scale = TRUE)
  var_p <- pca$sdev^2 / sum(pca$sdev^2)
  
  plot1 <- pca$x[, 1:topPC] %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = "PC", values_to = "projection") %>%
    mutate(PC = fct_relevel(PC, str_c("PC", 1:20))) %>%
    ggplot(aes(x = PC, y = projection)) +
    geom_beeswarm() + 
    labs(title = "PCA Projection Plot") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
  
  plot2 <- tibble(
    type = meta$Diagnosis,
    PC1 = pca$x[, PC1],
    PC2 = pca$x[, PC2]
  ) %>%
    ggplot(aes(x = PC1, y = PC2, col = type)) +
    geom_point() +
    labs(
      title = "Counts PCA",
      x = paste0("PC", PC1, " (", round(var_p[PC1] * 100), "% variance)"),
      y = paste0("PC", PC2, " (", round(var_p[PC2] * 100), "% variance)")
    )
  
  plot1 / plot2
}

# ===========================
# Differential Expression (DE) Functions
# ===========================

# DE表格
draw_table_DE <- function(dataf, slider) {
  data <- dataf %>%
    filter(log10(padj) < slider) %>%
    mutate(
      pvalue = formatC(pvalue, format = 'e'),
      padj = formatC(padj, format = 'e')
    )
  
  DT::datatable(
    data, 
    extensions = 'Buttons', 
    class = "display",
    options = list(
      paging = TRUE, searching = TRUE, fixedColumns = TRUE,
      autoWidth = TRUE, ordering = TRUE, dom = 'Bfrtip',
      buttons = c('copy', 'csv')
    )
  )
}

# Volcano Plot
volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
  data <- dataf %>%
    filter(!is.na(padj)) %>%
    mutate(colorf = log10(padj) < slider)
  
  ggplot(data, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)), color = colorf)) +
    geom_point() +
    scale_color_manual(values = c(color1, color2)) +
    labs(color = glue("p-value < 10^{slider}"))
}

# ===========================
# GSEA Functions
# ===========================

# 绘制Top GSEA结果
draw_top_gsea <- function(fgsea_result, num) {
  fgsea_result <- fgsea_result %>%
    mutate(pathway = forcats::fct_reorder(pathway, NES))
  
  ggplot(fgsea_result[1:num, ]) +
    geom_bar(aes(x = pathway, y = NES, fill = NES > 0), stat = 'identity') +
    scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
    coord_flip() +
    labs(title = "Top GSEA Results", y = "NES", x = "")
}

# GSEA表格
draw_table_gsea <- function(fgsea_result, pval_gsea, thre_gsea) {
  range <- switch(thre_gsea,
                  "Positive" = c(0, Inf),
                  "Negative" = c(-Inf, 0),
                  "All" = c(-Inf, Inf)
  )
  
  data <- fgsea_result %>%
    filter(log10(padj) < pval_gsea, NES > range[1], NES < range[2])
  
  DT::datatable(
    data, 
    extensions = 'Buttons', 
    class = "display",
    options = list(
      paging = TRUE, searching = TRUE, fixedColumns = TRUE,
      autoWidth = TRUE, ordering = TRUE, dom = 'Bfrtip',
      buttons = c('copy', 'csv')
    )
  )
}

# GSEA散点图
draw_plot2_gsea <- function(fgsea_result, pval_plot_gsea) {
  fgsea_result <- fgsea_result %>%
    mutate(passed = log10(padj) < pval_plot_gsea)
  
  ggplot(fgsea_result, aes(x = NES, y = padj, col = passed)) +
    geom_point() +
    scale_y_log10()
}
