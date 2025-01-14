library(tidyverse)
library(GEOquery)
library(Biobase)
library(biomaRt)
library(org.Hs.eg.db)
library(fgsea)
library(GSEABase)
library(data.table)
library(readr)


# input_file <- "GSE64810_series_matrix.txt.gz"
# output_file <- "GSE64810_series_matrix.txt.tsv"
# df <- read_delim(input_file, 
#                  delim = "\t",    
#                  col_names = TRUE) # 如果文件有表头，设为TRUE/#If the file has a header, set it to TRUE
# 
#将读入的dataframe写出为.tsv文/Write as. tsv file
# write_tsv(df, output_file)


#获取并处理GEO数据/#Obtain and process GEO data
gse <- GEOquery::getGEO('GSE64810', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
filtered_metadata <- dplyr::select(metadata, 1:2, 6, 8:9, 54:65)

#重命名列，导出样本元数据/Export sample metadata
colnames(filtered_metadata) <- c('sample_name', 'geo_accession', 'type', 'source_name', 'organism',
                                 'age_at_death', 'age_at_onset', 'cag', 'diagnosis', 'duration', 'hv_cortical_score',
                                 'hv_striatal_score', 'mrna_seq_reads', 'pmi', 'rin', 'tissue', 'vonsattel_grade')
write.csv(filtered_metadata, 'data/sample_info.csv')

# 定义run_gsea函数：对DESeq2结果进行基因集富集分析/Define run_gsea function: Perform gene set enrichment analysis on DESeq2 results
run_gsea <- function(deseq_res, min_size, max_size) {
# 提取、清理gene ID并转换为SYMBOL/Extract, clean up gene IDs and convert them to SYMBOL
  deseq_res <- dplyr::mutate(deseq_res, 
                             genes= str_extract(genes, '.*\\.') %>% str_replace('\\.', ''))
  
  get_hgnc <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=deseq_res$genes, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
  
  get_hgnc <- as_tibble(get_hgnc)
  
# 合并DE结果与符号注释，构建排名列表/Merge DE results with symbol annotations to construct a ranking list
  hgnc_res <- inner_join(deseq_res, get_hgnc, by=c("genes"="ENSEMBL"))
  
  rank_list <- hgnc_res %>% 
    dplyr::select(SYMBOL, log2FoldChange) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(log2FoldChange)) %>% 
    deframe()
  
# 加载GMT基因集并运行fgsea/Load GMT geneset and run fgsea
  pathways <- gmtPathways("data/c2.cp.v7.5.1.symbols.gmt")
  
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = rank_list,
                    minSize  = min_size,
                    maxSize  = max_size)
  
  return(fgseaRes)
}

# 读取DESeq2差异表达结果，运行GSEA分析并导出结果/Read the differential expression results of DESeq2, run GSEA analysis and export the results
deseqRes <- read_delim('data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.tsv',
                       delim = '\t', 
                       col_names = TRUE) %>%
  dplyr::rename('genes'='...1')

fgseaRes <- run_gsea(deseqRes, 
                     min_size = 15, 
                     max_size = 500)

fgseaRes <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

write_csv(fgseaRes, 'data/fgsea_results.csv')
