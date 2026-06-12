#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
suppressMessages(library(argparser, quietly = TRUE))

# Create a parser
p <- arg_parser("run KEGG enrichment analysis for Avena oat")
# Add command line arguments
p <- add_argument(p, "--gene", help = "input: gene list file", type = "character")
p <- add_argument(p, "--db", help = "input: org db name", type = "character")
p <- add_argument(p, "--rlib", help = "input: R library path", type = "character", default = "/mnt/e/software/R_library")
p <- add_argument(p, "--output", help = "output prefix", type = "character")
p <- add_argument(p, "--kegg_ko", help = "kegg to ko file", type = "character")
p <- add_argument(p, "--width", help = "output width", type = "numeric", default = 15)
p <- add_argument(p, "--height", help = "output height", type = "numeric", default = 15)
p <- add_argument(p, "--dpi", help = "output dpi", type = "numeric", default = 300)
p <- add_argument(p, "--pvalue", help = "pvalue cutoff", type = "numeric", default = 1)
p <- add_argument(p, "--qvalue", help = "qvalue cutoff", type = "numeric", default = 1)

# Parse the command line arguments
argv <- parse_args(p)
input_gene_list <- argv$gene
org_db <- argv$db
R_library <- argv$rlib
output <- argv$output
kegg_ko <- argv$kegg_ko
width <- argv$width
height <- argv$height
dpi <- argv$dpi
pvalueCutoff <- argv$pvalue
qvalueCutoff <- argv$qvalue


# 计程序运行时间
start_time <- Sys.time()
ko_annotaion <- kegg_ko
suppressMessages(library(tidyverse))
gene <- read_tsv(input_gene_list, col_names = F)
# 获取output文件夹路径
wd <- dirname(output)
wd <- paste(wd, "/", output, sep = "")
# 如果文件夹不存在新建
if (!dir.exists(wd)) {
    dir.create(wd)
}
# 设置工作路径
setwd(wd)
suppressMessages(library(clusterProfiler))
# 加载燕麦对应的org.db文件
if (org_db == "org.Asanfensan.eg.db") {
    suppressMessages(library(org.Asanfensan.eg.db, lib.loc = R_library))
    print("org.Asanfensan.eg.db is loaded!")
    pathway2gene <- AnnotationDbi::select(org.Asanfensan.eg.db,
        keys = keys(org.Asanfensan.eg.db),
        columns = c("Pathway", "KO")
    ) %>%
        na.omit() %>%
        dplyr::select(Pathway, GID) %>%
        mutate(Pathway = paste0("ko", Pathway))
    print(paste("Total pathway number:", length(unique(pathway2gene$Pathway)), sep = ""))
} else {
    # 终止程序
    stop("Please input correct org.db name!")
}
############## kegg enrichment analysis################
geneid <- gene$X1
print(paste("Total gene number:", length(geneid), sep = ""))

pathway2name <- read_tsv(ko_annotaion, col_names = F) %>%
    dplyr::select(Pathway = X1, GENENAME = X2) %>%
    distinct()

# pathway2name <- AnnotationDbi::select(org.Asativa.eg.db,
#                                       keys = keys(org.Asativa.eg.db),
#                                       columns = c("Pathway","GENENAME")) %>%
#   na.omit() %>%
#   dplyr::select(Pathway, GENENAME)

ekpsig <- enricher(as.factor(geneid),
    TERM2GENE = pathway2gene,
    TERM2NAME = pathway2name,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    pAdjustMethod = "BH",
    minGSSize = 1
)
ekp_results <- as.data.frame(ekpsig)
# head(ekp_results)
epk <- ekp_results[1:20, ]
# head(epk)
write.table(ekp_results, file = "All_kegg.txt", sep = "\t", quote = F, row.names = F)
# barplot(ekpsig, showCategory=20,color="pvalue",
#        font.size=16)
# library(see)
# barplot(ekpsig, showCategory=20,color="pvalue",
#        font.size=16)
# library(Cairo)
# Cairo.capabilities()
# 分割字符串为分子和分母

epk$GeneRatio_num <- sapply(strsplit(epk$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))


p <- ggplot(epk, aes(y = Description, x = GeneRatio_num, fill = p.adjust)) +
    geom_bar(stat = "identity", width = 0.8) + # 柱状图宽度设置
    scale_fill_gradient(low = "red", high = "blue") +
    labs(
        title = "KEGG Pathways Enrichment", # 设置标题、x轴和Y轴名称
        x = "Gene Ratio",
        y = "Pathway"
    ) +
    theme(
        axis.title.x = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 16)
    ) +
    theme_bw()
ggsave("keggplot.pdf", width = width, height = height)
ggsave("keggplot.png", width = width, height = height, dpi = dpi)
# dotplot(ekpsig,showCategory=20,color="pvalue",
#        font.size=16)
p <- ggplot(epk, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_colour_gradient(low = "green", high = "red") +
    labs(color = expression(p.adjust), size = "Gene number", x = "Gene Ratio", y = "Pathway", title = "KEGG Pathway enrichment")
ggsave("keggdot.pdf", width = width, height = height)
ggsave("keggdot.png", width = width, height = height, dpi = dpi)
# 更改ekp_results的列名

colnames(epk) <- c("id", "Description", "GeneRatio", "BgRatio", "pvalue", "padjust", "qvalue", "geneID", "Count")
# head(epk)
colnames(epk) <- c("id", "Description", "GeneRatio", "BgRatio", "RichFactor", "FoldEnrichment", "zScore", "pvalue", "padjust", "qvalue", "geneID", "Count", "GeneRatio_num")
write.table(epk, file = "kegg.tsv", sep = "\t", quote = F, row.names = F)
stop_time <- Sys.time()
print(paste("Kegg enrichment analysis is done! Time used:", stop_time - start_time, sep = ""))
