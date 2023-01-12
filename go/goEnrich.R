#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
suppressMessages(library(argparser, quietly=TRUE))
# Create a parser
p <- arg_parser("run Go and KEGG enrichment analysis for Avena oat")
# Add command line arguments
p <- add_argument(p, "--gene", help="input: gene list file", type="character")
p <- add_argument(p, "--db", help="input: org db name", type="character")
p <- add_argument(p, "--rlib", help="input: R library path", type="character",default="/mnt/e/software/R_library")
p <- add_argument(p, "--output", help="output prefix", type="character")
p <- add_argument(p, "--width", help="output width", type="numeric",default=15)
p <- add_argument(p, "--height", help="output height", type="numeric",default=15)
p <- add_argument(p, "--dpi", help="output dpi", type="numeric",default=300)


# Parse the command line arguments
argv <- parse_args(p)
input_gene_list <- argv$gene
org_db <- argv$db
R_library <- argv$rlib
output <- argv$output
width <- argv$width
height <- argv$height
dpi <- argv$dpi
#########################################################################################################
# 计程序运行时间
start_time <- Sys.time()
# 获取output文件夹路径
wd <- dirname(output)
wd <- paste(wd, "/",output, sep = "")
# 新建文件夹
dir.create(wd, showWarnings = FALSE)
# 设置工作路径
setwd(wd)
display_number = c(15, 15, 15)
# install.packages('./org.Asativa.eg.db_1.0.tar.gz',repos = NULL,lib = '/mnt/e/software/R_library')
#remove.packages('org.Asativa.eg.db')
# 加载燕麦对应的org.db文件
if (org_db == 'sfs'){
    suppressMessages(library(org.Asativa.eg.db, lib.loc = R_library))
    print("org.Asativa.eg.db is loaded!")
    org_db <- "org.Asativa.eg.db"
}else{
    # 终止程序
    stop("Please input correct org.db name!")
}
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))

gene <- read_tsv(paste('../',input_gene_list,sep=""),col_names = F)
geneid <- gene$X1
ego_BP <- enrichGO(OrgDb=org_db,
                   keyType = "GID", 
                   gene = geneid,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1, 
                   ont = "BP",
                   readable=FALSE)

# 打印程序运行过程
bp_count <- nrow(as.data.frame(ego_BP))
print(paste("BP enrichment analysis is done!, BP count is:",bp_count,sep=""))
if (bp_count < display_number[1]) {
  display_number[1] <- bp_count
}

ego_result_BP <- na.omit(as.data.frame(ego_BP))[1:display_number[1], ] 
# ego_result_BP <- ego_result_BP[order(ego_result_BP$Count),]


ego_CC <- enrichGO(OrgDb=org_db,
                   keyType = "GID", 
                   gene = geneid,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1, 
                   ont = "CC",
                   readable=FALSE)

cc_count <- nrow(as.data.frame(ego_CC))
print(paste("CC enrichment analysis is done!, CC count is:",cc_count,sep=""))
if (cc_count < display_number[2]) {
  display_number[2] <- cc_count
}                   
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ] 
# ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),]
ego_MF <- enrichGO(gene = geneid,                       #差异基因 vector 
                   keyType = "GID",                                   #差异基因的 ID 类型，需要是 OrgDb 支持的 
                   OrgDb = org_db,                               #对应的OrgDb 
                   ont = "MF",                                             #GO 分类名称，CC BP MF 
                   pvalueCutoff = 1,
                   qvalueCutoff = 1, 
                   readable = FALSE)

mf_count <- nrow(as.data.frame(ego_MF))
print(paste("MF enrichment analysis is done!, MF count is:",mf_count,sep=""))
if (mf_count < display_number[3]) {
  display_number[3] <- mf_count
}                   
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ] 
# ego_result_MF <- ego_result_MF[order(ego_result_MF$Count),]

go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                           Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                           GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                           type=factor(c(rep("biological process", display_number[1]), rep("cellular component",display_number[2]),rep("molecular function", display_number[3])), levels=c("molecular function", "cellular component", "biological process")),
                           gene=c(ego_result_BP$geneID, ego_result_CC$geneID, ego_result_MF$geneID)) %>% na.omit()
# 生成总表写出
print('GO enrichment analysis is done! write out all GO terms!')
all_BP <- as.data.frame(ego_BP) %>% na.omit()
all_CC <- as.data.frame(ego_CC) %>% na.omit()
all_MF <- as.data.frame(ego_MF) %>% na.omit()
# 合并三个表
all_GO <- rbind(all_BP, all_CC, all_MF)
# 写出
write_tsv(all_GO, 'all_GO.tsv')
write.table(go_enrich_df,file = "goplot_data.txt",sep = "\t",quote = F,row.names = F)
## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
labels=(go_enrich_df$Description)
names(labels) = rev(1:nrow(go_enrich_df))
## colors for bar // green, blue, orange
CPCOLS <- c( "#66C3A5", "#FD8D62","#8DA1CB")
print('GO enrichment analysis is done! plot GO terms!')
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50",size = 10)) +
  labs(title = "The Most Enriched GO Terms")
ggsave("goplot.pdf",device = cairo_pdf,width =width, height =height)
ggsave("goplot.png",width = width, height = height, dpi = dpi)
# 结束时间
stop_time <- Sys.time()
print(paste('GO enrichment analysis is done! Time used:',stop_time-start_time,sep=""))