options(stringsAsFactors = F)
suppressMessages(library(argparser, quietly=TRUE))

p <- arg_parser("make orgdb for avena")
p <- add_argument(p, "--genelist", help="gene list")
p <- add_argument(p, "--gene2go", help="gene2go")
p <- add_argument(p, "--gene2kegg", help="gene2kegg")
p <- add_argument(p, "--tax_id", help="tax_id")
p <- add_argument(p, "--genus", help="genus")
p <- add_argument(p, "--species", help="species")

argv <- parse_args(p)
all_gene_id <- argv$genelist
gene2go <- argv$gene2go
gene2kegg <- argv$gene2kegg
tax_id <- argv$tax_id
genus <- argv$genus
species <- argv$species

library(tidyverse)
library(clusterProfiler)
library(pkgbuild)


avena_gene <- read_tsv(all_gene_id,col_names = F,show_col_type=FALSE) %>%
  dplyr::select(GID = X1, GENENAME = X2) %>% na.omit()

gene2go <- read_delim(gene2go, "\t",
                      escape_double = FALSE, trim_ws = TRUE,col_names = F,show_col_type=FALSE) %>%
  dplyr::select(GID = X1, GO = X2) %>%
  dplyr::mutate(EVIDENCE = "IEA") %>% na.omit()

gene2ko <- read_delim(gene2kegg,
                      "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F,show_col_type=FALSE) %>%
  dplyr::select(GID = X1, KO = X2) %>%
  na.omit()
gene2pathway <- read_delim(gene2kegg,
                           "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F,show_col_type=FALSE) %>%
  dplyr::select(GID=X1, Pathway=X3) %>%
  na.omit()

gene2pathway[gene2pathway==""] <- NA
library(AnnotationForge)
#STEP6： 制作自己的Orgdb
# 查询物种的Taxonomy，例如要查avena
# https://www.ncbi.nlm.nih.gov/taxonomy/?term=avena
# tax_id = "4498"
# genus = "avena"
# species = "sativa"
#gene2go <- unique(gene2go)
gene2go<-gene2go[!duplicated(gene2go),]
gene2ko<-gene2ko[!duplicated(gene2ko),]
gene2pathway<-gene2pathway[!duplicated(gene2pathway),]
gene2gos <- dplyr::distinct(gene2go)
gene2kos<- dplyr::distinct(gene2ko)
# write.table(gene2gos,file = "gene2go.txt",sep = "\t",quote = F,row.names = F)
makeOrgPackage(gene_info=avena_gene,
               go=gene2gos,
               ko=gene2kos,
               maintainer = "yukaiquan <1962568272@qq.com>",
               author = "yukaiquan",
               pathway=gene2pathway,
               version="2.0",
               outputDir = '.',
               tax_id=tax_id,
               genus=genus,
               species=species,
               goTable="go")

