#下载TCGA数据
# if (!requireNamespace("BiocManager", quietly=TRUE))#   install.packages("BiocManager")# BiocManager::install("TCGAbiolinks")
 library(TCGAbiolinks)library(dplyr)library(DT)library(SummarizedExperiment)#下面填入要下载的癌症种类
request_cancer=c("PRAD","BLCA","KICH","KIRC","KIRP")for (i in request_cancer) {
  cancer_type=paste("TCGA",i,sep="-")
  print(cancer_type)
  #下载临床数据
  clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")
  write.csv(clinical,file = paste(cancer_type,"clinical.csv",sep = "-"))
  
  #下载rna-seq的counts数据
  query <- GDCquery(project = cancer_type, 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  
  GDCdownload(query, method = "api", files.per.chunk = 100)
  expdat <- GDCprepare(query = query)
  count_matrix=assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"Counts.csv",sep = "-"))
  
  #下载miRNA数据
  query <- GDCquery(project = cancer_type, 
                    data.category = "Transcriptome Profiling", 
                    data.type = "miRNA Expression Quantification", 
                    workflow.type = "BCGSC miRNA Profiling")
  
  GDCdownload(query, method = "api", files.per.chunk = 50)
  expdat <- GDCprepare(query = query)
  count_matrix=assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"miRNA.csv",sep = "-"))
  
  #下载Copy Number Variation数据
  query <- GDCquery(project = cancer_type, 
                    data.category = "Copy Number Variation", 
                    data.type = "Copy Number Segment")
  
  GDCdownload(query, method = "api", files.per.chunk = 50)
  expdat <- GDCprepare(query = query)
  count_matrix=assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"Copy-Number-Variation.csv",sep = "-"))
  
  #下载甲基化数据
  query.met <- GDCquery(project =cancer_type,
                        legacy = TRUE,
                        data.category = "DNA methylation")
  GDCdownload(query.met, method = "api", files.per.chunk = 300)
  expdat <- GDCprepare(query = query)
  count_matrix=assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"methylation.csv",sep = "-"))}

##################################################################################################
# clinical data praperation


##################################################################################################
# data QC


##################################################################################################
# DEG


##################################################################################################
# upsetR veen


##################################################################################################
# circlize mapping
library(tidyverse)
library(igraph)
library(ggraph)
library(ggplot2)

path <-
    paste0(
        'https://raw.githubusercontent.com/rfordatascience/tidytuesday/',
        'master/data/2019/2019-11-12/'
    )

data <- read_csv(paste0(path, 'loc_cran_packages.csv'))

# most popular programming languages from TIOBE Index (Nov. 2019) found in data
# (only languages with position <= 16 are considered)
popular_languages <- c(
    'Java', 'C', 'Python', 'C++', 'C#', 'Visual Basic', 'JavaScript', 'PHP', 'SQL', 'Ruby', 'Objective C++', 'Assembly', 'R'
)

# number of packages to display
number_of_pkgs <- 300

# find largest packages written in popular languages
top_packages <- data %>%
    filter(language %in% popular_languages) %>%
    group_by(pkg_name) %>%
    summarize(total_code = sum(code)) %>%
    arrange(desc(total_code)) %>%
    head(number_of_pkgs) %>%
    select(pkg_name, total_code)

# all popular languages per package
top_languages_per_pkg <- data %>%
    filter(
        pkg_name %in% top_packages$pkg_name,
        language %in% popular_languages
    ) %>%
    arrange(pkg_name, desc(code)) %>%
    group_by(pkg_name) %>%
    mutate(
        main = row_number() == 1, # main language of package should be opaque
        total_code = sum(code)
    ) %>%
    ungroup() %>%
    select(language, pkg_name, code, total_code, main)

# only following languages found in given packages
(top_languages <- top_languages_per_pkg %>%
        pull(language) %>%
        unique %>%
        sort)

top_language_colors <- c(
    '#efb306',
    '#eb990c',
    '#e8351e',
    '#cd023d',
    '#852f88',
    '#4e54ac',
    '#0f8096',
    '#7db954',
    '#17a769',
    '#000000'
)

names(top_language_colors) <- c(
    'Assembly',
    'C',
    'C++',
    'JavaScript',
    'Java',
    'R',
    'Python',
    'Ruby',
    'SQL',
    'All'
)

edges1 <- top_languages_per_pkg %>%
    transmute(from = language, to = pkg_name, total_code = code, main)

edges2 <- top_languages_per_pkg %>%
    count(language, wt = code, name = 'total_code') %>%
    transmute(
        from = '',
        to = language,
        total_code,
        main = TRUE
    )

edges <- bind_rows(edges1, edges2)

vertices1 <- top_languages_per_pkg %>%
    filter(main) %>%
    transmute(
        node = pkg_name, language, total_code, level = 1
    )

vertices2 <- edges2 %>%
    transmute(
        node = to, language = to, total_code, level = 2
    )

vertices3 <- tibble(
    node = '', language = NA, total_code = 0, level = 3
)

vertices <- bind_rows(vertices1, vertices2, vertices3) %>%
    mutate(
        radius = total_code**(1.8), # scaling circles
        language = factor(language, names(top_language_colors))
    ) %>%
    arrange(level, language, node)

graph <- graph_from_data_frame(edges, vertices = vertices)

# create custom layout by updating existing circle layout
layout <- create_layout(graph, layout = 'circle')

outer_circle <- layout %>%
    filter(level == 1) %>%
    mutate(language = factor(language, names(top_language_colors))) %>%
    arrange(language, desc(name)) %>%
    mutate(
        x = cos((row_number() - 1) / number_of_pkgs * 2 * pi),
        y = sin((row_number() - 1) / number_of_pkgs * 2 * pi)
    )

# positioning circle centers manually by specifying polar coords
angles <- c(3, 43, 119, 160, 178, 255, 350, 190, 340, 0)
radii <- c(0.8, 0.5, 0.6, 0.4, 0.65, 0.45, 0.6, 0.7, 0.38, 0)
centers <- tibble(
    x = radii * cos(angles / 180 * pi),
    y = radii * sin(angles / 180 * pi)
)
inner_circle <- bind_cols(centers, select(filter(layout, level != 1), -x, -y))

layout[] <- bind_rows(outer_circle, inner_circle) %>%
    arrange(ggraph.index)

ggraph(layout) +
    geom_edge_diagonal(
        aes(edge_color = node1.language, edge_alpha = as.factor(main)),
        edge_width = 0.3, show.legend = FALSE
    ) +
    geom_node_point(
        aes(size = radius, color = language),
        alpha = 0.6, show.legend = FALSE
    ) +
    geom_node_text(
        aes(
            x = 1.0175 * x,
            y = 1.0175 * y,
            label = name,
            angle = -((-node_angle(x, y) + 90) %% 180) + 90,
            filter = !(name %in% top_languages)
        ),
        size = 2, hjust = 'outward', family = 'Oswald'
    ) +
    geom_node_text(
        aes(
            x = x,
            y = y,
            label = name,
            filter = name %in% top_languages
        ),
        size = 6, hjust = 0.5, family = 'Oswald'
    ) +
    geom_node_text(
        aes(
            x = x,
            y = y - 0.045,
            label = ifelse(
                total_code > 1000,
                format(total_code, big.mark = ','),
                total_code
            ),
            filter = name %in% top_languages
        ),
        size = 3, hjust = 0.5, family = 'Oswald'
    ) +
    scale_edge_color_manual(values = top_language_colors) +
    scale_color_manual(values = top_language_colors) +
    scale_size_area(max_size = 150) +
    scale_edge_alpha_manual(values = c(0.15, 1)) +
    coord_fixed() +
    labs(
        title = 'LOC of Popular Programming Languages in 300 CRAN Packages',
        subtitle = 'considered are largest CRAN packages written in one (or more) of top 16 programming languages from TIOBE Index (Nov. 2019)',
        caption = '#tidytuesday 46|2019 spren9er'
    ) +
    theme_void() +
    theme(
        text = element_text(family = 'Oswald'),
        legend.position = c(0.645, 0.51),
        plot.title = element_text(
            face = 'bold', hjust = 0.5, size = 20, margin = margin(t = 45, b = 3)
        ),
        plot.subtitle = element_text(
            face = 'plain', hjust = 0.5, size = 13, margin = margin(t = 5, b = 3)),
        plot.caption = element_text(
            face = 'plain', color = '#dedede', size = 8, hjust = 1,
            margin = margin(b = 20)
        )
    )

ggsave(
    'images/tidytuesday_201946_cran_packages.png',
    width = 12, height = 12.5, dpi = 300
)

##################################################################################################
# complexheatmap
setwd("C://Users//Administrator//Desktop//paper//Figures//Figure5//UCEC/")
options(stringsAsFactors = F)

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(viridis)
library(seriation)
library(circlize)
library(ggpubr)
library(Cairo)
library(dplyr)

UCEC_entropy<-read.csv(file = "UCEC_signals_gene_entropy_matrix.csv",header = T,check.names = F)
UCEC_expression<-read.csv(file="",header = T,check.names = F)
UCEC_annotations<-read.csv(file = "UCEC_annotation.csv",header = T,check.names = F)
UCEC_important_genelist<-read.csv(file = "UCEC_SNE_genelist.csv",header = T,check.names = F)
UCEC_important_gene_expression<-merge(UCEC_important_genelist,UCEC_expression,by = 'Gene')

UCEC_entropy2<-UCEC_entropy %>% arrange(desc(Gene))
rownames(UCEC_entropy2)<-UCEC_entropy2$Gene
UCEC_entropy2<-UCEC_entropy2[,-1]

UCEC_important_gene_expression2<-UCEC_important_gene_expression %>% arrange(desc(Gene))
rownames(UCEC_important_gene_expression2)<-UCEC_important_gene_expression2$Gene
UCEC_important_gene_expression2<-UCEC_important_gene_expression2[,-1]

col_fun = colorRamp2(c(-4, -2, 0, 2,4), c("RoyalBlue4", "Orchid2", "white", "Tomato","red"))

#ESCA用下列注释
#ha_stage_left<-HeatmapAnnotation(AJCC = HNSC_annotations$AJCC,
#                                 col = list(AJCC = c("Stage I"="ForestGreen","Stage IIA"="LimeGreen","Stage IIB"="Yellow4","Stage IIIA"="IndianRed1","Stage IIIB"="Firebrick","Stage IV"="LightSalmon4")),
#                                 annotation_name_side = "left")
#ha_stage_right<-HeatmapAnnotation(AJCC = HNSC_annotations$AJCC,
#                                  col = list(AJCC = c("Stage I"="ForestGreen","Stage IIA"="LimeGreen","Stage IIB"="Yellow4","Stage IIIA"="IndianRed1","Stage IIIB"="Firebrick","Stage IV"="LightSalmon4")),
#                                 annotation_name_side = "right")

#HNSC用下列注释：
#ha_stage_left<-HeatmapAnnotation(AJCC = HNSC_annotations$AJCC,
#                                 col = list(AJCC = c("Stage I"="ForestGreen","Stage II"="LimeGreen","Stage III"="IndianRed1","Stage IV"="LightSalmon4")),
#                                 annotation_name_side = "left")
#ha_stage_right<-HeatmapAnnotation(AJCC = HNSC_annotations$AJCC,
#                                  col = list(AJCC = c("Stage I"="ForestGreen","Stage II"="LimeGreen","Stage III"="IndianRed1","Stage IV"="LightSalmon4")),
#                                  annotation_name_side = "right")
#READ用下列注释:
#ha_stage_left<-HeatmapAnnotation(AJCC = HNSC_annotations$AJCC,
#                                 col = list(AJCC = c("Stage I"="ForestGreen","Stage II"="LimeGreen","Stage III"="IndianRed1","Stage IV"="LightSalmon4")),
#                                 annotation_name_side = "left")
#ha_stage_right<-HeatmapAnnotation(AJCC = HNSC_annotations$AJCC,
#                                  col = list(AJCC = c("Stage I"="ForestGreen","Stage II"="LimeGreen","Stage III"="IndianRed1","Stage IV"="LightSalmon4")),
#                                  annotation_name_side = "right")

#UCEC用下列注释：
ha_stage_left<-HeatmapAnnotation(AJCC = UCEC_annotations$AJCC,
                                 col = list(AJCC = c("Stage IA"="#9DC8C8","Stage IB"="#58C9B9","Stage IC"="#519D9E","Stage IIA"="#D1B6E1","Stage IIB"="#D499B9","Stage IIIA"="#9055A2","Stage IIIB"="#8F2D56","Stage IV"="#49010F")),
                                 annotation_name_side = "left")
ha_stage_right<-HeatmapAnnotation(AJCC = UCEC_annotations$AJCC,
                                 col = list(AJCC = c("Stage IA"="#9DC8C8","Stage IB"="#58C9B9","Stage IC"="#519D9E","Stage IIA"="#D1B6E1","Stage IIB"="#D499B9","Stage IIIA"="#9055A2","Stage IIIB"="#8F2D56","Stage IV"="#49010F")),
                                 annotation_name_side = "right")


mat1<-as.matrix(UCEC_entropy2)
mat2<-as.matrix(UCEC_important_gene_expression2)

CairoPDF(file = "UCEC_entropy_vs_expression_all.pdf",width = 16,height = 12)
p1<-Heatmap(mat1,
        col = col_fun,
        show_column_names=F,
        show_row_names = T,
        cluster_rows = F,
        row_names_gp = gpar(fontsize = 2.5),
        cluster_columns = F,
        row_title = "Vital Genes",
        row_names_side = 'left',
        column_title = "UCEC Important Gene Entropy",
        column_title_side = 'top',
        heatmap_legend_param = list(title = 'Entropy'),
        bottom_annotation = ha_stage_left)

p2<-Heatmap(mat2,
            col = col_fun,
            show_column_names=F,
            show_row_names = T,
            row_names_gp = gpar(fontsize = 2.5),
            cluster_columns = F,
            cluster_rows = F,
            row_title = "Vital Genes",
            row_names_side = 'left',
            column_title = "UCEC Important Gene Expression",
            heatmap_legend_param = list(title = 'Expression'),
            column_title_side = 'top',
            bottom_annotation = ha_stage_right)

ht_list<-p1+p2

draw(ht_list,merge_legends=T,legend_border=T,gap=unit(0.4,"cm"))

dev.off()

##################################################################################################
# soft-clustering
setwd("C://Users//Administrator//Desktop//paper//Figures//Figure6//READ/")
options(stringsAsFactors = F)

library(Cairo)
library(dplyr)
library(Mfuzz)

UCEC_expression<-read.csv(file="genesall.csv",header = T,check.names = F,row.names = 1)

UCEC_expression<-data.matrix(UCEC_expression)

entropygenes<-read.table("sne.txt",header = T)
entropygenes<-entropygenes$Gene

expgenes<-read.table("exp.txt",header = T)
expgenes<-expgenes$Gene

eset <- new("ExpressionSet",exprs = UCEC_expression)
eset <- standardise(eset)
c <- 6
m <- mestimate(eset)
cl <- mfuzz(eset, c = c, m = m)
cl$size
CairoPDF("READmfuu.pdf",height = 8,width = 8)
mfuzz.plot2(
    eset1,
    cl1,
    mfrow=c(2,3),
    x11 = F,
    centre = T,
    colo = 'fancy',
    time.labels = c("Before","Critical","After"),)
dev.off()

genes<-as.data.frame(cl$cluster)
names(genes)<-"cluster"
genes1<-subset(genes,cluster==1)
genes1<-row.names(genes1)
genes3<-subset(genes,cluster==3)
genes3<-row.names(genes3)

cluster3sne<-as.data.frame(intersect(genes3,entropygenes))
cluster3exp<-as.data.frame(intersect(genes3,expgenes))
write.csv(cluster3sne,"cluster3sne.csv")
write.csv(cluster3exp,"cluster3exp.csv")


##################################################################################################
# WGCNA co-expression


##################################################################################################
# GSVA
setwd("C://Users//Administrator//Desktop//paper//GSVA//ESCA//")

options(stringsAsFactors = F)

library(GSVA)
library(GSEABase)
library(dplyr)
library(limma)

C2CP<-getGmt("C2CP.symbols.gmt")

#这里要改：
df<-read.csv("ESCA_tumor_exprssion_matrix.csv",header = T,check.names = F,row.names = 1)

df<-log(df)
df<-as.matrix(df)
cancerPath<-gsva(df, gset.idx.list=C2CP, kcdf="Gaussian", parallel.sz=2)
dim(cancerPath)
pm<-as.data.frame(cancerPath)

write.csv(pm,file = "ESCA_PATHWAT_MX.csv")

##################################################################################################
# KEGG annotation plot
setwd("")

library(ggplot2)
library(ggpubr)
library(Cairo)
darkpath<-read.csv("dark_gene_pathway.csv",header = T)
netpath<-read.csv("network_pathway.csv",header = T)

darkpath$p<- round(-log10(darkpath$P.value))
darkpath$count<-c(7,7,2,4,3,4,2,2,3,2,2,2)
popo1<-data.frame(darkpath$Term,darkpath$count,darkpath$p)
names(popo1)<-c("Term","Count","P")

netpath$p<- round(-log10(netpath$P.value))
netpath$count<-c(6,6,7,6,9,4,3,6,5,5,7,3)
popo2<-data.frame(netpath$Term,netpath$count,netpath$p)
names(popo2)<-c("Term","Count","P")


CairoPDF("darkpath.pdf",width = 9,height = 5)
ggplot(data = popo1)+
    geom_point(aes(x=P,y=Term,size=Count,color=P))+
    scale_color_gradient(low="#379392",high ="#E71D36")+
    labs(x="Log(p-value)",y="KEGG",title="Dark genes",color = expression(log10(p-value)))+
    scale_size_continuous(range=c(3,8))+
    theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
dev.off()


CairoPDF("networkpath.pdf",width = 9,height = 5)
ggplot(data = popo2)+
    geom_point(aes(x=Count,y=Term,size=P,color=Count))+
    scale_color_gradient(low="#379392",high ="#E71D36")+
    labs(x="Count",y="KEGG",title="Network genes",size=expression(-log[10](p-value)))+
    scale_size_continuous(range=c(3,8))+
    scale_x_continuous(breaks = c(2.5,3,4,5,6,7,8,9))+
    theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

##################################################################################################
# KEGG pathview
setwd("C://Users//Administrator//Desktop//paper//Figures//Figure6//READ")
options(stringsAsFactors = F)
library(pathview)

pathview(gene.data = readexp[,1],
            pathway.id = "04621", 
            gene.idtype ="entrez",
            species = "hsa", 
            out.suffix = "READ-p53", 
            kegg.native = T,
            key.pos = "topright") 

##################################################################################################
# Survival analysis
setwd("C://Users//Administrator//Desktop//paper//Figures/Figure6/ESCA/")

options(stringsAsFactors = F)

library(survival)
library(parallel)
library(survminer)
library(patchwork)
library(Cairo)
#cancer:便于复制
#ESCA
#HNSC
#READ
#UCEC
#Dead:0
#Alive:1

ESCA_entropy<-read.csv("ESCA_signals_gene_entropy_matrix.csv",header = T,check.names = F)
ESCA_expression<-read.csv("ESCA_exp_Clinical_information.csv",header = T,check.names = F)
splots <- list()

ESCA_surv1<-Surv(ESCA_entropy$time,ESCA_entropy$vital_status)
log_rank_p1<- apply(ESCA_entropy[7:length(names(ESCA_entropy))],2,function(values){
    group<-ifelse(values>=mean(na.omit(values)),"High","Low")
    kmfit2<- survival::survfit(ESCA_surv1~group,data = ESCA_entropy)
    data.survdiff<-survival::survdiff(ESCA_surv1~group)
    p.value = 1 - pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
})
log_rank_p1<-log_rank_p1[log_rank_p1<0.05]
gene_marker1<-as.data.frame(sort(log_rank_p))
gene_marker_name1<-row.names(gene_marker1)
gene_marker_name1
write.table(gene_marker1,file = 'ESCA_entropy_gene_marker.txt',sep = '\t')
# "PER2"   "DAPK3"  "HCST"   "RAET1E" "ULBP1"  "ULBP3" 
#请注意以下语句中，变量名和绘图中的title参数需要手动一次一次替换成对应的基因名，目前没找到能够直接一次性出图的办法
group1<-ifelse(ESCA_entropy$PER2>mean(ESCA_entropy$PER2),
              "High","Low")
sfit1<-survfit(Surv(time,vital_status)~group1,data = ESCA_entropy)
splots[[1]]<-ggsurvplot(sfit1,
           legend.title = NULL,
           conf.int = F,
           pval = T,
           title = 'Gene Entropy Level',
           risk.table = F,
           palette = c("#E7B800","#2E9FDF"),
           risk.table.col = 'strata',
           tables.height = 0.2,
           xlab = "Time In Days",
           pval.method = TRUE,
           surv.median.line = "hv",
           ncensor.plot = T)

ESCA_surv2<-Surv(ESCA_expression$time,ESCA_expression$vital_status)
log_rank_p2<- apply(ESCA_expression[7:length(names(ESCA_expression))],2,function(values){
    group<-ifelse(values>=median(na.omit(values)),"High","Low")
    kmfit2<- survival::survfit(ESCA_surv2~group,data = ESCA_expression)
    data.survdiff<-survival::survdiff(ESCA_surv2~group)
    p.value = 1 - pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
})
log_rank_p2<-log_rank_p2[log_rank_p2<0.05]
gene_marker2<-as.data.frame(sort(log_rank_p2))
gene_marker_name2<-row.names(gene_marker2)
gene_marker_name2
write.table(gene_marker2,file = 'ESCA_expression_gene_marker.txt',sep = '\t')
#注意，差异基因显著的不画图，只保存名称。
#相反，将上面熵基因对预后显著的在表达的层面画图以凸显表达不能用来预测预后。
group2<-ifelse(ESCA_expression$PER2>median(ESCA_expression$PER2),"High","Low")
sfit2<-survfit(Surv(time,vital_status)~group2,data = ESCA_expression)
splots[[2]]<-ggsurvplot(sfit2,
           legend.title = NULL,
           conf.int = F,
           pval = T,
           title = 'Gene Expression Level',
           risk.table = F,
           palette = c("#E7B800","#2E9FDF"),
           risk.table.col = 'strata',
           tables.height = 0.2,
           xlab = "Time In Days",
           pval.method = TRUE,
           surv.median.line = "hv",
           ncensor.plot = T)

CairoPDF(file = 'ESCA_SNEvsEXP_PER2.pdf',width = 14,height = 6,)
arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 1)
dev.off()

rm(list = ls())

##################################################################################################
# package list and version:
