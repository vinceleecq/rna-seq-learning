# 数据来自TCGA-LUAD（肺腺癌），随机抽取了20个样本，
# 其中10个为肿瘤样本Tumor，10个为正常样本Normal

# 学习内容：
# 数据读取和预处理
# PCA降维分析
# 使用DESeq2进行差异分析
# 对差异分析结果可视化,绘制火山图、热图
# 使用clusterProfiler进行GO、KEGG富集分析

## 安装R包----------------------------------------------------------------------
installed.packages("dplyr")
installed.packages("stringr")
installed.packages("ggplot2")
installed.packages("pheatmap")

BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")


## 加载R包----------------------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)
library(pheatmap)

library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)




# 数据读取和预处理--------------------------------------------------------------
exp_data = read.csv("./data/exp_data.csv",row.names = 1)
clin = read.csv("./data/clin.csv",row.names = 1)

dim(exp_data)
exp_data[1:3,1:4]
clin[1:3,1:4]
glimpse(clin)#使用dplyr包的函数查看临床数据的结构概览（列名、数据类型、前几个值）


# Note: read.csv 会将列名中的-自动转化为. ,-特殊符号不推荐出现在列名中
# 将clin数据框里面的样本名中的-也转化为.
rownames(clin) = str_replace_all(rownames(clin),pattern = "-",replacement = "\\.") #统一名称


# 取出样本分组信息
group_data = clin %>% dplyr::select(group)
head(group_data)


# 过滤低表达基因,保留那些至少在25%的样本中表达量大于10的基因
keep <- rowSums(exp_data > 10) >= round(20*0.25)
table(keep)
exp_data = exp_data[keep,]

exp_data_matrix = as.matrix(exp_data)



# PCA降维分析-------------------------------------------------------------------
# 在正式分析之前，我们可以先检查不同分组样本差异情况
# 基因表达矩阵往往有着上万个基因，光凭肉眼不可能看出样本之间的差异，
# 所以我们需要进行降维分析，将上万的基因（维度）降维成二维的形式 
# 通过降维分析，我们可以直观的看出样本间的差异。

# 使用count数据进行PCA

# count数据为均值-方差依赖，随着均值的增大，方差也在增大
plot(rowMeans(exp_data), rowSds(exp_data_matrix), 
     main='Raw counts: sd vs mean', 
     xlim=c(0,10000),
     ylim=c(0,5000))


# 直接用原始 count 做 PCA、聚类或可视化会会更容易受到高表达基因影响，
# 可以使用方差稳定转换或log进行处理。
# VST : variance stabilizing transformation
# vst首先对count数据做测序深度校正，然后在进行方差稳定化变换 输出结果适合做降维或者聚类
# ref:https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
exp_data_vst <- vst(exp_data_matrix)

plot(rowMeans(exp_data_vst), rowSds(exp_data_vst), 
     main='VST counts: sd vs mean')


# 执行PCA
# 执行主成分分析，获取样本在主成分上的坐标
pca <- prcomp(t(exp_data_vst)) 
d <- data.frame(PC1 = pca$x[,1],PC2 = pca$x[,2],PC3 = pca$x[,3],PC4 = pca$x[,4])
d <- cbind(d,group_data)
head(d)

# 统计每个主成分解释的方差百分比
# 方差占比越大说明降维后的分布与原始数据越接近。
percentVar <- round(100*(pca$sdev^2 / sum( pca$sdev^2 )),1)  # sdev--标准差（standard deviation）  sdev^2--方差（variance）
names(percentVar) <- paste0("PC",1:length(percentVar))

mycolors = mycolors <- c("#0072b5","#bc2c29")
p1 <- ggplot(d, aes(PC1,PC2)) +
  geom_point(aes(colour = group),size = 2) +
  xlab(paste0("PC1, Variance: ", percentVar["PC1"], "%")) +
  ylab(paste0("PC2, Variance: ", percentVar["PC2"], "%")) +
  scale_color_manual(values = mycolors)+
  theme_bw()
p1
ggsave("./pca_plot.pdf",plot = p1,width = 5, height = 4)


# 当然我们也可以使用TPM数据，同样过滤低表达,然后进行log处理之后做PCA



# 使用DESeq2进行差异分析--------------------------------------------------------
# 发现哪些基因在疾病组中上调，哪些基因在疾病组中下调，这种变化有没有显著（P<0.05）
# ref:https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


# 分组因子化，指定Normal为参考组
group_data$group = as.factor(group_data$group) %>% relevel(ref = "Normal")

# 构建 DESeqDataSet 对象
dds = DESeqDataSetFromMatrix(countData = exp_data, # 基因表达矩阵
                             colData = group_data, # 样本信息
                             design= ~ group)      # 告诉DESeq2 要比较什么因素


# 差异分析
dds = DESeq(dds)

# 查看保存的结果名称
resultsNames(dds)

# 提取差异分析结果
res = results( dds, name= "group_Tumor_vs_Normal")
head(res)

deg_res = as.data.frame(res)
deg_res = deg_res %>% 
    arrange(pvalue)
head(deg_res)


# 筛选差异基因
abs_logFC_cut = 2
p_cut = 0.05
deg_res = deg_res %>% 
  mutate(Change = case_when(
    log2FoldChange >= 2 & padj < 0.05 ~ "Up",
    log2FoldChange <= -2 & padj < 0.05 ~ "Down",
    T ~ "No_Sig"
  ))

deg_res$Change %>% table()
write.csv(deg_res,file = "./差异分析结果.csv")



# 差异分析结果可视化----------------------------------------------------------------

# 火山图
p <- ggplot(data = deg_res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Change),size = 1,alpha = 0.8) +
  geom_vline(xintercept = 2,color = "grey", linetype="longdash",size = 0.4) +
  geom_vline(xintercept = -2,color = "grey", linetype="longdash",size = 0.4) +
  geom_hline(yintercept = -log10(0.05) ,color = "grey",linetype="longdash", size = 0.4) +
  scale_color_manual(values = c("#0072b5", "grey", "#bc2c29"))+
  theme_bw()
p
ggsave(filename = "./火山图.pdf",width = 6,height = 5)

# 差异基因热图
diff_genes = deg_res %>% filter(Change %in% c("Up","Down")) %>% rownames()
exp_data_deg = exp_data_vst[diff_genes,]
dim(exp_data_deg)


# 对每个基因的表达水平进行scale:(x−mean(x))/sd(x) 从而让各基因在同一尺度下比较
# 避免热图中极值主导颜色
scale_data <- t(scale(t(exp_data_deg))) %>% as.data.frame()

pheatmap(scale_data,
         treeheight_row = 0,
         scale = "none",
         fontsize = 10,
         border_color = NA,
         show_colnames = T,
         show_rownames = F,
         cluster_cols = T,
         cluster_rows = T 
)


# 添加列注释
anno_col = group_data %>% arrange(group)
anno_colors = list(group = c( Normal = "#0072b5",Tumor = "#bc2c29" ))
main_colors = colorRampPalette(c("#0F88D4","white","pink"))(20)

pdf("./差异基因热图.pdf",width = 6,height = 5)
pheatmap(scale_data,
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         color = main_colors,
         treeheight_row = 0,
         scale = "none",
         fontsize = 10,
         border_color = NA,
         show_colnames = F,
         show_rownames = F,
         cluster_cols = T,
         cluster_rows = T 
)
dev.off()


# GO富集分析--------------------------------------------------------------------
# 分析这些差异基因主要参与哪些生物过程，涉及哪些通路。从而可以推断这些通路可能与疾病的发生相关
# ref:https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html


# GO
# GO Gene Ontology（基因本体论）是一个数据库
# 从三个方面对基因进行描述
# 生物过程（Biological Process, BP）描述基因产物参与的生物学目标或过程，如：如细胞分裂、代谢过程或信号传导。
# 分子功能（Molecular Function, MF）描述基因产物的生化活动或作用，如：酶活性、结合功能或转运功能
# 细胞组分（Cellular Component, CC）描述基因产物所存在或作用的细胞位置，如：细胞器、细胞结构及特定的生物膜区域等。


go_enrich <- enrichGO(gene        = diff_genes,
                      OrgDb         = "org.Hs.eg.db",  # 物种对应的基因注释数据库  # https://bioconductor.org/packages/release/BiocViews.html#___Organism 
                      keyType       = "SYMBOL",        # SYMBOL or ENSEMBL or ENTREZID
                      ont           = "ALL",           # ALL BP MF CC
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      minGSSize = 10,
                      readable = T)

go_enrich_df <- go_enrich@result
glimpse(go_enrich_df)
write.csv(go_enrich_df,file = "./GO_富集分析结果.csv")

p_facet <- dotplot(go_enrich,split="ONTOLOGY",x = "GeneRatio",showCategory = 10)+
  facet_grid(ONTOLOGY~., scale='free',space = "free_y")+
  scale_y_discrete(labels= function(x) str_wrap(x,width = 55) )+
  scale_size(range = c(3,8))

p_facet
ggsave("./GO_dotplot_facet.pdf",plot = p_facet,width = 8,height = 10)


# 只绘制BP
filter_go_BP <- go_enrich[go_enrich@result$ONTOLOGY == "BP",asis=T]
p_bp <- dotplot(filter_go_BP,showCategory = 10,x = "GeneRatio") +
  scale_size(range = c(3,8)) +
  scale_y_discrete(labels= function(x) str_wrap(x,width = 55))
p_bp
ggsave("./GO_dotplot_BP.pdf",plot = p_bp,width = 8,height = 6)




# KEGG富集分析------------------------------------------------------------------

entrez <- bitr(geneID = diff_genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
head(entrez)
entrez_gene <- entrez$ENTREZID

kegg_enrich <- enrichKEGG(
  gene = entrez_gene,
  organism  = "hsa",         # KEGG数据库人对应用hsa 其他物种参考https://www.genome.jp/kegg/catalog/org_list.html
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05,
)
kegg_enrich

# 把结果里的ENTREZID转化为SYMBOL格式方便阅读
kegg_enrich <- setReadable(kegg_enrich,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
kegg_enrich_df <- kegg_enrich@result

write.csv(kegg_enrich_df,file = "./KEGG_富集分析结果.csv")


p_kegg <- dotplot(kegg_enrich,showCategory = 10,x = "GeneRatio") +
  scale_size(range = c(3,8)) +
  scale_y_discrete(labels= function(x) str_wrap(x,width = 45))

p_kegg
ggsave("./KEGG_dotplot.pdf",plot = p_kegg,width = 8,height = 6)


p_kegg_bar <- 
  barplot(kegg_enrich,showCategory = 10)+
  scale_y_discrete(labels= function(x) str_wrap(x,width = 45) )+
  scale_size(range = c(3,8))
p_kegg_bar
ggsave("./KEGG_barplot.pdf",plot = p_kegg_bar,width = 8,height = 6)




