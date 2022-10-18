## http://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html

#cnetplot--visualize networks
x <- list(A = letters[1:10], B=letters[5:12], C=letters[sample(1:26, 15)])
x
p10 <- cnetplot(x)
p10


d <- setNames(rnorm(26), letters)
d
p11 <- cnetplot(x, foldChange=d) + 
  scale_color_gradient2(name='associated data', low='darkgreen', high='firebrick')
p11
cowplot::plot_grid(p10, p11, ncol=2, labels=LETTERS[1:2]) 

### bar plots
edo <- enrichDGN(de)
barplot(edo, showCategory=20) 

### dot plots
# gene set enrichment analysis
edo2 <- gseDO(geneList)
dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

### gene concept networks
# barplot() and dotplot() only displayed most significant or selected enriched terms
# to display which genes are involved in these significant terms. 
# genes may belong to multiple annotation categories
# cnetplot() depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network

## convert gene ID to Symbol
edox <- setReadable(edo, OrgDb='org.Hs.eg.db', keyType='ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
p1
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p2
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
p3
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

## node_label parameter supports 4 possible selections (“category”, “gene”, “all” and “none”),
p4 <- cnetplot(edox, node_label="category", 
               cex_label_category = 1.2) 
p4
p5 <- cnetplot(edox, node_label="gene", 
               cex_label_gene = 0.8) 
p5
p6 <- cnetplot(edox, node_label="all") 
p6
p7 <- cnetplot(edox, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
p7
cowplot::plot_grid(p4, p5, p6, p7, ncol=2, labels=LETTERS[1:4])

## also does heatmaps enriched categories
p8 <- heatplot(edox, showCategory=5)
p8
p9 <- heatplot(edox, foldChange=geneList, showCategory=5)
p9
cowplot::plot_grid(p8, p9, ncol=1, labels=LETTERS[1:2])

### treeplots
# treeplots are based on hierarchical clustering of enriched terms
# pairwise similarities of the enriched terms calculated by the pairwise_termsim() function
# by default uses Jaccard’s similarity index (JC) (set similarities)
# also use semantic similarity values if it is supported (e.g., GO, DO and MeSH)
edox2 <- pairwise_termsim(edox)
p12 <- treeplot(edox2)
p12
p13 <- treeplot(edox2, hclust_method = "average")
p13
aplot::plot_list(p12, p13, tag_levels='A')

### enrichment maps--networks of enriched terms
# Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. 
# In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify 
# functional module.
edo <- pairwise_termsim(edo)
p14 <- emapplot(edo)
p14
p15 <- emapplot(edo, cex_category=1.5)
p15
p16 <- emapplot(edo, layout="kk")
p16
p17 <- emapplot(edo, cex_category=1.5,layout="kk") 
p17
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

# Biological theme comparison
# enrichment maps can also be based on compareCluster output (from clusterProfiler)
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
xx <- pairwise_termsim(xx)                     
p18 <- emapplot(xx)
p18
p19 <- emapplot(xx, legend_n=2) 
p19
p20 <- emapplot(xx, pie="count")
p20
p21 <- emapplot(xx, pie="count", cex_category=1.5, layout="kk")
p21
cowplot::plot_grid(p18, p19, p20, p21, ncol=2, labels=LETTERS[1:4])

#upset plots
# alternative to cnetplot for visualizing the complex association between genes and gene sets
# emphasizes the gene overlapping among different gene sets.
upsetplot(edo)
# upsetplot(kk2) 
