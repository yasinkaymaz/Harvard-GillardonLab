
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")
load("~/data/Tasic2018.seurat.Robj")

source("code/functions.R")

detable <- DElimma(L5ALM.ITPT)
DoHeatmap(object = L5ALM.ITPT, use.scaled = T, genes.use = rownames(detable), slim.col.label = TRUE, remove.key = TRUE, group.label.loc="top", group.label.rot=T)

load("~/data/Tasic2018-reAnalysis_L5ITPT.RData")

pdf("output/test.pdf",width = 20,height = 90)
VlnPlot(object = Tasic2018, features.plot = c("Nrgn","Slc17a7","Rgs4","Arpp19","Ptn","Arpp21"), use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()

pdf("output/test2.pdf",width = 20,height = 90)
VlnPlot(object = Tasic2018, features.plot = c("S100a10", "Fam3c", "Deptor", "Cnih3", "Rgs4", "Vstm2l"), use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()

marker.it <- markers.L5IT.vs.Gaba_Glu %>%
  add_column(gene=row.names(markers.L5IT.vs.Gaba_Glu)) %>%
  select(-p_val)%>%
  filter(p_val_adj < 0.01) %>%
  filter(pct.1 > 0.5) %>% filter( ((pct.1 - pct.2)/max(pct.1,pct.2) > 0.5))





save(Tasic2018, file="~/data/Tasic2018.seurat.Robj")
save(Tasic2016, file="~/data/Tasic2016.seurat.Robj")



InteractiveTable(FilterTable(markers.L5IT.vs.Gaba_Glu))
VlnPlot(object = Tasic2018, features.plot = head(FilterTable(markers.L5IT.vs.Gaba_Glu)$gene), use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)




Genes=FilterTable(markers.L5IT.vs.Gaba_Glu)$gene
n_genes = length(Genes)
output = "output/Tasic2018-markers.L5IT.vs.Glu.pdf"
SeuratObj = Tasic2018
pdf(output,width = 20,height = n_genes*5)
VlnPlot(object = SeuratObj, features.plot = Genes, use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()

FilterTable(markers.L5IT.vs.Gaba_Glu)$gene

FilterTable(annotateGenes(markers.L5IT.vs.Gaba_Glu))

DoHeatmap(object = Tasic2016, genes.use = FilterTable(markers.L5IT.vs.Gaba_Glu)$gene, slim.col.label = TRUE, remove.key = TRUE,group.label.loc="top",group.label.rot=T)

