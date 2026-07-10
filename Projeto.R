library(scRNAseq)
library(tidyverse)
library(umap)
library(pheatmap)
library(randomForest)
library(caret)
library(matrixStats)
library(gridExtra)
library(mclust)
library(dendextend)
library(Rtsne)
library(kernlab)
library(e1071)
library(xgboost)
library(class)
library(RColorBrewer)
library(vegan)
library(pROC)
library(cluster)
library(aricode)
library(igraph)

sce<-ZeiselBrainData()
counts<-assay(sce,"counts")
genes_keep<-rowSums(counts>0)>=10
counts<-counts[genes_keep,]
totals<-colSums(counts)
norm<-t(t(counts)/totals)*10000
logcounts<-log(norm+1)
gene_vars<-rowVars(logcounts)
variable_genes<-order(gene_vars,decreasing=TRUE)[1:2000]
expr<-logcounts[variable_genes,]
expr_scaled<-t(scale(t(expr)))
expr_t<-t(expr_scaled)
cell_labels<-sce$level1class

cat("\n========== ANÁLISE DESCRITIVA DO BANCO ==========\n")
cat("\n--- DIMENSÕES ---\n")
cat(paste("Células totais:",ncol(counts),"\n"))
cat(paste("Genes totais:",nrow(counts),"\n"))
cat(paste("Genes filtrados (expr > 0 em >=10 células):",sum(genes_keep),"\n"))
cat(paste("Genes selecionados:",length(variable_genes),"\n"))

# contagem
cat(paste("Mediana de genes/célula:",median(colSums(counts>0)),"\n"))
cat(paste("Média de genes/célula:",round(mean(colSums(counts>0)),1),"\n"))
cat(paste("Mediana de UMIs/célula:",median(colSums(counts)),"\n"))
cat(paste("Média de UMIs/célula:",round(mean(colSums(counts)),1),"\n"))
table_tipos<-table(cell_labels)
print(table_tipos)
cat(paste("\nTotal de tipos celulares:",length(unique(cell_labels)),"\n"))

# esparsidade

## prop. global (Matriz Inteira)

prop_zeros_global<-sum(counts==0)/(nrow(counts)*ncol(counts))*100
cat(paste("Esparsidade Global (Matriz Completa):",round(prop_zeros_global,2),"%\n"))

## esparsidade por célula
zeros_por_celula<-colSums(counts==0)/nrow(counts)*100
df_zeros<-data.frame(prop_zeros=zeros_por_celula,tipo=cell_labels)

## por tipo
stats_zeros<-df_zeros%>%group_by(tipo)%>%summarise(Media=mean(prop_zeros),Mediana=median(prop_zeros))
print(stats_zeros)

p_zeros<-ggplot(df_zeros,aes(x=reorder(tipo,prop_zeros,FUN=median),y=prop_zeros,fill=tipo))+
  geom_boxplot(alpha=0.8,outlier.size=0.5,outlier.alpha=0.3)+
  coord_flip()+
  theme_minimal(base_size=16)+
  theme(plot.background=element_rect(fill="white",color=NA),
        panel.background=element_rect(fill="white",color=NA),
        panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),
        panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold",size=18,hjust=0.5,color="#2C3E50"),
        axis.title=element_text(face="bold",size=14,color="#34495E"),
        axis.text=element_text(color="#5D6D7E",size=12),
        legend.position="none")+
  labs(title="Distribuição de Esparsidade (Zeros) por Tipo",x=NULL,y="% de Genes Não Detectados")+
  scale_fill_viridis_d(option="turbo")
print(p_zeros)

df_qc<-data.frame(total_counts=colSums(counts),n_genes=colSums(counts>0),tipo=cell_labels)
p_umi<-ggplot(df_qc,aes(x=total_counts))+geom_histogram(bins=50,fill="#FF5BAE",color="#D20062",alpha=0.8)+
theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),
panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",
linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),
axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12))+
ggtitle("Distribuição de UMIs por Célula")+labs(x="Total de UMIs (escala log10)",y="Número de Células")+
scale_x_log10(labels=scales::label_number(big.mark="."))
print(p_umi)
p_genes<-ggplot(df_qc,aes(x=n_genes))+geom_histogram(bins=50,fill="#A1DD70",color="#799351",alpha=0.8)+
theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),
panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),
panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),
axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12))+
ggtitle("Genes Detectados por Célula")+labs(x="Número de Genes Detectados",y="Número de Células")+
scale_x_continuous(labels=scales::label_number(big.mark="."))
print(p_genes)
p_scatter<-ggplot(df_qc,aes(x=total_counts,y=n_genes,color=tipo))+geom_point(size=2,alpha=0.6)+theme_minimal(base_size=16)+
theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),
panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),
plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),
axis.text=element_text(color="#5D6D7E",size=12),legend.position="bottom",legend.title=element_text(face="bold",size=12),
legend.text=element_text(size=11))+ggtitle("Relação entre UMIs e Genes Detectados")+labs(x="Total de UMIs (escala log10)",
y="Genes Detectados",color="Tipo Celular")+scale_x_log10(labels=scales::label_number(big.mark="."))+
scale_y_continuous(labels=scales::label_number(big.mark="."))+scale_color_viridis_d(option="turbo")
print(p_scatter)
df_comp<-data.frame(tipo=cell_labels)%>%group_by(tipo)%>%summarise(n=n())%>%mutate(prop=n/sum(n)*100)
p_comp<-ggplot(df_comp,aes(x=reorder(tipo,n),y=n,fill=tipo))+geom_col()+geom_text(aes(label=n),hjust=-0.2,size=3)+coord_flip()+theme_minimal()+ggtitle("Composição do Dataset por Tipo Celular")+labs(x="Tipo Celular",y="Número de Células")+scale_fill_viridis_d()+theme(legend.position="none")+ylim(0,max(df_comp$n)*1.1)
print(p_comp)

# estatísticas por tipo celular
stats_por_tipo<-df_qc%>%group_by(tipo)%>%summarise(n_celulas=n(),media_UMIs=round(mean(total_counts),0),mediana_UMIs=median(total_counts),media_genes=round(mean(n_genes),0))
print(stats_por_tipo)

# variabilidade gênica
cat(paste(round(length(variable_genes)/nrow(logcounts)*100,1),"% genes mais variáveis selecionados\n"))
cat(paste("Range de variância: [",round(min(gene_vars[variable_genes]),3),"-",round(max(gene_vars[variable_genes]),3),"]\n"))
df_var_genes<-data.frame(variancia=sort(gene_vars,decreasing=TRUE)[1:2000],rank=1:2000)
p_var_genes<-ggplot(df_var_genes,aes(x=rank,y=variancia))+geom_line(color="#2C3E50",linewidth=1.5)+geom_hline(yintercept=median(df_var_genes$variancia),linetype="dotted",color="#95A5A6",linewidth=0.8)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),plot.subtitle=element_text(hjust=0,color="#7F8C8D",size=13),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="none")+labs(title="Perfil de Variância dos Genes Selecionados",subtitle=paste0("Mediana da variância: ",round(median(df_var_genes$variancia),3)),x="Ordem Decrescente de Variância",y="Variância (escala log)")+scale_x_continuous(breaks=c(1,500,1000,1500,2000),expand=c(0.01,0))+scale_y_continuous(expand=c(0,0),limits=c(0,max(df_var_genes$variancia)*1.05))
print(p_var_genes)

# correlação de pearson
top30_genes<-names(sort(gene_vars,decreasing=TRUE)[1:30])
expr_top30<-expr[top30_genes,]
cor_matrix<-cor(t(expr_top30))

#fazer sem cluster primeiro, depois na expressão diferencial acho que faz mais sentido
pheatmap(cor_matrix,cluster_rows=FALSE,cluster_cols=FALSE, 
         scale="none",show_rownames=TRUE,show_colnames=TRUE,
         fontsize_row=7,fontsize_col=7,color=colorRampPalette(c("#2E86AB","#F6F6F6","#E63946"))(100),
         main="Correlação de Pearson entre 30 Genes Mais Variáveis")

#violin plot dos 9 genes mais vairáveis
top9_var_genes<-names(sort(gene_vars,decreasing=TRUE)[1:9])
expr_top9<-t(logcounts[top9_var_genes,])
df_violin<-data.frame(expr_top9,Tipo=cell_labels)
df_violin_long<-melt(df_violin,id.vars="Tipo",variable.name="Gene",value.name="Expressao")

p_violin_desc<-ggplot(df_violin_long,aes(x=Tipo,y=Expressao,fill=Tipo))+
  geom_violin(scale="width",trim=FALSE,alpha=0.8,linewidth=0.2)+
  facet_wrap(~Gene,scales="free_y",ncol=3)+
  theme_minimal(base_size=11)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text=element_text(face="bold.italic",size=10),
        legend.position="bottom",
        legend.title=element_text(face="bold"),
        plot.title=element_text(face="bold",size=14,hjust=0.5))+
  labs(title="Perfil de Expressão dos 9 Genes de Maior Variância",x=NULL,y="Log-Expressão")+
  scale_fill_viridis_d(option="turbo",name="Tipo Celular")
print(p_violin_desc)

# ============================================================================
# EXPRESSÃO DIFERENCIAL: GENES MARCADORES
# ============================================================================
encontrar_marcadores<-function(expr,labels,tipo1,tipo2,method="wilcox"){
  cells_tipo1<-labels==tipo1
  cells_tipo2<-labels==tipo2
  expr_tipo1<-expr[,cells_tipo1]
  expr_tipo2<-expr[,cells_tipo2]
  mean1<-rowMeans(expr_tipo1)
  mean2<-rowMeans(expr_tipo2)
  log2fc<-(mean1-mean2)/log(2)
  
  pct1<-rowMeans(expr_tipo1>0)*100
  pct2<-rowMeans(expr_tipo2>0)*100
  
  if(method=="wilcox"){
    pvals<-sapply(1:nrow(expr),function(i){wilcox.test(expr_tipo1[i,],expr_tipo2[i,])$p.value})
  }else{
    pvals<-sapply(1:nrow(expr),function(i){t.test(expr_tipo1[i,],expr_tipo2[i,])$p.value})
  }
  
  df_result<-data.frame(gene=rownames(expr),
                        log2fc=log2fc,
                        pvalue=pvals,
                        padj=p.adjust(pvals,method="fdr"),
                        mean_tipo1=mean1,
                        mean_tipo2=mean2,
                        pct_tipo1=pct1,
                        pct_tipo2=pct2,
                        stringsAsFactors=FALSE)
  df_result<-df_result%>%arrange(padj,desc(abs(log2fc)))
  return(df_result)
}

cat("\n========== ANÁLISE DE EXPRESSÃO DIFERENCIAL ==========\n")
tipos_unicos<-unique(cell_labels)

cat("\nGenes marcadores mais expressivos por tipo celular:\n")
for(tipo in tipos_unicos){
  outros<-cell_labels!=tipo
  expr_tipo<-expr[,cell_labels==tipo]
  expr_outros<-expr[,outros]
  mean_tipo<-rowMeans(expr_tipo)
  mean_outros<-rowMeans(expr_outros)
  log2fc<-(mean_tipo-mean_outros)/log(2) 
  variable_genes_tipo<-names(sort(log2fc,decreasing=TRUE)[1:10])
  cat(paste0("\n",tipo,":\n"))
  cat(paste(variable_genes_tipo,collapse=", "),"\n")
}

# ============================================================================
# COMPARAÇÃO: PYRAMIDAL CA1 VS OLIGODENDROCYTES
# ============================================================================
tipo1<-"pyramidal CA1"
tipo2<-"oligodendrocytes"

marcadores_ca1_vs_oligo<-encontrar_marcadores(expr,cell_labels,tipo1,tipo2)
sig_markers<-marcadores_ca1_vs_oligo%>%filter(padj<0.05,abs(log2fc)>1)

cat(paste("\nGenes significativos (FDR<0.05, |log2FC|>1):",nrow(sig_markers),"\n"))
cat(paste("  CA1:",sum(sig_markers$log2fc>0),"  |  Oligo:",sum(sig_markers$log2fc<0),"\n\n"))

print(head(sig_markers%>%select(gene,log2fc,padj,pct_tipo1,pct_tipo2),20))

top10_ca1<-sig_markers%>%filter(log2fc>0)%>%head(10)
top10_oligo<-sig_markers%>%filter(log2fc<0)%>%arrange(log2fc)%>%head(10)

# ============================================================================
# VOLCANO PLOT
# ============================================================================
marcadores_ca1_vs_oligo$categoria<-case_when(
  marcadores_ca1_vs_oligo$padj<0.05 & marcadores_ca1_vs_oligo$log2fc>1~"Pyramidal CA1",
  marcadores_ca1_vs_oligo$padj<0.05 & marcadores_ca1_vs_oligo$log2fc<(-1)~"Oligodendrocytes",
  marcadores_ca1_vs_oligo$padj<0.05~"Significativo (FC<2)",
  TRUE~"Não significativo"
)

n_ca1<-sum(marcadores_ca1_vs_oligo$categoria=="Pyramidal CA1")
n_oligo<-sum(marcadores_ca1_vs_oligo$categoria=="Oligodendrocytes")

p_volcano <- ggplot(marcadores_ca1_vs_oligo, aes(x=log2fc, y=-log10(padj))) +
  geom_point(aes(color=categoria), alpha=0.7, size=1.8) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="#34495E", linewidth=0.6) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="#34495E", linewidth=0.6) +
  scale_color_manual(values=c("Pyramidal CA1"="#FA812F",
                              "Oligodendrocytes"="#4A90E2",
                              "Significativo (FC<2)"="#95A5A6",
                              "Não significativo"="gray80")) +
  scale_y_continuous(expand=expansion(mult=c(0.01, 0.3))) +
  coord_cartesian(clip="off") +
  theme_minimal(base_size=16) +
  theme(plot.background=element_rect(fill="white", color=NA),
        panel.background=element_rect(fill="white", color=NA),
        panel.grid.major=element_line(color="#E8E8E8", linewidth=0.4),
        panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold", size=18, hjust=0.5, color="#2C3E50"),
        plot.subtitle=element_text(hjust=0.5, color="#7F8C8D", size=13),
        axis.title=element_text(face="bold", size=14, color="#34495E"),
        axis.text=element_text(color="#5D6D7E", size=12),
        legend.position="right",
        legend.title=element_text(face="bold", size=12),
        legend.text=element_text(size=11)) +

    labs(title="Volcano Plot: Pyramidal CA1 vs Oligodendrocytes",
       subtitle=paste0(n_ca1, " genes Pyramidal CA1  •  ", n_oligo, " genes Oligodendrocytes"),
       x="Log2 Fold Change",
       y="-Log10 (FDR)",
       color="Categoria")
print(p_volcano)

# ============================================================================
# HEATMAP: TOP MARCADORES
# ============================================================================
genes_plot<-c(top10_ca1$gene,top10_oligo$gene)
cells_plot<-cell_labels%in%c(tipo1,tipo2)
expr_plot<-expr[genes_plot,cells_plot]
expr_plot_scaled<-t(scale(t(expr_plot)))

anno_col<-data.frame(Tipo=cell_labels[cells_plot])
rownames(anno_col)<-colnames(expr_plot_scaled)
anno_colors<-list(Tipo=c("pyramidal CA1"="#FA812F","oligodendrocytes"="#4A90E2"))

pheatmap(expr_plot_scaled,
         scale="none",
         show_colnames=FALSE,
         show_rownames=TRUE,
         annotation_col=anno_col,
         annotation_colors=anno_colors,
         cluster_cols=TRUE,
         cluster_rows=TRUE,
         clustering_distance_rows="correlation",
         clustering_distance_cols="euclidean",
         color=colorRampPalette(c("#2E86AB","#F6F6F6","#E63946"))(100),
         fontsize_row=8,
         fontsize=10,
         main="Genes Marcadores: Pyramidal CA1 vs Oligodendrocytes",
         border_color=NA)

# ============================================================================
# HEATMAP: EXPRESSÃO MÉDIA POR TIPO
# ============================================================================
mean_expr<-do.call(cbind,lapply(tipos_unicos,function(t){rowMeans(expr[,cell_labels==t,drop=FALSE])}))
colnames(mean_expr)<-tipos_unicos
top_var_genes<-names(sort(apply(mean_expr,1,var),decreasing=TRUE)[1:50])
mean_expr_scaled<-t(scale(t(mean_expr[top_var_genes,])))
p_heat2<-pheatmap(mean_expr_scaled,
                  scale="none",
                  clustering_distance_rows="correlation",
                  clustering_distance_cols="euclidean",
                  show_rownames=TRUE,
                  fontsize_row=6,
                  fontsize_col=9,
                  fontsize=10,
                  color=colorRampPalette(c("#1A237E","#FFFFFF","#D32F2F"))(100),
                  main="Expressão Média: 50 Genes Mais Variáveis por Tipo Celular",
                  border_color=NA,
                  angle_col=45,
                  cutree_rows=5,
                  cutree_cols=3)

# ============================================================================
# REDUÇÃO DE DIMENSIONALIDADE: PCA
# ============================================================================
set.seed(335705)
pca_result<-prcomp(expr_t)
var_exp<-summary(pca_result)$importance[2,1:10]*100
var_exp_cumsum<-cumsum(var_exp)
cat(paste("\nVariância explicada pelos 10 primeiros PCs:",round(var_exp_cumsum[10],1),"%\n"))
cat(paste("PC1 explica:",round(var_exp[1],1),"% | PC2 explica:",round(var_exp[2],1),"%\n"))
df_pca<-data.frame(PC1=pca_result$x[,1],PC2=pca_result$x[,2],tipo=cell_labels)
p1<-ggplot(df_pca,aes(x=PC1,y=PC2,color=tipo))+geom_point(size=2.5,alpha=0.7)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),
panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),
plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),
legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Projeção PCA dos Tipos Celulares")+
labs(x=paste0("PC1 (",round(var_exp[1],1),"%)"),y=paste0("PC2 (",round(var_exp[2],1),"%)"),color="Tipo Celular")+scale_color_viridis_d(option="turbo")
print(p1)
df_var<-data.frame(PC=factor(paste0("PC",1:10),levels=paste0("PC",1:10)),variancia=var_exp,cumulativa=var_exp_cumsum)
p2<-ggplot(df_var,aes(x=PC,y=variancia))+geom_col(fill="#134686",color="#1A252F",alpha=0.8,width=0.7)+geom_line(aes(y=cumulativa,group=1),color="#FA812F",linewidth=1.2)+
geom_point(aes(y=cumulativa),color="#FA812F",size=3)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),
panel.background=element_rect(fill="white",color=NA),panel.grid.major.y=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.major.x=element_blank(),
panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),
axis.text=element_text(color="#5D6D7E",size=12),axis.text.x=element_text(angle=0))+ggtitle("Variância Explicada por Componente Principal")+
labs(x="Componente Principal",y="Variância Explicada (%)")+scale_y_continuous(sec.axis=sec_axis(~.,name="Variância Cumulativa (%)",breaks=seq(0,50,10)))+
annotate("text",x=8,y=var_exp_cumsum[10]+2,label=paste0("Total (PC1-10): ",round(var_exp_cumsum[10],1),"%"),color="#FA812F",fontface="bold",size=4)
print(p2)

# ============================================================================
# REDUÇÃO: t-SNE
# ============================================================================
set.seed(335705)
tsne_result<-Rtsne(expr_t,dims=2,perplexity=30,max_iter=5000,pca=F)
df_tsne<-data.frame(tSNE1=tsne_result$Y[,1],tSNE2=tsne_result$Y[,2],tipo=cell_labels)
p_tsne<-ggplot(df_tsne,aes(x=tSNE1,y=tSNE2,color=tipo))+geom_point(size=2.5,alpha=0.7)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Projeção t-SNE dos Tipos Celulares")+labs(x="UMAP 1",y="UMAP 2",color="Tipo Celular")+scale_color_viridis_d(option="turbo")
print(p_tsne)

# ============================================================================
# REDUÇÃO: UMAP
# ============================================================================
set.seed(335705)
umap_result<-umap(expr_t,n_neighbors=15,min_dist=0.1)
df_umap<-data.frame(UMAP1=umap_result$layout[,1],UMAP2=umap_result$layout[,2],tipo=cell_labels)
p3<-ggplot(df_umap,aes(x=UMAP1,y=UMAP2,color=tipo))+geom_point(size=2.5,alpha=0.7)+theme_minimal(base_size=16)+
theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),
panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),
plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),
axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),
legend.text=element_text(size=10))+ggtitle("Projeção UMAP dos Tipos Celulares")+labs(x="UMAP 1",y="UMAP 2",color="Tipo Celular")+
scale_color_viridis_d(option="turbo")
print(p3)

# ============================================================================
# REDUÇÃO: Kernel PCA (RBF)
# ============================================================================
set.seed(335705)
kpca_result<-kpca(~.,data=as.data.frame(expr_t),kernel="rbfdot",kpar=list(sigma=0.00001),features=2)
df_kpca<-data.frame(kPC1=pcv(kpca_result)[,1],kPC2=pcv(kpca_result)[,2],tipo=cell_labels)
p_kpca<-ggplot(df_kpca,aes(x=kPC1,y=kPC2,color=tipo))+geom_point(size=2.5,alpha=0.7)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Projeção Kernel PCA (RBF) dos Tipos Celulares")+labs(x="Kernel PC 1",y="Kernel PC 2",color="Tipo Celular")+scale_color_viridis_d(option="turbo")
print(p_kpca)

# ============================================================================
# REDUÇÃO: Isomap
# ============================================================================
set.seed(335705)
isomap_result<-isomap(dist(expr_t),ndim=2,k=10)
df_isomap<-data.frame(Iso1=isomap_result$points[,1],Iso2=isomap_result$points[,2],tipo=cell_labels)
p_isomap<-ggplot(df_isomap,aes(x=Iso1,y=Iso2,color=tipo))+geom_point(size=2.5,alpha=0.7)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Projeção Isomap dos Tipos Celulares")+labs(x="Isomap 1",y="Isomap 2",color="Tipo Celular")+scale_color_viridis_d(option="turbo")
print(p_isomap)

# ============================================================================
# COMPARAÇÃO: Avaliação Quantitativa dos Métodos
# ============================================================================
silhouette_pca<-silhouette(as.numeric(factor(cell_labels)),dist(df_pca[,1:2]))
sil_mean_pca<-mean(silhouette_pca[,3])
cat(paste("PCA:        ",round(sil_mean_pca,3),"\n"))
silhouette_tsne<-silhouette(as.numeric(factor(cell_labels)),dist(df_tsne[,1:2]))
sil_mean_tsne<-mean(silhouette_tsne[,3])
cat(paste("t-SNE:      ",round(sil_mean_tsne,3),"\n"))
silhouette_umap<-silhouette(as.numeric(factor(cell_labels)),dist(df_umap[,1:2]))
sil_mean_umap<-mean(silhouette_umap[,3])
cat(paste("UMAP:       ",round(sil_mean_umap,3),"\n"))
silhouette_kpca<-silhouette(as.numeric(factor(cell_labels)),dist(df_kpca[,1:2]))
sil_mean_kpca<-mean(silhouette_kpca[,3])
cat(paste("Kernel PCA: ",round(sil_mean_kpca,3),"\n"))
silhouette_isomap<-silhouette(as.numeric(factor(df_isomap$tipo)),dist(df_isomap[,1:2]))
sil_mean_isomap<-mean(silhouette_isomap[,3])
cat(paste("Isomap:     ",round(sil_mean_isomap,3),"\n"))
df_pca_comp<-df_pca%>%dplyr::select(PC1,PC2,tipo)%>%dplyr::rename(Dim1=PC1,Dim2=PC2)%>%dplyr::mutate(metodo="PCA")
df_tsne_comp<-df_tsne%>%dplyr::select(tSNE1,tSNE2,tipo)%>%dplyr::rename(Dim1=tSNE1,Dim2=tSNE2)%>%dplyr::mutate(metodo="t-SNE")
df_umap_comp<-df_umap%>%dplyr::select(UMAP1,UMAP2,tipo)%>%dplyr::rename(Dim1=UMAP1,Dim2=UMAP2)%>%dplyr::mutate(metodo="UMAP")
df_all<-rbind(df_pca_comp,df_tsne_comp,df_umap_comp)
p_all_red<-ggplot(df_all,aes(x=Dim1,y=Dim2,color=tipo))+geom_point(size=0.8,alpha=0.5)+facet_wrap(~metodo,scales="free")+theme_minimal(base_size=14)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),strip.text=element_text(face="bold",size=12,color="#2C3E50"),legend.position="bottom",legend.title=element_text(face="bold",size=11))+ggtitle("Comparação dos Principais Métodos de Redução")+scale_color_viridis_d(option="turbo")+labs(color="Tipo Celular")
print(p_all_red)

# ============================================================================
# REDUÇÃO (SOBRE 100 PCs)
# ============================================================================
pcs_input<-pca_result$x[,1:100]

set.seed(335705)
tsne_result<-Rtsne(pcs_input,dims=2,perplexity=30,max_iter=5000,check_duplicates=FALSE)
df_tsne<-data.frame(tSNE1=tsne_result$Y[,1],tSNE2=tsne_result$Y[,2],tipo=cell_labels)
p_tsne<-ggplot(df_tsne,aes(x=tSNE1,y=tSNE2,color=tipo))+geom_point(size=2.5,alpha=0.7)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Projeção t-SNE (sobre 100 PCs)")+labs(x="UMAP 1",y="UMAP 2",color="Tipo Celular")+scale_color_viridis_d(option="turbo")
print(p_tsne)

set.seed(335705)
umap_result<-umap::umap(pcs_input,n_neighbors=15,min_dist=0.1)
df_umap<-data.frame(UMAP1=umap_result$layout[,1],UMAP2=umap_result$layout[,2],tipo=cell_labels)
p3<-ggplot(df_umap,aes(x=UMAP1,y=UMAP2,color=tipo))+geom_point(size=2.5,alpha=0.7)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Projeção UMAP (sobre 100 PCs)")+labs(x="UMAP 1",y="UMAP 2",color="Tipo Celular")+scale_color_viridis_d(option="turbo")
print(p3)

set.seed(335705)
kpca_result<-kpca(~.,data=as.data.frame(pcs_input),kernel="rbfdot",kpar=list(sigma=0.1),features=2)
df_kpca<-data.frame(kPC1=pcv(kpca_result)[,1],kPC2=pcv(kpca_result)[,2],tipo=cell_labels)
p_kpca<-ggplot(df_kpca,aes(x=kPC1,y=kPC2,color=tipo))+geom_point(size=2.5,alpha=0.7)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Projeção Kernel PCA (sobre 100 PCs)")+labs(x="Kernel PC 1",y="Kernel PC 2",color="Tipo Celular")+scale_color_viridis_d(option="turbo")
print(p_kpca)

set.seed(335705)
isomap_result<-isomap(dist(pcs_input),ndim=2,k=10)
df_isomap<-data.frame(Iso1=isomap_result$points[,1],Iso2=isomap_result$points[,2],tipo=cell_labels)
p_isomap<-ggplot(df_isomap,aes(x=Iso1,y=Iso2,color=tipo))+geom_point(size=2.5,alpha=0.7)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Projeção Isomap (sobre 100 PCs)")+labs(x="Isomap 1",y="Isomap 2",color="Tipo Celular")+scale_color_viridis_d(option="turbo")
print(p_isomap)

# ============================================================================
# CLUSTERING: K-Means
# ============================================================================
set.seed(335705)
k_real<-length(unique(cell_labels))
cat(paste("Número real de tipos celulares:",k_real,"\n"))
kmeans_result<-kmeans(expr_t,centers=k_real,nstart=25)
df_umap_km<-data.frame(UMAP1=df_umap$UMAP1,UMAP2=df_umap$UMAP2,cluster_kmeans=as.factor(kmeans_result$cluster))
p_kmeans<-ggplot(df_umap_km,aes(x=UMAP1,y=UMAP2,color=cluster_kmeans))+geom_point(size=2.5,alpha=0.7)+stat_ellipse(level=0.95,linewidth=1)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("K-means Clustering projetado no UMAP")+labs(x="UMAP 1",y="UMAP 2",color="Cluster")+scale_color_viridis_d(option="turbo")
print(p_kmeans)
ari_kmeans<-adjustedRandIndex(cell_labels,kmeans_result$cluster)
nmi_kmeans<-NMI(cell_labels,kmeans_result$cluster)
cat(paste("\nARI (Adjusted Rand Index):",round(ari_kmeans,3),"\n"))
cat(paste("NMI (Normalized Mutual Information):",round(nmi_kmeans,3),"\n"))
table_kmeans<-table(Tipo_Real=cell_labels,Cluster_Kmeans=kmeans_result$cluster)
print(table_kmeans)

# ============================================================================
# CLUSTERING: Hierarchical
# ============================================================================
set.seed(335705)
dist_matrix<-dist(expr_t)
hc_result<-hclust(dist_matrix,method="ward.D2")
hc_clusters<-cutree(hc_result,k=k_real)
df_umap_hc<-data.frame(UMAP1=df_umap$UMAP1,UMAP2=df_umap$UMAP2,cluster_hc=as.factor(hc_clusters))
p_hc<-ggplot(df_umap_hc,aes(x=UMAP1,y=UMAP2,color=cluster_hc))+geom_point(size=2.5,alpha=0.7)+stat_ellipse(level=0.95,linewidth=1)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Hierarchical Clustering projetado no UMAP")+labs(x="UMAP 1",y="UMAP 2",color="Cluster")+scale_color_viridis_d(option="turbo")
print(p_hc)
ari_hc<-adjustedRandIndex(cell_labels,hc_clusters)
nmi_hc<-NMI(cell_labels,hc_clusters)
cat(paste("Hierarchical - ARI:",round(ari_hc,3),"| NMI:",round(nmi_hc,3),"\n"))
table_hc<-table(Tipo_Real=cell_labels,Cluster_Hierarquico=hc_clusters)
print(table_hc)

# ============================================================================
# CLUSTERING: Gaussian Mixture Model
# ============================================================================

# Excessivamente custoso de rodar, rodei apenas nos primeiros 100 PCS e deixo o código comentado aqui.

# set.seed(335705)
# gmm_result<-Mclust(expr_t,G=8)
# df_umap$cluster_gmm<-as.factor(gmm_result$classification)
# p_gmm<-ggplot(df_umap,aes(x=UMAP1,y=UMAP2,color=cluster_gmm))+geom_point(size=1.5,alpha=0.6)+stat_ellipse(level=0.95,size=1)+theme_minimal()+ggtitle("Gaussian Mixture Model (k=8)")+scale_color_viridis_d()
# print(p_gmm)
# table_gmm<-table(cell_labels,gmm_result$classification)
# print("Tabela: Tipo Real vs GMM")
# print(table_gmm)
# print(paste("BIC do GMM:",round(gmm_result$bic,2)))
# ari_gmm<-adjustedRandIndex(cell_labels,gmm_result$classification)
# nmi_gmm<-NMI(cell_labels,gmm_result$classification)
# cat(paste("GMM - ARI:",round(ari_gmm,3),"| NMI:",round(nmi_gmm,3),"| BIC:",round(gmm_result$bic,2),"\n"))
# 

set.seed(335705)
pcs_gmm<-pca_result$x[,1:100]
gmm_result<-Mclust(pcs_gmm,G=k_real)
df_tsne_gmm<-data.frame(tSNE1=tsne_result$Y[,1],tSNE2=tsne_result$Y[,2],cluster_gmm=as.factor(gmm_result$classification))
p_gmm<-ggplot(df_tsne_gmm,aes(x=tSNE1,y=tSNE2,color=cluster_gmm))+geom_point(size=2.5,alpha=0.7)+stat_ellipse(level=0.95,linewidth=1)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Gaussian Mixture Model projetado no t-SNE")+labs(x="UMAP 1",y="UMAP 2",color="Cluster")+scale_color_viridis_d(option="turbo")
print(p_gmm)
ari_gmm<-adjustedRandIndex(cell_labels,gmm_result$classification)
nmi_gmm<-NMI(cell_labels,gmm_result$classification)
cat(paste("GMM - ARI:",round(ari_gmm,3),"| NMI:",round(nmi_gmm,3),"| BIC:",round(gmm_result$bic,2),"\n"))
table_gmm<-table(Tipo_Real=cell_labels,Cluster_GMM=gmm_result$classification)
print(table_gmm)

# ============================================================================
# CLUSTERING: Spectral Clustering
# ============================================================================
set.seed(335705)
sc_result<-specc(as.matrix(expr_t),centers=k_real)
sc_clusters<-as.vector(sc_result)
df_umap_spec<-data.frame(UMAP1=df_umap$UMAP1,UMAP2=df_umap$UMAP2,cluster_spec=as.factor(sc_clusters))
p_spec<-ggplot(df_umap_spec,aes(x=UMAP1,y=UMAP2,color=cluster_spec))+geom_point(size=2.5,alpha=0.7)+stat_ellipse(level=0.95,linewidth=1)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Spectral Clustering projetado no UMAP")+labs(x="UMAP 1",y="UMAP 2",color="Cluster")+scale_color_viridis_d(option="turbo")
print(p_spec)
ari_spec<-adjustedRandIndex(cell_labels,sc_clusters)
nmi_spec<-NMI(cell_labels,sc_clusters)
cat(paste("Spectral - ARI:",round(ari_spec,3),"| NMI:",round(nmi_spec,3),"\n"))
table_spec<-table(Tipo_Real=cell_labels,Cluster_Spectral=sc_clusters)
print(table_spec)

# ============================================================================
# CLUSTERING: Louvain (Graph-based)
# ============================================================================
set.seed(335705)
# Construir grafo KNN
k<-15
dists<-as.matrix(dist(expr_t))
knn_idx<-t(apply(dists,1,function(x)order(x)[2:(k+1)]))
edges<-NULL
for(i in 1:nrow(knn_idx)){
  for(j in 1:k){
    edges<-rbind(edges,c(i,knn_idx[i,j]))
  }
}
graph<-graph_from_edgelist(edges,directed=FALSE)
graph<-simplify(graph)
louvain_result<-cluster_louvain(graph)
louvain_clusters<-as.numeric(membership(louvain_result))
df_umap_louvain<-data.frame(UMAP1=df_umap$UMAP1,UMAP2=df_umap$UMAP2,cluster_louvain=as.factor(louvain_clusters))
p_louvain<-ggplot(df_umap_louvain,aes(x=UMAP1,y=UMAP2,color=cluster_louvain))+geom_point(size=2.5,alpha=0.7)+stat_ellipse(level=0.95,linewidth=1)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Louvain Clustering projetado no UMAP")+labs(x="UMAP 1",y="UMAP 2",color="Cluster")+scale_color_viridis_d(option="turbo")
print(p_louvain)
ari_louvain<-adjustedRandIndex(cell_labels,louvain_clusters)
nmi_louvain<-NMI(cell_labels,louvain_clusters)
n_clusters_louvain<-length(unique(louvain_clusters))
cat(paste("Louvain - ARI:",round(ari_louvain,3),"| NMI:",round(nmi_louvain,3),"| Clusters:",n_clusters_louvain,"\n"))
table_louvain<-table(Tipo_Real=cell_labels,Cluster_Louvain=louvain_clusters)
print(table_louvain)

# ============================================================================
# CLUSTERING (SOBRE 100 PCs)
# ============================================================================
pcs_input<-pca_result$x[,1:100]
k_real<-length(unique(cell_labels))
cat(paste("Número real de tipos celulares:",k_real,"\n"))

# ============================================================================
# CLUSTERING: K-Means
# ============================================================================
set.seed(335705)
kmeans_result<-kmeans(pcs_input,centers=k_real,nstart=25)
df_umap_km<-data.frame(UMAP1=umap_result$layout[,1],UMAP2=umap_result$layout[,2],cluster_kmeans=as.factor(kmeans_result$cluster))
p_kmeans<-ggplot(df_umap_km,aes(x=UMAP1,y=UMAP2,color=cluster_kmeans))+geom_point(size=2.5,alpha=0.7)+stat_ellipse(level=0.95,linewidth=1)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("K-means (100 PCs) projetado no UMAP")+labs(x="UMAP 1",y="UMAP 2",color="Cluster")+scale_color_viridis_d(option="turbo")
print(p_kmeans)
ari_kmeans<-adjustedRandIndex(cell_labels,kmeans_result$cluster)
nmi_kmeans<-NMI(cell_labels,kmeans_result$cluster)
cat(paste("\nARI (Adjusted Rand Index):",round(ari_kmeans,3),"\n"))
cat(paste("NMI (Normalized Mutual Information):",round(nmi_kmeans,3),"\n"))
table_kmeans<-table(Tipo_Real=cell_labels,Cluster_Kmeans=kmeans_result$cluster)
print(table_kmeans)

# ============================================================================
# CLUSTERING: Hierarchical
# ============================================================================
set.seed(335705)
dist_matrix<-dist(pcs_input)
hc_result<-hclust(dist_matrix,method="ward.D2")
hc_clusters<-cutree(hc_result,k=k_real)
df_umap_hc<-data.frame(UMAP1=umap_result$layout[,1],UMAP2=umap_result$layout[,2],cluster_hc=as.factor(hc_clusters))
p_hc<-ggplot(df_umap_hc,aes(x=UMAP1,y=UMAP2,color=cluster_hc))+geom_point(size=2.5,alpha=0.7)+stat_ellipse(level=0.95,linewidth=1)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Hierarchical (100 PCs) projetado no UMAP")+labs(x="UMAP 1",y="UMAP 2",color="Cluster")+scale_color_viridis_d(option="turbo")
print(p_hc)
ari_hc<-adjustedRandIndex(cell_labels,hc_clusters)
nmi_hc<-NMI(cell_labels,hc_clusters)
cat(paste("Hierarchical - ARI:",round(ari_hc,3),"| NMI:",round(nmi_hc,3),"\n"))
table_hc<-table(Tipo_Real=cell_labels,Cluster_Hierarquico=hc_clusters)
print(table_hc)

# ============================================================================
# CLUSTERING: Gaussian Mixture Model
# ============================================================================
set.seed(335705)
gmm_result<-Mclust(pcs_input,G=k_real)
df_umap_gmm<-data.frame(UMAP1=umap_result$layout[,1],UMAP2=umap_result$layout[,2],cluster_gmm=as.factor(gmm_result$classification))
p_gmm<-ggplot(df_umap_gmm,aes(x=UMAP1,y=UMAP2,color=cluster_gmm))+geom_point(size=2.5,alpha=0.7)+stat_ellipse(level=0.95,linewidth=1)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("GMM (100 PCs) projetado no UMAP")+labs(x="UMAP 1",y="UMAP 2",color="Cluster")+scale_color_viridis_d(option="turbo")
print(p_gmm)
ari_gmm<-adjustedRandIndex(cell_labels,gmm_result$classification)
nmi_gmm<-NMI(cell_labels,gmm_result$classification)
cat(paste("GMM - ARI:",round(ari_gmm,3),"| NMI:",round(nmi_gmm,3),"| BIC:",round(gmm_result$bic,2),"\n"))
table_gmm<-table(Tipo_Real=cell_labels,Cluster_GMM=gmm_result$classification)
print(table_gmm)

# ============================================================================
# CLUSTERING: Spectral Clustering
# ============================================================================
set.seed(335705)
sc_result<-specc(pcs_input,centers=k_real)
sc_clusters<-as.vector(sc_result)
df_umap_spec<-data.frame(UMAP1=umap_result$layout[,1],UMAP2=umap_result$layout[,2],cluster_spec=as.factor(sc_clusters))
p_spec<-ggplot(df_umap_spec,aes(x=UMAP1,y=UMAP2,color=cluster_spec))+geom_point(size=2.5,alpha=0.7)+stat_ellipse(level=0.95,linewidth=1)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Spectral (100 PCs) projetado no UMAP")+labs(x="UMAP 1",y="UMAP 2",color="Cluster")+scale_color_viridis_d(option="turbo")
print(p_spec)
ari_spec<-adjustedRandIndex(cell_labels,sc_clusters)
nmi_spec<-NMI(cell_labels,sc_clusters)
cat(paste("Spectral - ARI:",round(ari_spec,3),"| NMI:",round(nmi_spec,3),"\n"))
table_spec<-table(Tipo_Real=cell_labels,Cluster_Spectral=sc_clusters)
print(table_spec)

# ============================================================================
# CLUSTERING: Louvain (Graph-based)
# ============================================================================
set.seed(335705)
k<-15
dists<-as.matrix(dist(pcs_input))
knn_idx<-t(apply(dists,1,function(x)order(x)[2:(k+1)]))
edges<-NULL
for(i in 1:nrow(knn_idx)){
  for(j in 1:k){
    edges<-rbind(edges,c(i,knn_idx[i,j]))
  }
}
graph<-graph_from_edgelist(edges,directed=FALSE)
graph<-simplify(graph)
louvain_result<-cluster_louvain(graph)
louvain_clusters<-as.numeric(membership(louvain_result))
df_umap_louvain<-data.frame(UMAP1=umap_result$layout[,1],UMAP2=umap_result$layout[,2],cluster_louvain=as.factor(louvain_clusters))
p_louvain<-ggplot(df_umap_louvain,aes(x=UMAP1,y=UMAP2,color=cluster_louvain))+geom_point(size=2.5,alpha=0.7)+stat_ellipse(level=0.95,linewidth=1)+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="right",legend.title=element_text(face="bold",size=12),legend.text=element_text(size=10))+ggtitle("Louvain (100 PCs) projetado no UMAP")+labs(x="UMAP 1",y="UMAP 2",color="Cluster")+scale_color_viridis_d(option="turbo")
print(p_louvain)
ari_louvain<-adjustedRandIndex(cell_labels,louvain_clusters)
nmi_louvain<-NMI(cell_labels,louvain_clusters)
n_clusters_louvain<-length(unique(louvain_clusters))
cat(paste("Louvain - ARI:",round(ari_louvain,3),"| NMI:",round(nmi_louvain,3),"| Clusters:",n_clusters_louvain,"\n"))

# ============================================================================
# COMPARAÇÃO: Métricas de Clustering
# ============================================================================
df_metricas<-data.frame(Metodo=c("K-means","Hierarchical","Spectral","Louvain"),ARI=c(ari_kmeans,ari_hc,ari_spec,ari_louvain),NMI=c(nmi_kmeans,nmi_hc,nmi_spec,nmi_louvain))
print(df_metricas)
p_ari<-ggplot(df_metricas,aes(x=reorder(Metodo,ARI),y=ARI,fill=Metodo))+geom_col(alpha=0.8)+coord_flip()+theme_minimal(base_size=16)+theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),panel.grid.major.x=element_line(color="#E8E8E8",linewidth=0.4),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=18,hjust=0,color="#2C3E50"),axis.title=element_text(face="bold",size=14,color="#34495E"),axis.text=element_text(color="#5D6D7E",size=12),legend.position="none")+ggtitle("Adjusted Rand Index por Método")+labs(x="Método de Clustering",y="ARI")+scale_fill_viridis_d(option="turbo")+geom_text(aes(label=round(ARI,3)),hjust=-0.2,size=5,fontface="bold")
print(p_ari)

# ============================================================================
# PREPARAÇÃO (para evitar data leakage)
# ============================================================================
set.seed(335705)
trainIndex<-createDataPartition(cell_labels,p=0.7,list=FALSE)
vars_tr<-rowVars(logcounts[,trainIndex])
genes_sel<-order(vars_tr,decreasing=TRUE)[1:2000]
tr_raw<-t(logcounts[genes_sel,trainIndex])
te_raw<-t(logcounts[genes_sel,-trainIndex])
tr_s<-scale(tr_raw)
te_s<-scale(te_raw,center=attr(tr_s,"scaled:center"),scale=attr(tr_s,"scaled:scale"))
treino<-data.frame(tr_s,tipo=as.factor(cell_labels[trainIndex]))
teste<-data.frame(te_s,tipo=as.factor(cell_labels[-trainIndex]))
pca_clean<-prcomp(tr_s)
pcs_treino<-pca_clean$x[,1:100]
pcs_teste<-predict(pca_clean,newdata=te_s)[,1:100]
print(paste("Split Rigoroso -> Treino:",nrow(treino),"| Teste:",nrow(teste)))

# ============================================================================
# 1. RANDOM FOREST
# ============================================================================
set.seed(335705)
rf_model<-randomForest(tipo~.,data=treino,ntree=500,importance=TRUE)
pred_rf<-predict(rf_model,teste)
pred_rf_prob<-predict(rf_model,teste,type="prob")
conf_rf<-confusionMatrix(pred_rf,teste$tipo)
print("========== Random Forest ==========")
print(paste("Acurácia Global:",round(conf_rf$overall['Accuracy'],5)))
print(cbind(conf_rf$byClass[,c("Precision","Sensitivity","F1")], N=table(teste$tipo)))

# Importância (RF)
imp_rf<-data.frame(gene=rownames(importance(rf_model)),importancia=importance(rf_model)[,"MeanDecreaseGini"])%>%arrange(desc(importancia))%>%head(20)
p_imp_rf<-ggplot(imp_rf,aes(x=reorder(gene,importancia),y=importancia))+geom_col(fill="#E74C3C",alpha=0.85,width=0.7)+coord_flip()+theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold",size=14,hjust=0.5),axis.title.x=element_text(face="italic",size=11),axis.title.y=element_text(face="bold",size=11),axis.text.y=element_text(face="italic",size=10),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank())+labs(title="20 Genes de Maior Importância - Random Forest",x="Gene",y=expression(paste("Importância (",italic("Mean Decrease Gini"),")")))+scale_y_continuous(expand=expansion(mult=c(0,0.05)))
print(p_imp_rf)

# Matriz Confusão (RF)
melt_rf<-melt(conf_rf$table)
colnames(melt_rf)<-c("Real","Predito","Freq")
melt_rf<-melt_rf%>%left_join(melt_rf%>%group_by(Predito)%>%summarise(Total=sum(Freq)),by="Predito")%>%mutate(Prop=Freq/Total*100)
p_conf_rf<-ggplot(melt_rf,aes(x=Predito,y=Real,fill=Prop))+geom_tile(color="white")+scale_fill_gradient2(low="#FFF5F0",mid="#FB6A4A",high="#67000D",midpoint=50,limits=c(0,100),name="%")+geom_text(aes(label=Freq),size=3,fontface="bold")+theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold",size=14,hjust=0.5),axis.text.x=element_text(angle=45,hjust=1,face="italic"),axis.text.y=element_text(face="italic"),panel.grid=element_blank())+labs(title="Matriz de Confusão - Random Forest",x="Predito",y="Real")
print(p_conf_rf)

# AUC (RF)
auc_rf<-sapply(levels(teste$tipo),function(t) auc(ifelse(teste$tipo==t,1,0),pred_rf_prob[,t]))
print(data.frame(Tipo=names(auc_rf),AUC=auc_rf))
cat(paste("AUC Média RF:",round(mean(auc_rf),4),"\n"))

# --- RF (100 PCs) ---
set.seed(335705)
rf_pca<-randomForest(tipo~.,data=data.frame(pcs_treino,tipo=treino$tipo),ntree=500)
pred_rf_pca<-predict(rf_pca,data.frame(pcs_teste,tipo=teste$tipo))
print(paste("Acurácia RF (100 PCs):",round(confusionMatrix(pred_rf_pca,teste$tipo)$overall['Accuracy'],5)))

# ============================================================================
# 2. SVM
# ============================================================================
set.seed(335705)
svm_model<-svm(tipo~.,data=treino,kernel="radial",cost=1,gamma=0.001,probability=TRUE)
pred_svm<-predict(svm_model,teste)
pred_svm_prob<-attr(predict(svm_model,teste,probability=TRUE),"probabilities")
conf_svm<-confusionMatrix(pred_svm,teste$tipo)
print("========== SVM ==========")
print(paste("Acurácia Global:",round(conf_svm$overall['Accuracy'],5)))
print(cbind(conf_svm$byClass[,c("Precision","Sensitivity","F1")], N=table(teste$tipo)))

# Matriz Confusão (SVM)
melt_svm<-melt(conf_svm$table)
colnames(melt_svm)<-c("Real","Predito","Freq")
melt_svm<-melt_svm%>%left_join(melt_svm%>%group_by(Predito)%>%summarise(Total=sum(Freq)),by="Predito")%>%mutate(Prop=Freq/Total*100)
p_conf_svm<-ggplot(melt_svm,aes(x=Predito,y=Real,fill=Prop))+geom_tile(color="white")+scale_fill_gradient2(low="#F2F0F7",mid="#9E9AC8",high="#54278F",midpoint=50,limits=c(0,100),name="%")+geom_text(aes(label=Freq),size=3,fontface="bold")+theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold",size=14,hjust=0.5),axis.text.x=element_text(angle=45,hjust=1,face="italic"),axis.text.y=element_text(face="italic"),panel.grid=element_blank())+labs(title="Matriz de Confusão - SVM",x="Predito",y="Real")
print(p_conf_svm)

# AUC (SVM)
auc_svm<-sapply(levels(teste$tipo),function(t) auc(ifelse(teste$tipo==t,1,0),pred_svm_prob[,t]))
print(data.frame(Tipo=names(auc_svm),AUC=auc_svm))
cat(paste("AUC Média SVM:",round(mean(auc_svm),4),"\n"))

# --- SVM (100 PCs) ---
set.seed(335705)
svm_pca<-svm(tipo~.,data=data.frame(pcs_treino,tipo=treino$tipo),kernel="radial",cost=10)
pred_svm_pca<-predict(svm_pca,data.frame(pcs_teste,tipo=teste$tipo))
print(paste("Acurácia SVM (100 PCs):",round(confusionMatrix(pred_svm_pca,teste$tipo)$overall['Accuracy'],5)))

# ============================================================================
# 3. XGBOOST
# ============================================================================
set.seed(335705)
lab_tr<-as.numeric(treino$tipo)-1
lab_te<-as.numeric(teste$tipo)-1
dtr<-xgb.DMatrix(data=as.matrix(treino[,-ncol(treino)]),label=lab_tr)
dte<-xgb.DMatrix(data=as.matrix(teste[,-ncol(teste)]),label=lab_te)
xgb_model<-xgboost(data=dtr,nrounds=100,max_depth=6,eta=0.3,objective="multi:softprob",num_class=length(unique(lab_tr)),verbose=0)
prob_xgb<-predict(xgb_model,dte,reshape=TRUE)
pred_xgb_idx<-max.col(prob_xgb)-1
pred_xgb<-factor(pred_xgb_idx,levels=0:(length(unique(lab_tr))-1),labels=levels(teste$tipo))
conf_xgb<-confusionMatrix(pred_xgb,teste$tipo)
print("========== XGBoost ==========")
print(paste("Acurácia Global:",round(conf_xgb$overall['Accuracy'],5)))
print(cbind(conf_xgb$byClass[,c("Precision","Sensitivity","F1")], N=table(teste$tipo)))

# Importância (XGB)
imp_xgb<-xgb.importance(model=xgb_model)
p_imp_xgb<-ggplot(head(imp_xgb,20),aes(x=reorder(Feature,Gain),y=Gain))+geom_col(fill="#27AE60",alpha=0.85,width=0.7)+coord_flip()+theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold",size=14,hjust=0.5),axis.title.x=element_text(face="italic",size=11),axis.title.y=element_text(face="bold",size=11),axis.text.y=element_text(face="italic",size=10),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank())+labs(title="20 Genes de Maior Importância - XGBoost",x="Gene",y=expression(paste("Importância (",italic("Gain"),")")))+scale_y_continuous(expand=expansion(mult=c(0,0.05)))
print(p_imp_xgb)

# Matriz Confusão (XGB)
melt_xgb<-melt(conf_xgb$table)
colnames(melt_xgb)<-c("Real","Predito","Freq")
melt_xgb<-melt_xgb%>%left_join(melt_xgb%>%group_by(Predito)%>%summarise(Total=sum(Freq)),by="Predito")%>%mutate(Prop=Freq/Total*100)
p_conf_xgb<-ggplot(melt_xgb,aes(x=Predito,y=Real,fill=Prop))+geom_tile(color="white")+scale_fill_gradient2(low="#F7FCF5",mid="#74C476",high="#00441B",midpoint=50,limits=c(0,100),name="%")+geom_text(aes(label=Freq),size=3,fontface="bold")+theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold",size=14,hjust=0.5),axis.text.x=element_text(angle=45,hjust=1,face="italic"),axis.text.y=element_text(face="italic"),panel.grid=element_blank())+labs(title="Matriz de Confusão - XGBoost",x="Predito",y="Real")
print(p_conf_xgb)

# AUC (XGB)
auc_xgb<-sapply(1:length(levels(teste$tipo)),function(i) auc(ifelse(lab_te==(i-1),1,0),prob_xgb[,i]))
names(auc_xgb)<-levels(teste$tipo)
print(data.frame(Tipo=names(auc_xgb),AUC=auc_xgb))
cat(paste("AUC Média XGB:",round(mean(auc_xgb),4),"\n"))

# --- XGB (100 PCs) ---
set.seed(335705)
dtr_pca<-xgb.DMatrix(data=pcs_treino,label=lab_tr)
dte_pca<-xgb.DMatrix(data=pcs_teste,label=lab_te)
xgb_pca<-xgboost(data=dtr_pca,nrounds=100,max_depth=6,eta=0.3,objective="multi:softprob",num_class=length(unique(lab_tr)),verbose=0)
prob_xgb_pca<-predict(xgb_pca,dte_pca,reshape=TRUE)
pred_xgb_pca<-factor(max.col(prob_xgb_pca)-1,levels=0:(length(unique(lab_tr))-1),labels=levels(teste$tipo))
print(paste("Acurácia XGB (100 PCs):",round(confusionMatrix(pred_xgb_pca,teste$tipo)$overall['Accuracy'],5)))

# ============================================================================
# 4. NAIVE BAYES
# ============================================================================
set.seed(335705)
nb_model<-naiveBayes(tipo~.,data=treino)
pred_nb<-predict(nb_model,teste)
pred_nb_prob<-predict(nb_model,teste,type="raw")
conf_nb<-confusionMatrix(pred_nb,teste$tipo)
print("========== Naive Bayes ==========")
print(paste("Acurácia Global:",round(conf_nb$overall['Accuracy'],5)))
print(cbind(conf_nb$byClass[,c("Precision","Sensitivity","F1")], N=table(teste$tipo)))

# Matriz Confusão (NB)
melt_nb<-melt(conf_nb$table)
colnames(melt_nb)<-c("Real","Predito","Freq")
melt_nb<-melt_nb%>%left_join(melt_nb%>%group_by(Predito)%>%summarise(Total=sum(Freq)),by="Predito")%>%mutate(Prop=Freq/Total*100)
p_conf_nb<-ggplot(melt_nb,aes(x=Predito,y=Real,fill=Prop))+geom_tile(color="white")+scale_fill_gradient2(low="#E0F2F1",mid="#4DB6AC",high="#004D40",midpoint=50,limits=c(0,100),name="%")+geom_text(aes(label=Freq),size=3,fontface="bold")+theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold",size=14,hjust=0.5),axis.text.x=element_text(angle=45,hjust=1,face="italic"),axis.text.y=element_text(face="italic"),panel.grid=element_blank())+labs(title="Matriz de Confusão - Naive Bayes",x="Predito",y="Real")
print(p_conf_nb)

# AUC (NB)
auc_nb<-sapply(levels(teste$tipo),function(t) auc(ifelse(teste$tipo==t,1,0),pred_nb_prob[,t]))
print(data.frame(Tipo=names(auc_nb),AUC=auc_nb))
cat(paste("AUC Média NB:",round(mean(auc_nb),4),"\n"))

# --- NB (100 PCs) ---
set.seed(335705)
nb_pca<-naiveBayes(tipo~.,data=data.frame(pcs_treino,tipo=treino$tipo))
pred_nb_pca<-predict(nb_pca,data.frame(pcs_teste,tipo=teste$tipo))
print(paste("Acurácia NB (100 PCs):",round(confusionMatrix(pred_nb_pca,teste$tipo)$overall['Accuracy'],5)))

# ============================================================================
# 5. KNN
# ============================================================================
set.seed(335705)
knn_model<-knn3(tipo~.,data=treino,k=5)

# Predição de Classes (para Matriz de Confusão)
pred_knn<-predict(knn_model,teste,type="class")
conf_knn<-confusionMatrix(pred_knn,teste$tipo)
print("========== KNN ==========")
print(paste("Acurácia Global:",round(conf_knn$overall['Accuracy'],5)))
print(cbind(conf_knn$byClass[,c("Precision","Sensitivity","F1")], N=table(teste$tipo)))

# Matriz Confusão (KNN)
melt_knn<-melt(conf_knn$table)
colnames(melt_knn)<-c("Real","Predito","Freq")
melt_knn<-melt_knn%>%left_join(melt_knn%>%group_by(Predito)%>%summarise(Total=sum(Freq)),by="Predito")%>%mutate(Prop=Freq/Total*100)
p_conf_knn<-ggplot(melt_knn,aes(x=Predito,y=Real,fill=Prop))+geom_tile(color="white")+scale_fill_gradient2(low="#FFF5EB",mid="#FDBE85",high="#A63603",midpoint=50,limits=c(0,100),name="%")+geom_text(aes(label=Freq),size=3,fontface="bold")+theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold",size=14,hjust=0.5),axis.text.x=element_text(angle=45,hjust=1,face="italic"),axis.text.y=element_text(face="italic"),panel.grid=element_blank())+labs(title="Matriz de Confusão - KNN",x="Predito",y="Real")
print(p_conf_knn)

# AUC (KNN) - Usando o mesmo modelo
pred_knn_prob<-predict(knn_model,teste,type="prob")
auc_knn<-sapply(levels(teste$tipo),function(t) auc(ifelse(teste$tipo==t,1,0),pred_knn_prob[,t]))
print(data.frame(Tipo=names(auc_knn),AUC=auc_knn))
cat(paste("AUC Média KNN:",round(mean(auc_knn),4),"\n"))

# --- KNN (100 PCs) ---
set.seed(335705)
# Treinando novo modelo nos PCs
knn_pca_model<-knn3(x=pcs_treino,y=treino$tipo,k=5)
pred_knn_pca<-predict(knn_pca_model,pcs_teste,type="class")
print(paste("Acurácia KNN (100 PCs):",round(confusionMatrix(pred_knn_pca,teste$tipo)$overall['Accuracy'],5)))

# ============================================================================
# COMPARAÇÕES FINAIS
# ============================================================================
df_acc<-data.frame(Metodo=c("Random Forest","SVM","XGBoost","Naive Bayes","KNN"),Original=c(conf_rf$overall['Accuracy'],conf_svm$overall['Accuracy'],conf_xgb$overall['Accuracy'],conf_nb$overall['Accuracy'],conf_knn$overall['Accuracy']))
p_acc<-ggplot(df_acc,aes(x=reorder(Metodo,Original),y=Original))+geom_col(fill="#2980B9",alpha=0.85,width=0.7)+coord_flip()+theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold",size=14,hjust=0.5),axis.title.x=element_text(face="bold",size=11),axis.title.y=element_text(face="bold",size=11),axis.text.y=element_text(face="italic",size=10),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank())+labs(title="Comparação de Acurácia por Método",x="Método",y="Acurácia Global")+scale_y_continuous(expand=expansion(mult=c(0,0.05)))+geom_text(aes(label=round(Original,4)),hjust=-0.1,size=3.5)
print(p_acc)
print(df_acc)

# Validação Cruzada
control<-trainControl(method="cv",number=5)
set.seed(335705)
df_cv<-treino
model_rf_cv<-train(tipo~.,data=df_cv,method="rf",trControl=control,ntree=100)
model_svm_cv<-train(tipo~.,data=df_cv,method="svmRadial",trControl=control)
model_xgb_cv<-train(tipo~.,data=df_cv,method="xgbTree",trControl=control,verbose=FALSE)
model_nb_cv<-train(tipo~.,data=df_cv,method="naive_bayes",trControl=control)
model_knn_cv<-train(tipo~.,data=df_cv,method="knn",trControl=control)
resumps<-resamples(list(RF=model_rf_cv,SVM=model_svm_cv,XGB=model_xgb_cv,NB=model_nb_cv,KNN=model_knn_cv))
print(summary(resumps))
print(dotplot(resumps,metric="Accuracy",main="CV Accuracy (5-Fold)"))