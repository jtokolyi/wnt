# Step 1: Load necessary libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

count_matrix <- read.csv("/storage/RNASeq2025/06_salmon/alp.gene.counts.matrix", header = TRUE, row.names = 1, sep = "\t")
count_matrix <- count_matrix[,order(names(count_matrix))]

metadata <- data.frame(
    condition = factor(c("ALP","ALP","ALP","ALP",
                         "DMSO","DMSO","DMSO","DMSO")),
    row.names = colnames(count_matrix))

dds <- DESeqDataSetFromMatrix(countData = round(count_matrix), colData = metadata,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)>=10) >=3,]

dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)

sig_genes <- rownames(res_lrt)[which(res_lrt$padj < 0.05)]

rld <- rlog(dds, blind = FALSE)
rld_sig <- rld[sig_genes, ]

plotPCA(rld_sig, intgroup = "condition")

#############################################################################################
trin <- read.csv("/storage/RNASeq2025/07_trinotate/alp_cdhit/Trinotate_alp_cdhit.xls",sep="\t")
wnt.trin <- trin[grep("GO:0060070|GO:0060828|GO:0090090|GO:0090263", trin$gene_ontology_BLASTP),]
wnt.trin$sprot_gene <- sapply(strsplit(sapply(strsplit(wnt.trin$sprot_Top_BLASTP_hit,split="^",fixed=T),"[",1),split="_"),"[",1)
wnt.trin <- wnt.trin[order(wnt.trin$transcript_id),]

wnt_hynames <- read_excel("/storage/RNASeq2025/08_deseq2/wnt_hydra_gene_names_to_find.xlsx")
wnt_hynames <- wnt_hynames[order(wnt_hynames$transcript_id),]
wnt.trin$hydra.names <-wnt_hynames$short_names
    
blast <- read.delim("/storage/RNASeq2025/07_trinotate/alp_cdhit/alp_cdhit_hydra_proteome_blast.tsv",
  header = FALSE,  sep = "\t", stringsAsFactors = FALSE)
colnames(blast) <- c("prot_id","subject_id","description","evalue","bitscore")
merged_df <- merge(wnt.trin,blast[which(!duplicated(blast$prot_id)),],by = "prot_id",all.x = TRUE)
merged_df$Gene_Name_Hydra <- sub(".*\\bGN=([^ ]+).*", "\\1", merged_df$description)
write.table(merged_df, file="/storage/RNASeq2025/07_trinotate/alp_cdhit/wnt_genes_table_with_hydra_names.csv",col.names=T, row.names=F, sep="\t", quote=F)

wnt_count_matrix <- count_matrix[row.names(count_matrix)%in%wnt.trin$X.gene_id,]
dds_wnt <- DESeqDataSetFromMatrix(countData = round(wnt_count_matrix), colData = metadata,
                              design = ~ condition)
dds_wnt$condition <- relevel(dds_wnt$condition, "DMSO")
dds_wnt <- dds_wnt[rowSums(counts(dds_wnt)>=10) >=3,]
res_wnt <- results(DESeq(dds_wnt))%>%data.frame()%>%
  tibble::rownames_to_column(var="geneID") %>% 
  left_join(wnt.trin[c("X.gene_id", "hydra.names")],join_by("geneID"=="X.gene_id"))%>%
  group_by(geneID) %>% slice(1) %>% ungroup()

res_wnt$gene_type <- ifelse(res_wnt$padj>0.01 | abs(res_wnt$log2FoldChange)<1 , "Not affected",
                         ifelse(res_wnt$log2FoldChange < -1 , "Down-regulated", "Up-regulated"))
cols <- c("Up-regulated" = "#CA0020", "Down-regulated" = "#0571B0", "Not affected" = "white")
sizes <- c("Up-regulated" = 2, "Down-regulated" = 2, "Not affected" = 1)
alphas <- c("Up-regulated" = 0.5, "Down-regulated" = 0.5, "Not affected" = 0.5)

vplot <- ggplot(res_wnt[-which(is.na(res_wnt$padj)),], aes(y=-log10(padj), x=log2FoldChange,
                             size=gene_type, fill=gene_type, color=gene_type))+
  geom_point(shape=21, col="black", alpha=0.5)+theme_bw()+xlim(-5.5,5.5)+ylim(-5,45)+
  ggrepel::geom_label_repel(
    data = res_wnt %>% filter(padj < 0.01 & abs(log2FoldChange) > 1),
    aes(label = hydra.names,color=gene_type),fill="white", 
    size = 3,
    max.overlaps = 20,
    min.segment.length = 0,
    box.padding=0.4,
    point.padding=0.5,
    segment.curvature=0.1,
    nudge_x = res_wnt %>% 
    filter(padj < 0.01 & abs(log2FoldChange) > 1) %>% 
    pull(log2FoldChange) %>% 
    sapply(function(x) ifelse(x > 0, 0.25, 0)),
    nudge_y=1,
    force=5,
    max.iter=5000,
    direction = "both",
    segment.color = "grey50"
  ) +
  labs(title="cWnt DEGs following ALP treatment")+
  geom_vline(xintercept=c(-1,1), linetype="dashed",alpha=0.5)+
  geom_hline(yintercept = -log10(0.01), linetype="dashed",alpha=0.5)+
  scale_size_manual(name="", values=sizes)+
  scale_fill_manual(name="", values=cols)+
  scale_color_manual(name="", values=cols)+
  xlab(bquote(~Log[2] ~ "fold change"))+ylab(bquote(~-Log[10] ~ P[adj]))

ragg::agg_png(file="/storage/RNASeq2025/08_deseq2/wnt_ALP.png", height=1500, width=2500, units="px", scaling=0.9, res=300)
vplot
dev.off()

#############################################################################################
library(topGO); library(rrvgo); library(GOSemSim); library(GO.db); library(igraph); library(ggraph); library(tidygraph); library(ggplot2)

gene2go <- trin %>%
  dplyr::select(gene_id = X.gene_id, GO_raw = gene_ontology_BLASTP) %>%
  filter(!is.na(GO_raw) & GO_raw != "") %>%
  mutate(GO_list = strsplit(as.character(GO_raw), "`")) %>%
  unnest(GO_list) %>%
  mutate(GO = sub("\\^.*", "", GO_list)) %>%
    dplyr::select(gene_id, GO) %>%
    filter(GO!=".")

colData(dds)$condition <- relevel(colData(dds)$condition, "DMSO")

dds <- DESeq(dds)
res <- results(dds)

up_genes <- rownames(res[which(res$padj < 0.05 & res$log2FoldChange > 0), ])
down_genes <- rownames(res[which(res$padj < 0.05 & res$log2FoldChange < -0),])

up.geneList <- factor(as.integer(rownames(res) %in% up_genes))
names(up.geneList) <- rownames(res)

up.GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = up.geneList,
              annot = annFUN.gene2GO,
              gene2GO = split(gene2go$GO, gene2go$gene_id))
resultFisher <- runTest(up.GOdata, algorithm = "classic", statistic = "fisher")
up.allGO <- usedGO(up.GOdata)
up.allResults <- GenTable(up.GOdata,classicFisher = resultFisher, orderBy = "classicFisher",
                       topNodes = length(up.allGO))
up.allResults$classicFisher <- as.numeric(as.character(gsub("< ", "", up.allResults$classicFisher)))
up.enrich <- up.allResults %>% filter(classicFisher < 0.001) 

go_terms <- up.enrich$GO.ID
scores <- -log10(as.numeric(up.enrich$classicFisher))
names(scores) <- go_terms
simMatrix <- calculateSimMatrix(go_terms, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")
sim_df <- as.data.frame(as.table(simMatrix))
colnames(sim_df) <- c("from", "to", "weight")
sim_df <- sim_df[sim_df$from != sim_df$to & sim_df$weight > 0.5, ]
g <- graph_from_data_frame(sim_df, directed = FALSE)
term_cluster <- setNames(reducedTerms$parentTerm, reducedTerms$go)
term_label   <- setNames(reducedTerms$term, reducedTerms$go)
term_score   <- setNames(reducedTerms$score, reducedTerms$go)
V(g)$cluster <- term_cluster[V(g)$name]
V(g)$label   <- term_label[V(g)$name]
V(g)$score   <- term_score[V(g)$name]
cluster_sizes <- table(V(g)$cluster)
valid_clusters <- names(cluster_sizes[cluster_sizes > 1])
g <- induced_subgraph(g, vids = V(g)[cluster %in% valid_clusters])
g_up <- g

nodes <- data.frame(
  id     = names(term_label),
  label  = unname(term_label),
  score  = unname(term_score),
  cluster = unname(term_cluster)
)

write.csv(sim_df, "edges_GO_up.csv", row.names = FALSE)
write.csv(nodes, "nodes_GO_up.csv", row.names = FALSE)



