

#####loadpackages-------------------------


library(tidyverse)
library(limma)
library(edgeR)
library(fgsea)
library(enrichR)
library(msigdbr)
library(here)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# =========================
# DIRECTORIES
# =========================

# LIMMA results from RNA-seq script
limma_dir <- here("results/beV_RNAseq_monocytes_limma_01")

# Output directory for downstream analysis (e.g. fgsea, plots)
outdir <- here("results/beV_RNAseq_monocytes_enrichment_02")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# =========================
# LOAD LIMMA RESULTS
# =========================

limmaRes <- readRDS(
  file.path(limma_dir, "limmaRes.rds")
)


# =========================
# LOAD BATCH-CORRECTED EXPRESSION
# =========================

expr_corrected <- readRDS(
  file.path(limma_dir, "expr_corrected.rds")
)


# =========================
# LOAD METADATA
# =========================

metadata <- readRDS(
  file.path(limma_dir, "metadata.rds")
)



#####Load MSigdbr------------
msig <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(x = msig$gene_symbol, f = msig$gs_name)
head(pathways)


fgsea_all <- list()

for (coef_current in unique(limmaRes$coef)) {
  
  ranks <- limmaRes |>
    filter(coef == coef_current)
  
  rank_vector <- ranks$logFC
  names(rank_vector) <- ranks$Gene.name
  rank_vector <- sort(rank_vector, decreasing = TRUE)
  
  fgseaRes <- fgsea(
    pathways = pathways,
    stats = rank_vector,
    minSize = 5,
    maxSize = 500
  ) |>
    as_tibble() |>
    mutate(coef = coef_current)
  
  fgsea_all[[coef_current]] <- fgseaRes
}

fgsea_all_combined <- bind_rows(fgsea_all)

show(fgsea_all_combined)


topPathways <- fgsea_all_combined |>
  filter(coef == "OMVvsHp") |>
  filter(padj < 0.05) |>
  slice_head(n = 50) |>
  slice_max(order_by = NES, n = 20) |>
  bind_rows(
    fgsea_all_combined |>
      filter(coef == "OMVvsHp") |>
      filter(padj < 0.05) |>
      slice_head(n = 50) |>
      slice_min(order_by = NES, n = 20)
  ) |>
  ungroup()

selected_pathways <- unique(topPathways$pathway)


topPathways <- fgsea_all_combined |>
  filter(pathway %in% selected_pathways)


topPathways <- topPathways |>
  filter(coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  arrange(coef, NES)

ggfgsea_combined <- ggplot(topPathways,
                           aes(x = coef, y = pathway, size = -log10(padj), color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(
    low = "#0571b0",
    mid = "#F7F7F7",
    high = "#ca0020",
    midpoint = 0
  ) +
  scale_size_continuous(range = c(3,5))+
  labs(
    title = "FGSEA: Top Pathways (Up & Down)",
    x = "Coefficient",
    y = "Pathway",
    size = "-log10(padj)",
    color = "NES"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  )

ggsave(
  filename = file.path(outdir, "enrichment_fgsea_combined_new.png"),
  plot = ggfgsea_combined,
  width = 12,
  height = 10,
  dpi = 300
)


#######filteringgenesinpathways-----------------


genes <- topPathways |>
  filter(pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") |>
  pull(leadingEdge) |>
  str_remove_all('c\\(|\\)|"|\\n') |>
  str_split(",\\s*") |>
  unlist()
head(genes)

goi_pathway <- limmaRes |>
  filter(coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  filter(Gene.name%in%genes)


a <- ggplot(goi_pathway, aes(x = coef, y = Gene.name)) +
  geom_point(aes(size = -log10(adj.P.Val), color = logFC)) +
  scale_color_gradient2(
    low = "#0571b0",
    mid = "#F7F7F7",
    high = "#ca0020",
    midpoint = 0
  )+
  theme_minimal() +
  labs(
    x = "Coefficient",
    y = "Gene",
    size = "-log10(adj.P.Val)",
    color = "logFC"
  )

ggsave(
  filename = file.path(outdir, "enrichment_fgsea_Epithilial_Mesenchymal_Pathway.png"),
  plot = a,
  width = 12,
  height = 10,
  dpi = 300
)

# -----------plotting expression data for these genes----

dat.list <- list()
for(gg in genes){
  dat.list[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat <- bind_rows(dat.list, .id = "Gene.name")
dat$treatment <- factor(dat$treatment, levels = c("mock", "Hp", "OMV"))

dat <- dat |>
  mutate(
    sample_treatment_order = factor(
      paste0(donor, ".", treatment),
      levels = c(
        "D1.mock", "D2.mock", "D3.mock",
        "D1.Hp",   "D2.Hp",   "D3.Hp",
        "D1.OMV",  "D2.OMV",  "D3.OMV"
      )
    )
  )


b <- ggplot(dat, aes(x = donor, y = Gene.name, fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#0571b0",
    mid = "#F7F7F7",
    high = "#ca0020",
    midpoint = 0
  ) +
  facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Gene",
    fill = "Scaled expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(1, "lines")  # adds space between mock / Hp / OMV
  )

ggsave(
  filename = file.path(outdir, "Expression_Mesenchymal_Pathway.png"),
  plot = b,
  width = 12,
  height = 10,
  dpi = 300
)


#------------------------------------------ enrichR on upregulated genes------, 

ENRICHR.DBS <- c(
                 "Reactome_2022",
                 "GO_Molecular_Function_2023",
                 "GO_Biological_Process_2023",
                 "CellMarker_2024")


coef_list <- unique(limmaRes$coef)

enrich_upregulated <- list()

for (coef_current in coef_list) {
  
  tt <- limmaRes |>
    filter(coef == coef_current)
  
  goi <- tt |>
    filter(adj.P.Val < 0.05, logFC > 1) |>
    pull(Gene.name)
  
  if (length(goi) < 5) next
  
  results <- enrichr(goi, ENRICHR.DBS)
  
  combined <- do.call(rbind, lapply(names(results), function(db) {
    df <- results[[db]]
    df$Database <- db
    df$Coefficient <- coef_current
    return(df)
  }))
  
  enrich_upregulated[[coef_current]] <- combined
}

final_combined_upregulated <- bind_rows(enrich_upregulated)


#select the pathways based on significance 
combined_sig <- final_combined_upregulated |>
  filter(Adjusted.P.value < 0.05) |>
  arrange(Adjusted.P.value)

#selecting top 5 enriched pathways per data base
top_per_db <- combined_sig |>
  filter(Coefficient %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  group_by(Database, Coefficient) |>
  slice_head(n = 5) |>
  ungroup() 


#plotting for each database
p <- ggplot(top_per_db,
            aes(x = Coefficient,                   # numeric x-axis
                y = reorder(Term, Coefficient),   # pathways on y-axis, ordered
                size = -log10(Adjusted.P.value),
                color = log2(Odds.Ratio))) +
  geom_point(alpha = 0.8) +                 # semi-transparent points
  scale_size_continuous(range = c(3,6)) + 
  scale_color_gradient2(
    low = "#F7F7F7",
    high = "#ca0020",
  )+
  labs(
    title = "Top Pathways per Database",
    x = "Coefficient",
    y = "Pathway / Term",
    size = "-log10(Adjusted.P.value)",
    color = "log2(Odds.Ratio)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),   # pathway names bigger
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25)
  )

ggsave(
  filename = file.path(outdir, "Pathways_enrichment_upregulatedgenes_databases.png"),
  plot = p,
  width = 12,
  height = 10,
  dpi = 300
)


##--------------------Downregulated genes------enrichR_____

coef_list <- unique(limmaRes$coef)

enrich_downregulated <- list()

for (coef_current in coef_list) {
  
  tt <- limmaRes |>
    filter(coef == coef_current)
  
  goi <- tt |>
    filter(adj.P.Val < 0.05, logFC < -1) |>
    pull(Gene.name)
  
  if (length(goi) < 5) next
  
  results <- enrichr(goi, ENRICHR.DBS)
  
  combined <- do.call(rbind, lapply(names(results), function(db) {
    df <- results[[db]]
    df$Database <- db
    df$Coefficient <- coef_current
    return(df)
  }))
  
  enrich_downregulated[[coef_current]] <- combined
}

final_combined_downregulated <- bind_rows(enrich_downregulated)

#select the arrange based on significance 
combined_sig <- final_combined_downregulated |>
  filter(Adjusted.P.value < 0.05) |>
  arrange(Adjusted.P.value)

#selecting top 5 enriched pathways per data base
top_per_db <- combined_sig |>
  filter(Coefficient %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  group_by(Database, Coefficient) |>
  slice_head(n = 5) |>
  ungroup() 

#plotting for each database
f <- ggplot(top_per_db,
            aes(x = Coefficient,                   # numeric x-axis
                y = reorder(Term, Coefficient),   # pathways on y-axis, ordered
                size = -log10(Adjusted.P.value),
                color = log2(Odds.Ratio))) +
  geom_point(alpha = 0.8) +                 # semi-transparent points
  scale_size_continuous(range = c(3,6)) +
  scale_color_gradient2(
    low = "#F7F7F7",
    high = "#ca0020",
  )+
  labs(
    title = "Top Pathways per Database",
    x = "Coefficient",
    y = "Pathway / Term",
    size = "-log10(Adjusted P Value)",
    color = "log2(Odds Ratio)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),   # pathway names bigger
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25)
  )
ggsave(
  filename = file.path(outdir, "Pathways_enrichment_downregulatedgenes_databases.png"),
  plot = f,
  width = 12,
  height = 10,
  dpi = 300
)


#list of genes regulating chemokines and their expression----------

genes_by_coef <- final_combined_upregulated |>
  filter(str_detect(Term, "1990869")) |>
  dplyr::select(Database, Coefficient, Genes)


gene_list <- str_split(genes_by_coef$Genes, ";") |>
  unlist() |>
  unique()


goi_pathway_chemokine <- limmaRes |>
  filter(Gene.name %in% gene_list) |>
  filter(coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))


i <- ggplot(goi_pathway_chemokine, aes(x = coef, y = Gene.name)) +
  geom_point(aes(size = -log10(adj.P.Val), color = logFC)) +
  scale_color_gradient2(
    low = "#0571b0",
    mid = "#F7F7F7",
    high = "#ca0020",
    midpoint = 0
  )+
  theme_minimal() +
  labs(
    x = "Coefficient",
    y = "Gene",
    size = "-log10(adj.P.Val)",
    color = "logFC"
  )

ggsave(
  filename = file.path(outdir, "Expression_of_genes_cellular_response_to_chemokine.png"),
  plot = i,
  width = 12,
  height = 10,
  dpi = 300
)

#------------------------expression list of the genes selected batch corrected--------

dat.list_chemo <- list()
for(gg in gene_list){
  dat.list_chemo[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat_chemo <- bind_rows(dat.list_chemo, .id = "Gene.name")
dat_chemo$treatment <- factor(dat_chemo$treatment, levels = c("mock", "Hp", "OMV"))



dat_chemo <- dat_chemo |>
  mutate(
    sample_treatment_order = factor(
      paste0(donor, ".", treatment),
      levels = c(
        "D1.mock", "D2.mock", "D3.mock",
        "D1.Hp",   "D2.Hp",   "D3.Hp",
        "D1.OMV",  "D2.OMV",  "D3.OMV"
      )
    )
  )


k <- ggplot(dat_chemo, aes(x = donor, y = Gene.name, fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#0571b0",
    mid = "#F7F7F7",
    high = "#ca0020",
    midpoint = 0
  ) +
  facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Gene",
    fill = "Scaled expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(1, "lines")  # adds space between mock / Hp / OMV
  )

ggsave(
  filename = file.path(outdir, "Expression of genes regulating cellular response to chemokines batch corrected.png"),
  plot = k,
  width = 12,
  height = 10,
  dpi = 300
)


#####fgsea_with_additional_pathways_----

# Hallmark (clean signals)
msig_h <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

# KEGG (mechanistic pathways)
msig_kegg <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:KEGG_LEGACY"
)

# GO Biological Process (broad + sensitive)
msig_go <- msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "BP"
)

msig_combined <- bind_rows(msig_h, msig_kegg, msig_go)
pathways <- split(x = msig_combined$gene_symbol, f = msig_combined$gs_name)
head(pathways)

keywords <- c("ACTIN", "CYTOSKELETON", "TLR", "TOLL", "INTEGRIN")

pathways_filtered <- pathways[
  grepl(paste(keywords, collapse = "|"), names(pathways), ignore.case = TRUE)
]

names(pathways_filtered)

fgsea_all <- list()

for (coef_current in unique(limmaRes$coef)) {
  
  ranks <- limmaRes |>
    filter(coef == coef_current)
  
  rank_vector <- ranks$logFC
  names(rank_vector) <- ranks$Gene.name
  rank_vector <- sort(rank_vector, decreasing = TRUE)
  
  fgseaRes <- fgsea(
    pathways = pathways_filtered,
    stats = rank_vector,
    minSize = 5,
    maxSize = 500
    #nperm = 10000
  ) |>
    as_tibble() |>
    mutate(coef = coef_current)
  
  fgsea_all[[coef_current]] <- fgseaRes
}

fgsea_all_combined <- bind_rows(fgsea_all)

show(fgsea_all_combined)
topPathways <- fgsea_all_combined |>
  group_by(coef) |>
  filter(padj < 0.05) |>
  slice_head(n=50) |>
  slice_max(order_by = NES, n = 10) |>   
  bind_rows(
    fgsea_all_combined |>
      group_by(coef) |>
      filter(padj < 0.05) |>
      slice_head(n=50) |>
      slice_min(order_by = NES, n = 10)  
  ) |>
  ungroup()

topPathways <- topPathways |>
  filter(coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  arrange(coef, NES)

ggfgsea_combined <- ggplot(topPathways,
                           aes(x = coef, y = pathway, size = -log10(padj), color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(
    low = "#0571b0",
    mid = "#F7F7F7",
    high = "#ca0020",
    midpoint = 0
  ) +
  labs(
    title = "FGSEA: Top Pathways (Up & Down) for selected pathways",
    x = "Coefficient",
    y = "Pathway",
    size = "-log10(padj)",
    color = "NES"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  )

ggsave(
  filename = file.path(outdir, "Genes regulating cytoskeleton, TLR, actin.png"),
  plot = ggfgsea_combined,
  width = 12,
  height = 10,
  dpi = 300
)


###--------------filter out genes in the KEGG_TOLL_pathway--------

genes_toll <- topPathways |>
  filter(pathway == "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY") |>
  pull(leadingEdge) |>
  str_remove_all('c\\(|\\)|"|\\n') |>
  str_split(",\\s*") |>
  unlist()
show(genes_toll)

goi_pathway <- limmaRes |>
  filter(coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  filter(Gene.name%in%genes_toll)


a <- ggplot(goi_pathway, aes(x = coef, y = Gene.name)) +
  geom_point(aes(size = -log10(adj.P.Val), color = logFC)) +
  scale_color_gradient2(
    low = "#0571b0",
    mid = "#F7F7F7",
    high = "#ca0020",
    midpoint = 0
  ) +
  theme_minimal() +
  labs(
    x = "Coefficient",
    y = "Gene",
    size = "-log10(adj.P.Val)",
    color = "logFC"
  )+
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),  # X-axis label
    axis.title.y = element_text(size = 20, face = "bold")   # Y-axis label
  )


ggsave(
  filename = file.path(outdir, "enrichment_fgsea_Toll_pathway.png"),
  plot = a,
  width = 12,
  height = 10,
  dpi = 300
)

####Plotting the xpression of genes in TOLL pathway-------------------------------

dat.list <- list()
for(gg in genes_toll){
  dat.list[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat <- bind_rows(dat.list, .id = "Gene.name")
dat$treatment <- factor(dat$treatment, levels = c("mock", "Hp", "OMV"))

dat <- dat |>
  mutate(
    sample_treatment_order = factor(
      paste0(donor, ".", treatment),
      levels = c(
        "D1.mock", "D2.mock", "D3.mock",
        "D1.Hp",   "D2.Hp",   "D3.Hp",
        "D1.OMV",  "D2.OMV",  "D3.OMV"
      )
    )
  )


b <- ggplot(dat, aes(x = sample_treatment_order, y = Gene.name, fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#0571b0",
    mid = "#F7F7F7",
    high = "#ca0020",
    midpoint = 0
  ) +
  facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Gene",
    fill = "Scaled expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 20),   # facet labels (mock/Hp/OMV)
    panel.spacing.x = unit(1, "lines")
  )


ggsave(
  filename = file.path(outdir, "Expression_TLR_pathway_genes_batch_corrected.png"),
  plot = b,
  width = 12,
  height = 10,
  dpi = 300
)



####podosomegeneslogfcplots_and_heatmaps-----------------


gene_list <- c(
  "ACTR2","ACTR3","ARPC1B","ARPC2","ARPC3","ARPC4","ARPC5",
  "RAC2","ACTB","ACTG1",
  "ACTN1","ACTN4","FLNA","GSN",
  "CAPG","CAPZA1","CAPZA2","CAPZB","WDR1",
  "ITGB2","ITGAL","ITGB5","TLN1","VCL","PXN","CD44",
  "DCTN1","VAPB","LSP1"
)

# -------------------------
# Filter RNA data
# -------------------------

rna_df <- limmaRes |>
  as.data.frame() |>
  filter(
    Gene.name %in% gene_list,
    coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV")
  )

# -------------------------
# Convert to wide format
# -------------------------

rna_matrix <- rna_df |>
  dplyr::select(Gene.name, coef, logFC) |>
  pivot_wider(
    names_from = coef,
    values_from = logFC
  )

# -------------------------
# RNA logFC comparison plot
# (you can change axes depending on what comparison you want)
# -------------------------

p_rna <- ggplot(
  rna_matrix,
  aes(
    x = treatmentHp,
    y = treatmentOMV
  )
) +
  
  geom_point(color = "darkgreen", size = 3) +
  
  geom_text(
    aes(label = Gene.name),
    vjust = -0.5,
    size = 3
  ) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    title = "RNA-seq Podosome Genes",
    x = "logFC (Treatment Hp)",
    y = "logFC (Treatment OMV)"
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

# -------------------------
# Save plot
# -------------------------



ggsave(
  filename = file.path(outdir, "LogFC_podosome__plot.png"),
  plot = p_rna,
  width = 10,
  height = 10,
  dpi = 300
)
####expression---------


dat.list <- list()
for(gg in gene_list){
  dat.list[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat <- bind_rows(dat.list, .id = "Gene.name")
dat$treatment <- factor(dat$treatment, levels = c("mock", "Hp", "OMV"))

dat <- dat |>
  mutate(
    sample_treatment_order = factor(
      paste0(donor, ".", treatment),
      levels = c(
        "D1.mock", "D2.mock", "D3.mock",
        "D1.Hp",   "D2.Hp",   "D3.Hp",
        "D1.OMV",  "D2.OMV",  "D3.OMV"
      )
    )
  )


vc <- ggplot(dat, aes(x = donor, y = Gene.name, fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#0571b0",
    mid = "#F7F7F7",
    high = "#ca0020",
    )+
  facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Gene",
    fill = "Scaled expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 20),   # facet labels (mock/Hp/OMV)
    panel.spacing.x = unit(1, "lines")
  )

ggsave(
    filename = file.path(outdir, "Expression_genes_regulating_podosomes_batch_corrected.png"),
    plot = vc,
    width = 10,
    height = 10,
    dpi = 300
  )
