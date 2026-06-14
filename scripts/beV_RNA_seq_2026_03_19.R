#load packages

require(tidyverse)
require(limma)
require(patchwork)
require(ComplexHeatmap)
require(pheatmap)
require(dplyr)
library(circlize)
library(edgeR)
library(ggplot2)
library(tibble)
library(stringr)
library(fgsea)
library(msigdbr)
library(forcats)
library(tidyr)
library(enrichR)

#load counts and saving as RDS
counts <- read.csv("D:/PLUS/bEV/raw_data/RNAseq_raw_counts 1.csv", row.names = 1)

# grouping by gene names and take average of genes when there are mutiple - 

df_mean <- counts |>
  group_by(Gene.name) |>
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
counts <- df_mean |> na.omit()
counts <- as.data.frame(counts)
rownames(counts) <- counts$Gene.name
counts$Gene.name <- NULL
#create and load metadata
## Sample column names (replace with colnames(your_data)[1:9])
samples <- colnames(counts)

## Build metadata by splitting on "."
metadata <- data.frame(
  sample    = samples,
  donor     = sub("\\..*$", "", samples),       # text before '.'
  treatment = sub("^.*\\.", "", samples)        # text after '.'
)



#For correlation matrix you need a datamatrix with only number hence removing the first column

corMT <- cor(counts, method = "spearman")
diag(corMT) <- NA   # makes diagonal NAs 
png("D:/PLUS/bEV/results/heatmapcorMT.png",
    width = 1600, height = 1600, res = 300)

pheatmap(corMT,
         scale = "none",
         col = colorRampPalette(c("blue", "white", "red"))(100))

dev.off()



#Dimentionality reduction
#convert the correlation matrix t a 2d distance graph - just to visualize in space if they are segregated. 

mds_plot <- data.frame(cmdscale(dist(2 - corMT), eig = TRUE, k = 2)$points) |>
  add_column(treatment = metadata$treatment) |>
  rownames_to_column("sample") |>
  ggplot(aes(x = X1, y = X2)) + 
  geom_point(aes(color = treatment), size = 3) +
  geom_text(aes(label = metadata$donor), vjust = 1.4) +
  theme_bw()
ggsave("D:/PLUS/bEV/results/mds_plot.png",
       plot = mds_plot,
       width = 8,
       height = 6,
       dpi = 300)



# Relevel treatment so "mock" is baseline

metadata <- metadata |>                      
  mutate(treatment = fct_relevel(treatment, "mock"))

# Create design matrix
mm <- model.matrix(~ donor + treatment, data = metadata)

# Set row names to sample IDs so they appear on y-axis
# Replace 'sample' with the actual column name in metadata containing sample names
rownames(mm) <- metadata$sample

# Save heatmap as PNG
png("D:/PLUS/bEV/results/design_matrix.png",
    width = 1200, height = 800, res = 300)

pheatmap(mm)
dev.off()


#Normalization DGE
#counts matrix ---> create edgeR object ----> remove low expression genes---->normalize library sizes---> ready for DGE

dge <- DGEList(counts = counts)  #edgeR uses tmm normalization
show(dge)
keep <- filterByExpr(dge, design = mm)
dge <- dge[keep, ]
dge <- calcNormFactors(dge)
show(dge)


#Voom is used for transforming counts data to log2 expression values - this is to make the data normal so that linear modeling can be performed. 
dataVoom <- voom(dge, design = mm, plot = TRUE)
show(dataVoom)

#defining a statistical model for DGE
fit <- lmFit(dataVoom, design = mm)
fit <- eBayes(fit) #uses bayesian moderation 
show(coef(fit))

limmaRes <- list() # start an empty list
for(coefx in colnames(coef(fit))){ # run a loop for each coefficient/name coefx chosen by us
  print(coefx)
  # topTable returns the statistics of our genes. We then store the result of each coefficient in a list.
  # The rownames (ENSEMBL Gene IDs) are stored in the column with the name "ensg"
  limmaRes[[coefx]] <- topTable(fit, coef=coefx, number = Inf) |>
    rownames_to_column("Gene.name")
}

limmaRes <- bind_rows(limmaRes, .id = "coef") # bind_rows combines the results and stores the name of the coefficient in the column "coef"
show(limmaRes)
limmaRes <- filter(limmaRes, coef != "(Intercept)") # then we keep all results except for the intercept
show(limmaRes)
colnames(coef(fit))
contrast.mt <- cbind(OMVvsHp = c(0,0,0,-1,1))
row.names(contrast.mt) <- colnames(coef(fit))


# Contrast fit similar to the original limma fit
limmaFit.contrast <- contrasts.fit(fit,contrast.mt)
limmaFit.contrast <- eBayes(limmaFit.contrast)

# Extract results for this contrast coefficient
limmaRes.contrast <- topTable(limmaFit.contrast, coef=colnames(contrast.mt),number = Inf) |>
  rownames_to_column("Gene.name") |>
  mutate(coef=colnames(contrast.mt))

# add them to the full table
limmaRes <- rbind(limmaRes.contrast, limmaRes) # add this coefficient to the result table
table(limmaRes$coef)


# data visulaization ------------------------------------------------------

# ---- Volcano plot (adj.P.Val) ----
volcano_adj <- limmaRes |>
  ggplot(aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(alpha = 0.3) +
  facet_grid(cols = vars(coef))

ggsave("D:/PLUS/bEV/results/volcano_adjPval.png",
       plot = volcano_adj, width = 10, height = 6, dpi = 300)


# ---- Volcano hex plot ----
volcano_hex <- limmaRes |>
  ggplot(aes(x = logFC, y = -log10(P.Value))) + 
  geom_hex() +
  facet_grid(cols = vars(coef))

ggsave("D:/PLUS/bEV/results/volcano_hex.png",
       plot = volcano_hex, width = 10, height = 6, dpi = 300)


# ---- P-value distribution ----
pval_hist <- limmaRes |>
  ggplot(aes(x = P.Value, fill = factor(floor(AveExpr)))) +
  geom_histogram() +
  facet_grid(cols = vars(coef))

ggsave("D:/PLUS/bEV/results/pvalue_distribution.png",
       plot = pval_hist, width = 10, height = 6, dpi = 300)


# # ---- Filter lowly expressed genes ----
# limmaRes_f <- limmaRes |> 
#   filter(AveExpr > -5)
# 
# 
# # ---- P-value distribution after filtering ----
# pval_hist_filtered <- limmaRes_f |>
#   ggplot(aes(x = P.Value, fill = factor(floor(AveExpr)))) + 
#   geom_histogram() +
#   facet_grid(cols = vars(coef))
# 
# ggsave("D:/PLUS/bEV/results/pvalue_distribution_filtered.png",
#        plot = pval_hist_filtered, width = 10, height = 6, dpi = 300)
# 
# 
# # ---- Volcano plot after filtering ----
# volcano_filtered <- limmaRes_f |>
#   ggplot(aes(x = logFC, y = -log10(adj.P.Val))) + 
#   geom_point(alpha = 0.3) +
#   facet_grid(cols = vars(coef))
# 
# ggsave("D:/PLUS/bEV/results/volcano_filtered.png",
#        plot = volcano_filtered, width = 10, height = 6, dpi = 300)
# 
# 
# # ---- Continue analysis with filtered results ----
# limmaRes <- limmaRes_f


# ---- Extract significant genes ----
limmaResSig <- limmaRes |> 
  filter(adj.P.Val < 0.05)

count(limmaResSig, coef)

# -------------------------
# Load your counts CSV (gmap)
#gmap <- read.csv("/Users/aarathyrg/Documents/PLUS/bEV/raw_data/RNAseq_raw_counts 1.csv")

#gmap <- gmap |>
 # rename(ensg = ID)

# Inspect columns
#head(gmap)
# Ensure there is a column for ENSEMBL IDs (say "ensg") and "gene.name"
#colnames(gmap)
#gmap <- gmap[, c("ensg","Gene.name")]
#head(gmap)
#limmaRes$gene <- gmap$Gene.name[match(limmaRes$ensg, gmap$ensg)]
#limmaResSig$gene <- gmap$Gene.name[match(limmaResSig$ensg, gmap$ensg)]
# -------------------------



#-------------OMVvsHp and OMVvsmock enriched and absent in HP vs mock-----------------

################removing batch effect--------------

design2 <- model.matrix(~ treatment, data = metadata)
expr_corrected <- removeBatchEffect(dataVoom$E, batch = metadata$donor, 
                                    design = design2)



# OMV vs mock (significant)
deg_OMV_mock_up <- limmaRes |>
  filter(coef == "treatmentOMV",
         adj.P.Val < 0.05,
         logFC > 1)
# OMV vs mock (significant)
deg_OMV_mock_down <- limmaRes |>
  filter(coef == "treatmentOMV",
         adj.P.Val < 0.05,
         logFC < -1)

# OMV vs Hp (significant)
deg_OMV_Hp <- limmaRes |>
  filter(coef == "OMVvsHp",
         adj.P.Val < 0.05,
         abs(logFC) > 1)

# Hp vs mock (NOT significant)
deg_Hp_mock_nonsig <- limmaRes |>
  filter(coef == "treatmentHp",
         adj.P.Val > 0.05)

# Combine logic
genes_OMV_specific_up <- deg_OMV_mock_up |>
  inner_join(deg_OMV_Hp, by = "Gene.name") |>
  inner_join(deg_Hp_mock_nonsig, by = "Gene.name") |>
  pull(Gene.name) |>
  unique()

genes_OMV_specific_down <- deg_OMV_mock_down |>
  inner_join(deg_OMV_Hp, by = "Gene.name") |>
  inner_join(deg_Hp_mock_nonsig, by = "Gene.name") |>
  pull(Gene.name) |>
  unique()


#upregulated genes specific to OMV---------------
dat.list_OMV <- list()
for(gg in genes_OMV_specific_up){
  dat.list_OMV[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat_OMV <- bind_rows(dat.list_OMV, .id = "Gene.name")

dat_OMV <- dat_OMV |>
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


d <- ggplot(dat_OMV, aes(x = donor, y = Gene.name, fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Gene",
    fill = "Scaled expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 20),   # bigger x label
    axis.title.y = element_text(size = 20),   # bigger y label
    plot.title = element_text(size = 24, face = "bold"),  # big heading
    strip.text = element_text(size = 16),  # facet labels
    panel.spacing.x = unit(1, "lines")
  )

ggsave("D:/PLUS/bEV/results/Expression_of_OMV_specific_upregulated_DEGs_batchcorrected.png",
       plot = d, width = 12, height = 15, dpi = 300)



###---downregulated genes specific to OMV-------------
dat.list_OMV <- list()
for(gg in genes_OMV_specific_down){
  dat.list_OMV[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat_OMV <- bind_rows(dat.list_OMV, .id = "Gene.name")

dat_OMV <- dat_OMV |>
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


d <- ggplot(dat_OMV, aes(x = donor, y = Gene.name, fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Gene",
    fill = "Scaled expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 20),   # bigger x label
    axis.title.y = element_text(size = 20),   # bigger y label
    plot.title = element_text(size = 24, face = "bold"),  # big heading
    strip.text = element_text(size = 16),  # facet labels
    panel.spacing.x = unit(1, "lines")
  )

ggsave("D:/PLUS/bEV/results/Expression_of_OMV_specific_downregulated_DEGs_batchcorrected.png",
       plot = d, width = 12, height = 15, dpi = 300)






#Enrichment analysis using fgsea-----------------------

output_dir <- "D:/PLUS/bEV/results/fgseaplots"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}



# Load MSigdbr
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
    low = "blue", mid = "white", high = "red", midpoint = 0
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

ggsave("D:/PLUS/bEV/results/fgseaplots/enrichment_fgsea_combined_new.png",
       plot = ggfgsea_combined, width = 12, height = 10, dpi = 300)


# #-------Filter genes in 3 interesting pathways----------
# 
# library(purrr)
# 
# genes <- topPathways |>
#   filter(pathway == "HALLMARK_MITOTIC_SPINDLE") |>
#   pull(leadingEdge) |>
#   str_remove_all('c\\(|\\)|"|\\n') |>
#   str_split(",\\s*") |>
#   unlist()
# head(genes)
# 
# goi_pathway <- limmaRes |>
#   filter(Gene.name%in%genes)
# 
# 
# a <- ggplot(goi_pathway, aes(x = coef, y = Gene.name)) +
#   geom_point(aes(size = -log10(adj.P.Val), color = logFC)) +
#   scale_color_gradient2(
#     low = "blue",      # downregulated
#     mid = "white",     # ~0
#     high = "red",      # upregulated
#     midpoint = 0
#   )+
#   theme_minimal() +
#   labs(
#     x = "Coefficient",
#     y = "Gene",
#     size = "-log10(adj.P.Val)",
#     color = "logFC"
#   )
#   
#   ggsave("D:/PLUS/bEV/results/enrichment_fgsea_Mitotic_spindle_Pathway.png",
#          plot = a, width = 12, height = 15, dpi = 300)
#   
#   
#  # -----------plotting expression data for these genes----
#   
#   dat.list <- list()
#   for(gg in genes){
#     dat.list[[gg]] <- metadata |>
#       mutate(Exp=scale(dataVoom$E[gg,])) |>
#       remove_rownames()
#   }
#   dat <- bind_rows(dat.list, .id = "Gene.name")
#   
#   dat <- dat |>
#     mutate(
#       sample_treatment_order = factor(
#         paste0(donor, ".", treatment),
#         levels = c(
#           "D1.mock", "D2.mock", "D3.mock",
#           "D1.Hp",   "D2.Hp",   "D3.Hp",
#           "D1.OMV",  "D2.OMV",  "D3.OMV"
#         )
#       )
#     )
#   
#   
#   b <- ggplot(dat, aes(x = donor, y = Gene.name, fill = Exp)) +
#     geom_tile() +
#     scale_fill_gradient2(
#       low = "blue", mid = "white", high = "red", midpoint = 0
#     ) +
#     facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
#     theme_minimal() +
#     labs(
#       x = "Sample",
#       y = "Gene",
#       fill = "Scaled expression"
#     ) +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       panel.spacing.x = unit(1, "lines")  # adds space between mock / Hp / OMV
#     )
#   
#   ggsave("D:/PLUS/bEV/results/Expression_Mitogic_pathway_genes.png",
#          plot = b, width = 12, height = 15, dpi = 300)
  
#####-------Filter genes in 3 interesting pathways----------
  
  library(purrr)
  
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
      low = "blue",      # downregulated
      mid = "white",     # ~0
      high = "red",      # upregulated
      midpoint = 0
    )+
    theme_minimal() +
    labs(
      x = "Coefficient",
      y = "Gene",
      size = "-log10(adj.P.Val)",
      color = "logFC"
    )
  
  ggsave("D:/PLUS/bEV/results/fgseaplots/enrichment_fgsea_Epithilial_Mesenchymal_Pathway.png",
         plot = a, width = 12, height = 15, dpi = 300)
  
  
  # -----------plotting expression data for these genes----
  
  dat.list <- list()
  for(gg in genes){
    dat.list[[gg]] <- metadata |>
      mutate(Exp=scale(expr_corrected[gg,])) |>
      remove_rownames()
  }
  dat <- bind_rows(dat.list, .id = "Gene.name")
  
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
      low = "blue", mid = "white", high = "red", midpoint = 0
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
  
  ggsave("D:/PLUS/bEV/results/fgseaplots/Expression_Epithilial_Mesenchymal_Pathway_genes_batch_corrected.png",
         plot = b, width = 12, height = 15, dpi = 300)
  
  #
 ####----------------------
  
  
  # dat.list <- list()
  # for(gg in genes){
  #   dat.list[[gg]] <- metadata |>
  #     mutate(Exp=scale(expr_corrected[gg,])) |>
  #     remove_rownames()
  # }
  # dat <- bind_rows(dat.list, .id = "Gene.name")
  # 
  # dat <- dat |>
  #   mutate(
  #     sample_treatment_order = factor(
  #       paste0(donor, ".", treatment),
  #       levels = c(
  #         "D1.mock", "D2.mock", "D3.mock",
  #         "D1.Hp",   "D2.Hp",   "D3.Hp",
  #         "D1.OMV",  "D2.OMV",  "D3.OMV"
  #       )
  #     )
  #   )
  # 
  # 
  # b <- ggplot(dat, aes(x = donor, y = Gene.name, fill = Exp)) +
  #   geom_tile() +
  #   scale_fill_gradient2(
  #     low = "blue", mid = "white", high = "red", midpoint = 0
  #   ) +
  #   facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
  #   theme_minimal() +
  #   labs(
  #     x = "Sample",
  #     y = "Gene",
  #     fill = "Scaled expression"
  #   ) +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1),
  #     panel.spacing.x = unit(1, "lines")  # adds space between mock / Hp / OMV
  #   )
  # 
  # ggsave("D:/PLUS/bEV/results/Expression_Mitogic_pathway_genes_batch_corrected.png",
  #        plot = b, width = 12, height = 15, dpi = 300)
  # 
  # 
######actin/cytoskeleton/Intergrin/TLR---------  

cyto_terms <- fgsea_all_combined |>
    filter(str_detect(pathway, regex("actin|cytoskeleton", ignore_case = TRUE)))
  
show(cyto_terms)

integrin_terms <- fgsea_all_combined |>
  filter(str_detect(pathway, regex("integrin|GPCRs", ignore_case = TRUE)))

show(integrin_terms)

TLR_terms <- fgsea_all_combined |>
  filter(str_detect(pathway, regex("TLR|Toll", ignore_case = TRUE)))

show(TLR_terms)


######_____Filter out top 100 DEGs from LimmaRes and plotting their expression-----
  

# top_100 <- limmaRes |>
#   filter(coef == "OMVvsHp")|>
#   filter(adj.P.Val<0.05) |>
#   arrange(desc(abs(logFC))) |>
#   slice_head(n=100)|>
#   filter(abs(logFC)>1)|>
#   pull(Gene.name)
# 
# dat.list_100 <- list()
# for(gg in top_100){
#   dat.list_100[[gg]] <- metadata |>
#     mutate(Exp=scale(expr_corrected[gg,])) |>
#     remove_rownames()
# }
# dat_100 <- bind_rows(dat.list_100, .id = "Gene.name")
# 
# dat_100 <- dat_100 |>
#   mutate(
#     sample_treatment_order = factor(
#       paste0(donor, ".", treatment),
#       levels = c(
#         "D1.mock", "D2.mock", "D3.mock",
#         "D1.Hp",   "D2.Hp",   "D3.Hp",
#         "D1.OMV",  "D2.OMV",  "D3.OMV"
#       )
#     )
#   )
# 
# 
# c <- ggplot(dat_100, aes(x = sample_treatment_order, y = Gene.name, fill = Exp)) +
#   geom_tile() +
#   scale_fill_gradient2(
#     low = "blue", mid = "white", high = "red", midpoint = 0
#   ) +
#   theme_minimal() +
#   labs(
#     x = "Sample",
#     y = "Gene",
#     fill = "Scaled expression"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# ggsave("D:/PLUS/bEV/results/Expression of top 100 DEGs.png",
#        plot = c, width = 12, height = 15, dpi = 300)
# 



# ###-------------- Upregulated and downregulated---------- 
# 
# dat.list_OMV <- list()
# for(gg in genes_OMV_specific_up){
#   dat.list_OMV[[gg]] <- metadata |>
#     mutate(Exp=scale(dataVoom$E[gg,])) |>
#     remove_rownames()
# }
# dat_OMV <- bind_rows(dat.list_OMV, .id = "Gene.name")
# 
# dat_OMV <- dat_OMV |>
#   mutate(
#     sample_treatment_order = factor(
#       paste0(donor, ".", treatment),
#       levels = c(
#         "D1.mock", "D2.mock", "D3.mock",
#         "D1.Hp",   "D2.Hp",   "D3.Hp",
#         "D1.OMV",  "D2.OMV",  "D3.OMV"
#       )
#     )
#   )
# 
# 
# e <- ggplot(dat_OMV, aes(x = donor, y = Gene.name, fill = Exp)) +
#   geom_tile() +
#   scale_fill_gradient2(
#     low = "blue", mid = "white", high = "red", midpoint = 0
#   ) +
#   facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
#   theme_minimal() +
#   labs(
#     x = "Sample",
#     y = "Gene",
#     fill = "Scaled expression"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.spacing.x = unit(1, "lines")  # adds space between mock / Hp / OMV
#   )
# 
# ggsave("D:/PLUS/bEV/results/Expression of OMV specific DEGs_upregulated.png",
#        plot = e, width = 12, height = 15, dpi = 300)
# 
# 
# ####OMVspecific degs down batchcorrected_________________
# dat.list_OMV <- list()
# for(gg in genes_OMV_specific_up){
#   dat.list_OMV[[gg]] <- metadata |>
#     mutate(Exp=scale(dataVoom$E[gg,])) |>
#     remove_rownames()
# }
# dat_OMV <- bind_rows(dat.list_OMV, .id = "Gene.name")
# 
# dat_OMV <- dat_OMV |>
#   mutate(
#     sample_treatment_order = factor(
#       paste0(donor, ".", treatment),
#       levels = c(
#         "D1.mock", "D2.mock", "D3.mock",
#         "D1.Hp",   "D2.Hp",   "D3.Hp",
#         "D1.OMV",  "D2.OMV",  "D3.OMV"
#       )
#     )
#   )
# 
# 
# e <- ggplot(dat_OMV, aes(x = donor, y = Gene.name, fill = Exp)) +
#   geom_tile() +
#   scale_fill_gradient2(
#     low = "blue", mid = "white", high = "red", midpoint = 0
#   ) +
#   facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
#   theme_minimal() +
#   labs(
#     x = "Sample",
#     y = "Gene",
#     fill = "Scaled expression"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.spacing.x = unit(1, "lines")  # adds space between mock / Hp / OMV
#   )
# 
# ggsave("D:/PLUS/bEV/results/Expression of OMV specific DEGs_downregulated.png",
#        plot = e, width = 12, height = 15, dpi = 300)
# 
# 
# 
# 



#------------------------------------------ enrichR on upregulated genes------, 
output_dir <- "D:/PLUS/bEV/results/enrichmentplots"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


ENRICHR.DBS <- c("GO_Biological_Process_2021",
                 "TRRUST_Transcription_Factors_2019",
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
  scale_colour_gradient2(high = "red", low = "blue")+
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

ggsave("D:/PLUS/bEV/results/enrichmentplots/Pathways_enrichment_upregulatedgenes_databases.png",
       plot = p, width = 10, height = 15, dpi = 300)

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
  scale_colour_gradient2(high = "red", low = "blue")+
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
ggsave("D:/PLUS/bEV/results/enrichmentplots/Pathways_enrichment_downregulatedgenes_databases.png",
       plot = f, width = 10, height = 15, dpi = 300)



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
    low = "blue",      # downregulated
    mid = "white",     # ~0
    high = "red",      # upregulated
    midpoint = 0
  )+
  theme_minimal() +
  labs(
    x = "Coefficient",
    y = "Gene",
    size = "-log10(adj.P.Val)",
    color = "logFC"
  )

ggsave("D:/PLUS/bEV/results/enrichmentplots/Expression_of_genes_cellular_response_to_chemokine.png",
       plot = i, width = 12, height = 15, dpi = 300)

#------------------------expression list of the genes selected batch corrected--------

dat.list_chemo <- list()
for(gg in gene_list){
  dat.list_chemo[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat_chemo <- bind_rows(dat.list_chemo, .id = "Gene.name")

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
    low = "blue", mid = "white", high = "red", midpoint = 0
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

ggsave("D:/PLUS/bEV/results/enrichmentplots/Expression of genes regulating cellular response to chemokines batch corrected.png",
       plot = k, width = 12, height = 15, dpi = 300)





# #------------------------expression list of the genes selected--------
# 
# dat.list_chemo <- list()
# for(gg in gene_list){
#   dat.list_chemo[[gg]] <- metadata |>
#     mutate(Exp=scale(dataVoom$E[gg,])) |>
#     remove_rownames()
# }
# dat_chemo <- bind_rows(dat.list_chemo, .id = "Gene.name")
# 
# dat_chemo <- dat_chemo |>
#   mutate(
#     sample_treatment_order = factor(
#       paste0(donor, ".", treatment),
#       levels = c(
#         "D1.mock", "D2.mock", "D3.mock",
#         "D1.Hp",   "D2.Hp",   "D3.Hp",
#         "D1.OMV",  "D2.OMV",  "D3.OMV"
#       )
#     )
#   )
# 
# 
# j <- ggplot(dat_chemo, aes(x = donor, y = Gene.name, fill = Exp)) +
#   geom_tile() +
#   scale_fill_gradient2(
#     low = "blue", mid = "white", high = "red", midpoint = 0
#   ) +
#   facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
#   theme_minimal() +
#   labs(
#     x = "Sample",
#     y = "Gene",
#     fill = "Scaled expression"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.spacing.x = unit(1, "lines")  # adds space between mock / Hp / OMV
#   )
# 
# ggsave("D:/PLUS/bEV/results/Expression of genes regulating cellular response to chemokine.png",
#        plot = j, width = 12, height = 15, dpi = 300)
# 
# 



####################------------ upregulated vs downregulated across coefficients-----------


# #Downregulated-----------
# down_genes <- limmaRes |>
#   filter(adj.P.Val < 0.05, logFC < -1) |>
#   count(coef, name = "downregulated")
# 
# # Upregulated
# up_genes <- limmaRes |>
#   filter(adj.P.Val < 0.05, logFC > 1) |>
#   count(coef, name = "upregulated")
# 
# gene_counts <- full_join(up_genes, down_genes, by = "coef")
# 
# l <- ggplot(gene_counts) +
#   geom_col(aes(x = coef, y = upregulated, fill = "Upregulated"),
#            width = 0.3, position = position_nudge(x = -0.15)) +
#   geom_col(aes(x = coef, y = downregulated, fill = "Downregulated"),
#            width = 0.3, position = position_nudge(x = 0.15)) +
#   geom_hline(yintercept = 0, color = "black") +
#   scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
#   theme_minimal(base_size = 16) +  # increases text sizes
#   labs(
#     x = "Coefficient",
#     y = "Number of significant genes",
#     fill = "Direction",
#     title = "Upregulated (red) and Downregulated (blue) genes per coefficient"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
#     axis.text.y = element_text(size = 14),
#     axis.title.x = element_text(size = 16),
#     axis.title.y = element_text(size = 16),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 12),
#     plot.title = element_text(size = 18, face = "bold")
#   )
# 
# # Save the plot
# ggsave("D:/PLUS/bEV/results/Upregulated_vs_Downregulated_across_coefs.png",
#        plot = l, width = 12, height = 8, dpi = 300)

##-----50 up and downregulated in OMV vs HP------


top50_up_genes <- limmaRes |>
  filter(coef == "OMVvsHp", adj.P.Val < 0.05, logFC > 0) |>
  arrange(desc(logFC)) |>
  slice_head(n = 50) |>
  pull(Gene.name)

top50_down_genes <- limmaRes |>
  filter(coef == "OMVvsHp",adj.P.Val < 0.05, logFC < 0) |>
  arrange(logFC) |>
  slice_head(n = 50) |>
  pull(Gene.name)

gene_list_OMVvsHp <- c(top50_up_genes, top50_down_genes)

dat.list_OMVvsHp <- list()
for(gg in top50_up_genes){
  dat.list_OMVvsHp[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat_OMVvsHp <- bind_rows(dat.list_OMVvsHp, .id = "Gene.name")

dat_OMVvsHp <- dat_OMVvsHp |>
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


q <- ggplot(dat_OMVvsHp, aes(x = donor, y = Gene.name, fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Gene",
    fill = "Scaled expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 14),   # facet labels (mock/Hp/OMV)
    panel.spacing.x = unit(1, "lines")
  )

ggsave("D:/PLUS/bEV/results/Expression_of_upregulated_genes_significantly_changed_in_OMVvsHp_batchcorrected.png",
       plot = q, width = 12, height = 15, dpi = 300)



##downregulated-----

dat.list_OMVvsHp <- list()
for(gg in top50_down_genes){
  dat.list_OMVvsHp[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat_OMVvsHp <- bind_rows(dat.list_OMVvsHp, .id = "Gene.name")

dat_OMVvsHp <- dat_OMVvsHp |>
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


r <- ggplot(dat_OMVvsHp, aes(x = donor, y = Gene.name, fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Gene",
    fill = "Scaled expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 14),   # facet labels (mock/Hp/OMV)
    panel.spacing.x = unit(1, "lines")
  )

ggsave("D:/PLUS/bEV/results/Expression_of_downregulated_genes_significantly_changed_in_OMVvsHp_batchcorrected.png",
       plot = r, width = 12, height = 15, dpi = 300)


#######enrichR_on_OMV_specific_gene_set----------------


coef_list <- unique(limmaRes$coef)

enrich_upregulated <- list()

for (coef_current in coef_list) {
  
  tt <- limmaRes |>
    filter(coef == coef_current)
  
  results <- enrichr(genes_OMV_specific_down, ENRICHR.DBS)
  
  combined <- do.call(rbind, lapply(names(results), function(db) {
    df <- results[[db]]
    df$Database <- db
    df$Coefficient <- coef_current
    return(df)
  }))
  
  enrich_upregulated[[coef_current]] <- combined
}

final_combined_upregulated_OMV_up <- bind_rows(enrich_upregulated)



#select the arrange based on significance 
combined_sig_OMV_up <- final_combined_upregulated |>
  filter(Adjusted.P.value < 0.05) |>
  arrange(Adjusted.P.value)

#selecting top 5 enriched pathways per data base
top_per_db_OMV_up <- combined_sig_OMV_up |>
  filter(Coefficient %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  group_by(Database, Coefficient) |>
  slice_head(n = 5) |>
  ungroup() 

#plotting for each database
p <- ggplot(top_per_db_OMV_up,
            aes(x = Coefficient,                   # numeric x-axis
                y = reorder(Term, Coefficient),   # pathways on y-axis, ordered
                size = -log10(Adjusted.P.value),
                color = log2(Odds.Ratio))) +
  geom_point(alpha = 0.8) +                 # semi-transparent points
  scale_size_continuous(range = c(3,6)) + 
  scale_colour_gradient2(high = "red", low = "blue")+
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

ggsave("D:/PLUS/bEV/results/Pathways_enrichment_OMV_specific_downregulated_databases.png",
       plot = p, width = 10, height = 15, dpi = 300)






# 
# #######list_pathways_with_selective terms_
# 
# 
# #####enrichR_on_terms_specific_for pathways----
# 
# term_enrich <- final_combined_upregulated |>
#   filter(grepl("integrin|TLR|toll|actin|cytoskeleton", Term, ignore.case = TRUE)) |>
#   dplyr::select(Database, Coefficient, Term, Adjusted.P.value, Combined.Score, Genes)
# show(term_enrich)
# 
# #select the arrange based on significance 
# combined_sig_term <- term_enrich |>
#   filter(Adjusted.P.value < 0.05) |>
#   arrange(Adjusted.P.value)
# 
# #selecting top 5 enriched pathways per data base
# top_per_db_term <- combined_sig_term |>
#   group_by(Database, Coefficient) |>
#   slice_head(n = 5) |>
#   ungroup() 
# 
# #plotting for each database
# s <- ggplot(top_per_db_term,
#             aes(x = Coefficient,                   # numeric x-axis
#                 y = reorder(Term, Coefficient),   # pathways on y-axis, ordered
#                 size = -log10(Adjusted.P.value),
#                 color = log2(Odds.Ratio))) +
#   geom_point(alpha = 0.8) +                 # semi-transparent points
#   scale_size_continuous(range = c(3,6)) +
#   scale_colour_gradient2(high = "red", low = "blue")+
#   labs(
#     title = "Top Pathways per Database",
#     x = "Coefficient",
#     y = "Pathway / Term",
#     size = "-log10(Adjusted P Value)",
#     color = "log2(Odds Ratio)"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_text(size = 12),   # pathway names bigger
#     axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
#     axis.title.x = element_text(size = 18),
#     axis.title.y = element_text(size = 18),
#     plot.title = element_text(size = 22, face = "bold"),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 12),
#     strip.text = element_text(size = 14),
#     panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
#     panel.grid.minor = element_line(color = "grey90", linewidth = 0.25)
#   )
# ggsave("D:/PLUS/bEV/results/Pathways_enrichment_upregulatedgenes_actin_cytoskeletal_TLR_Integrin_regulation.png",
#        plot = s, width = 20, height = 15, dpi = 300)





#####fgsea_with_additional_pathways_----


library(msigdbr)
library(dplyr)

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
    low = "blue", mid = "white", high = "red", midpoint = 0
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

ggsave("D:/PLUS/bEV/results/enrichment_fgsea_combined_actin_TLR_Integrin.png",
       plot = ggfgsea_combined, width = 12, height = 10, dpi = 300)
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
    low = "blue",      # downregulated
    mid = "white",     # ~0
    high = "red",      # upregulated
    midpoint = 0
  )+
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

ggsave("D:/PLUS/bEV/results/enrichment_fgsea_Toll_pathway.png",
       plot = a, width = 12, height = 15, dpi = 300)



####Plotting the xpression of genes in TOLL pathway-------------------------------

dat.list <- list()
for(gg in genes_toll){
  dat.list[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat <- bind_rows(dat.list, .id = "Gene.name")

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
    low = "blue", mid = "white", high = "red", midpoint = 0
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


ggsave("D:/PLUS/bEV/results/Expression_TLR_pthway_genes_batch_corrected.png",
       plot = b, width = 12, height = 15, dpi = 300)

#########Filtering TLR-2/4 upstream and downstream genes----------

# # Upstream: receptors and adaptors
# tlr2_4_upstream <- c(
#   "TLR2", "TLR4",       # receptors
#   "CD14", "LY96",       # co-receptors / LPS recognition
#   "MYD88", "TIRAP",     # adaptor proteins
#   "IRAK1", "IRAK4",     # early kinases
#   "TRAF6"               # E3 ubiquitin ligase
# )
# 
# # Downstream: transcription factors, cytokines, chemokines, effectors
# tlr2_4_downstream <- c(
#   "MAP3K7", "MAP3K8",        # kinases
#   "TAB1","TAB2","TAB3",      # scaffolds
#   "CHUK", "IKBKB", "IKBKG",  # IKK complex
#   "RIPK1", "NFKBIA",         # NF-κB regulators
#   "TNF", "IL1B", "IL6", "IFNB1",       # cytokines
#   "CXCL8","CXCL10","CXCL11","CXCL9",  # chemokines
#   "CCL3","CCL4","CCL5",                 # chemokines
#   "STAT1","IRF7","JUN","FOS",          # transcription factors
#   "CD40","CD80"                         # co-stimulatory molecules
# )
# 
# 
# tlr2_4_genes <- data.frame(
#   Gene.name = unique(c(tlr2_4_upstream, tlr2_4_downstream))
# ) |>
#   mutate(layer = case_when(
#     Gene.name %in% tlr2_4_upstream ~ "Upstream",
#     Gene.name %in% tlr2_4_downstream ~ "Downstream",
#     TRUE ~ "Other"
#   ))
# 
# goi_toll_select <- limmaRes |>
#   filter(coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
#   filter(Gene.name%in%tlr2_4_genes$Gene.name)
# 
# show(goi_toll_select)
# 
# goi_toll_select <- goi_toll_select |>
#   left_join(tlr2_4_genes, by = "Gene.name")
# 
# # Quick check
# head(goi_toll_select)
# table(goi_toll_select$layer)  # counts of Upstream vs Downstream
# 
# q <- ggplot(goi_toll_select, aes(x = coef, y = reorder(Gene.name, logFC))) +
#   geom_point(aes(size = -log10(adj.P.Val + 1e-10), color = logFC)) +
#   facet_wrap(~ layer, scales = "free_y") +  # separates Upstream vs Downstream
#   scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
#   theme_minimal() +
#   labs(
#     x = "Coefficient",
#     y = "Gene",
#     size = "-log10(adj.P.Val)",
#     color = "logFC"
#   )+
#   theme(
#     axis.title.x = element_text(size = 20, face = "bold"),  # X-axis label
#     axis.title.y = element_text(size = 20, face = "bold")   # Y-axis label
#   )
# 
# ggsave("D:/PLUS/bEV/results/Expression_TLR_pthway_genes_upstream_vs_downstream.png",
#        plot = q, width = 12, height = 15, dpi = 300)
# 
# ####expression_of_these_sets_of_genes------------
# 
# dat.list <- list()
# for(gg in tlr2_4_genes$Gene.name){
#   dat.list[[gg]] <- metadata |>
#     mutate(Exp=scale(expr_corrected[gg,])) |>
#     remove_rownames()
# }
# dat <- bind_rows(dat.list, .id = "Gene.name")
# 
# dat <- dat |> 
#   left_join(tlr2_4_genes, by = "Gene.name")  # adds the 'layer' column
# 
# dat <- dat |>
#   mutate(
#     sample_treatment_order = factor(
#       paste0(donor, ".", treatment),
#       levels = c(
#         "D1.mock", "D2.mock", "D3.mock",
#         "D1.Hp",   "D2.Hp",   "D3.Hp",
#         "D1.OMV",  "D2.OMV",  "D3.OMV"
#       )
#     )
#   )
# 
# 
# t <- ggplot(dat, aes(x = sample_treatment_order, y = reorder(Gene.name, Exp), fill = Exp)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
#   facet_wrap(~ layer, scales = "free_y") +  # Upstream vs Downstream
#   theme_minimal() +
#   labs(
#     x = "Sample",
#     y = "Gene",
#     fill = "Scaled expression"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
#     axis.text.y = element_text(size = 18),
#     axis.title.x = element_text(size = 20),
#     axis.title.y = element_text(size = 20),
#     strip.text = element_text(size = 20),   # facet labels (mock/Hp/OMV)
#     panel.spacing.x = unit(1, "lines")
#   )
# 
# ggsave("D:/PLUS/bEV/results/Expression_TLR2_4_pathway_selected_genes_batch_corrected.png",
#        plot = t, width = 12, height = 15, dpi = 300)

#######Filtering Intergin upstream and downstream genes------

# Upstream: receptors and adaptors
tlr2_4_upstream <- c(
  "TLR2", "TLR4",       # receptors
  "CD14", "LY96",       # co-receptors / LPS recognition
  "MYD88", "TIRAP",     # adaptor proteins
  "IRAK1", "IRAK4",     # early kinases
  "TRAF6"               # E3 ubiquitin ligase
)

# Downstream: transcription factors, cytokines, chemokines, effectors
tlr2_4_downstream <- c(
  "MAP3K7", "MAP3K8",        # kinases
  "TAB1","TAB2","TAB3",      # scaffolds
  "CHUK", "IKBKB", "IKBKG",  # IKK complex
  "RIPK1", "NFKBIA",         # NF-κB regulators
  "TNF", "IL1B", "IL6", "IFNB1",       # cytokines
  "CXCL8","CXCL10","CXCL11","CXCL9",  # chemokines
  "CCL3","CCL4","CCL5",                 # chemokines
  "STAT1","IRF7","JUN","FOS",          # transcription factors
  "CD40","CD80"                         # co-stimulatory molecules
)

###Crosstalk regulators-----------

tlr_regulators <- c(
  "WNT5A",
  "SRC",
  "SYK",
  "CDC42",
  "RAC1",
  "RHOA",
  "PIK3CA",
  "PIK3R1",
  "AKT1",
  "AKT2"
)


tlr2_4_genes <- data.frame(
  Gene.name = unique(c(tlr2_4_upstream, tlr2_4_downstream, tlr_regulators))
) |>
  mutate(layer = case_when(
    Gene.name %in% tlr2_4_upstream ~ "Upstream",
    Gene.name %in% tlr2_4_downstream ~ "Downstream",
    Gene.name %in% tlr_regulators ~ "Crosstalk_Regulators",
    TRUE ~ "Other"
  ))

goi_toll_select <- limmaRes |>
  filter(coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  filter(Gene.name%in%tlr2_4_genes$Gene.name)

show(goi_toll_select)

goi_toll_select <- goi_toll_select |>
  left_join(tlr2_4_genes, by = "Gene.name")

# Quick check
head(goi_toll_select)
table(goi_toll_select$layer)  # counts of Upstream vs Downstream

q <- ggplot(goi_toll_select, aes(x = coef, y = reorder(Gene.name, logFC))) +
  geom_point(aes(size = -log10(adj.P.Val), color = logFC)) +
  facet_wrap(~ layer, scales = "free_y") +  # separates Upstream vs Downstream
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    x = "Coefficient",
    y = "Gene",
    size = "-log10(adj.P.Val)",
    color = "logFC"
  )+
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),  # X-axis label
    axis.title.y = element_text(size = 30, face = "bold"),  # Y-axis label
    axis.text.y = element_text(size = 20),                 # gene names
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  # tilt X labels
    strip.text = element_text(size = 18)    # facet headings
  )

ggsave("D:/PLUS/bEV/results/Expression_TLR_pthway_genes_upstream_vs_downstream.png",
       plot = q, width = 15, height = 15, dpi = 300)

####expression_of_these_sets_of_genes------------

dat.list <- list()
for(gg in tlr2_4_genes$Gene.name){
  dat.list[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat <- bind_rows(dat.list, .id = "Gene.name")

dat <- dat |> 
  left_join(tlr2_4_genes, by = "Gene.name")  # adds the 'layer' column

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


t <- ggplot(dat, aes(x = sample_treatment_order, y = reorder(Gene.name, Exp), fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_wrap(~ layer, scales = "free_y") +  # Upstream vs Downstream
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Gene",
    fill = "Scaled expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 16),   # facet labels (mock/Hp/OMV)
    panel.spacing.x = unit(1, "lines")
  )

ggsave("D:/PLUS/bEV/results/Expression_TLR2_4_pathway_selected_genes_batch_corrected.png",
       plot = t, width = 12, height = 15, dpi = 300)

####WNT5a_TLR_papaer_genes

wnt5a_tlr_genes <- c(
  "WNT5A","WNT3A",
  "TLR2","TLR4","LY96",
  "MYD88",
  "NFKB1","RELA","NFKBIA",
  "MAPK1","MAPK3","MAPK8","MAPK14",
  "JUN","FOS",
  "IL10","TNF","IL6","CXCL8","CCL2",
  "S100A9",
  "PIK3CA","PIK3R1",
  "ROR2","RYK",
  "FZD2","FZD5","FZD7",
  "DVL1","DVL2","DVL3",
  "CAMK2A","CAMK2B","CAMK2D",
  "RHOA","RAC1"
)


goi_wnt5a_tlr_select <- limmaRes |>
  filter(coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  filter(Gene.name%in%wnt5a_tlr_genes)

show(goi_wnt5a_tlr_select)

# Quick check
head(goi_wnt5a_tlr_select)


q <- ggplot(goi_wnt5a_tlr_select, aes(x = coef, y = reorder(Gene.name, logFC))) +
  geom_point(aes(size = -log10(adj.P.Val), color = logFC)) +
  #facet_wrap(~ layer, scales = "free_y") +  # separates Upstream vs Downstream
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    x = "Coefficient",
    y = "Gene",
    size = "-log10(adj.P.Val)",
    color = "logFC"
  )+
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),  # X-axis label
    axis.title.y = element_text(size = 30, face = "bold"),  # Y-axis label
    axis.text.y = element_text(size = 20),                 # gene names
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  # tilt X labels
    strip.text = element_text(size = 18)    # facet headings
  )

ggsave("D:/PLUS/bEV/results/Expression_WNt5a_TLR_pthway_genes.png",
       plot = q, width = 12, height = 15, dpi = 300)

####expression_of_these_sets_of_genes------------

dat.list <- list()
for(gg in goi_wnt5a_tlr_select$Gene.name){
  dat.list[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat <- bind_rows(dat.list, .id = "Gene.name")

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


t <- ggplot(dat, aes(x = donor, y = Gene.name, fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0
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


ggsave("D:/PLUS/bEV/results/Expression_Wnt5a_TLR2_4_pathway_selected_genes_batch_corrected.png",
       plot = t, width = 12, height = 15, dpi = 300)


#########intergrin_pathway_genes-------------

integrin_genes <- c(
  # Integrin receptors (α and β subunits)
  "ITGA1","ITGA2","ITGA3","ITGA4","ITGA5","ITGA6","ITGA7","ITGA8","ITGA9","ITGA10",
  "ITGA11","ITGAV","ITGAX","ITGAD","ITGAE","ITGAL","ITGAM","ITGAP","ITGA2B",
  "ITGB1","ITGB2","ITGB3","ITGB4","ITGB5","ITGB6","ITGB7","ITGB8",
  
  # Integrin ligands / ECM components
  "FN1","LAMA1","LAMA2","LAMA3","LAMA4","LAMA5","LAMC1","LAMC2",
  "COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL6A1",
  "COL6A2","COL6A3","COL7A1","THBS1","THBS2","THBS3","THBS4","VCAM1","ICAM1",
  
  # Adaptor proteins / scaffolds
  "TALIN1","TALIN2","VINC","FERMT1","FERMT2","FERMT3","PAXILLIN","PXN","FAK","PTK2",
  
  # Early kinases / signaling nodes
  "SRC","LYN","FAK","PTK2B","SYK","PIK3CA","PIK3CB","PIK3CD","PIK3R1","AKT1","AKT2","AKT3",
  "MAPK1","MAPK3","MAPK8","MAPK14","RAC1","RHOA","RHOC","CDC42","RAP1A","RAP1B",
  
  # Transcription factors / downstream effectors
  "JUN","FOS","NFKB1","RELA","CREB1","ELK1","SRF","YAP1","TAZ","EP300","SMAD3",
  
  # Cytokines / chemokines regulated downstream
  "IL6","IL8","CCL2","CCL5","CXCL1","CXCL8"
)

goi_intergrin_select <- limmaRes |>
  filter(coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  filter(Gene.name%in%integrin_genes)

show(goi_intergrin_select)

# Quick check
head(goi_intergrin_select)


u <- ggplot(goi_intergrin_select, aes(x = coef, y = reorder(Gene.name, logFC))) +
  geom_point(aes(size = -log10(adj.P.Val), color = logFC)) +
  #facet_wrap(~ layer, scales = "free_y") +  # separates Upstream vs Downstream
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    x = "Coefficient",
    y = "Gene",
    size = "-log10(adj.P.Val)",
    color = "logFC"
  )+
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),  # X-axis label
    axis.title.y = element_text(size = 30, face = "bold"),  # Y-axis label
    axis.text.y = element_text(size = 14),                 # gene names
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  # tilt X labels
    strip.text = element_text(size = 18)    # facet headings
  )

ggsave("D:/PLUS/bEV/results/Coefficients_Integrin_pathway_genes.png",
       plot = u, width = 12, height = 15, dpi = 300)

####expression_of_these_sets_of_genes------------

dat.list <- list()
for(gg in goi_intergrin_select$Gene.name){
  dat.list[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat <- bind_rows(dat.list, .id = "Gene.name")

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


v <- ggplot(dat, aes(x = donor, y = Gene.name, fill = Exp)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0
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
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 20),   # facet labels (mock/Hp/OMV)
    panel.spacing.x = unit(1, "lines")
  )


ggsave("D:/PLUS/bEV/results/Expression_Integrin_pathway_selected_genes_batch_corrected.png",
       plot = v, width = 12, height = 15, dpi = 300)



######Actin_regulating_genes---------------------------

actin_adhesion_genes <- c(
  # Core actin isoforms
  "ACTB", "ACTG1", "ACTA1", "ACTA2", "ACTC1", "ACTG2",
  
  # Actin nucleation / polymerization
  "ARP2","ARP3","ARPC1A","ARPC1B","ARPC2","ARPC3","ARPC4","ARPC5",
  "WASP","WASL","WAVE1","WAVE2","WAVE3",
  "DIAPH1","DIAPH2","DIAPH3",
  "FMNL1","FMNL2","FMNL3",
  "SPIRE1","SPIRE2",
  
  # Actin depolymerization / severing / capping
  "CFL1","CFL2","DSTN",
  "GELSOLIN","CAPZA1","CAPZA2","CAPZB","TWF1","TWF2",
  
  # Crosslinkers / bundlers
  "FILAMIN_A","FILAMIN_B","FILAMIN_C",
  "FSCN1","FSCN2","FSCN3",
  "TJP1","ACTN1","ACTN2","ACTN4",
  
  # Small GTPases (Rho family)
  "RHOA","RHOB","RHOC","RAC1","RAC2","RAC3","CDC42","RHOJ","RHOQ",
  
  # GEFs/GAPs regulating Rho/Rac/CDC42
  "TIAM1","VAV1","VAV2","VAV3",
  "ARHGEF1","ARHGEF2","ARHGEF6","ARHGAP1","ARHGAP5","ARHGAP22","ARHGAP21",
  "ARHGEF7","DOCK1","DOCK5",  # CDC42-specific GEFs
  "ARHGAP31","ARHGAP32","ARHGAP33",
  
  # Filopodia / lamellipodia regulators (CDC42 downstream)
  "EVL","MTSS1","MTSS2","MENA","MYO10","MYO15A","MYO7A","IRSp53","WIPF1",
  
  # Focal adhesion / adhesion signaling genes
  "TLN1","TLN2","VINC","PXN","FERMT1","FERMT2","FERMT3",
  "ITGB1","ITGB3","ITGA2","ITGA5",
  "FAK","PTK2B","SRC","ILK","PAXILLIN","ZYX","VASP","CRK","CRKL",
  
  # CDC42 pathway kinases
  "PAK1","PAK2","PAK3","MRCK","LIMK1","LIMK2","ROCK1","ROCK2","MAPK1","MAPK3","MAPK8","MAPK14",
  
  # Actin-associated motor proteins
  "MYH9","MYH10","MYH11","MYH14",
  "MYO1A","MYO1B","MYO1C","MYO1D","MYO1E","MYO1F",
  
  # Non-canonical Wnt / Rho-CDC42 interaction
  "ROR2","RYK","FZD2","FZD5","FZD7","DVL1","DVL2","DVL3",
  
  # Calcium / CAMK signaling (actin reorganization)
  "CAMK2A","CAMK2B","CAMK2D",
  
  # Polarity and scaffold proteins
  "PAR3","PAR6","CDC42EP1","CDC42EP2","CDC42EP3","CDC42EP4"
)
goi_act_select <- limmaRes |>
  filter(coef %in% c("OMVvsHp", "treatmentHp", "treatmentOMV"))|>
  filter(Gene.name%in%actin_adhesion_genes)

show(goi_act_select)

# Quick check
head(goi_act_select)


x <- ggplot(goi_act_select, aes(x = coef, y = reorder(Gene.name, logFC))) +
  geom_point(aes(size = -log10(adj.P.Val), color = logFC)) +
  #facet_wrap(~ layer, scales = "free_y") +  # separates Upstream vs Downstream
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    x = "Coefficient",
    y = "Gene",
    size = "-log10(adj.P.Val)",
    color = "logFC"
  )+
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),  # X-axis label
    axis.title.y = element_text(size = 30, face = "bold"),  # Y-axis label
    axis.text.y = element_text(size = 10),                 # gene names
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  # tilt X labels
    strip.text = element_text(size = 18)    # facet headings
  )

ggsave("D:/PLUS/bEV/results/Coefficients_actin_upstream_and_downstream_genes.png",
       plot = x, width = 12, height = 15, dpi = 300)

####expression_of_these_sets_of_genes------------

dat.list <- list()
for(gg in goi_act_select$Gene.name){
  dat.list[[gg]] <- metadata |>
    mutate(Exp=scale(expr_corrected[gg,])) |>
    remove_rownames()
}
dat <- bind_rows(dat.list, .id = "Gene.name")

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
    low = "blue", mid = "white", high = "red", midpoint = 0
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
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 20),   # facet labels (mock/Hp/OMV)
    panel.spacing.x = unit(1, "lines")
  )


ggsave("D:/PLUS/bEV/results/Expression_genes_upstream_and_downstream_actin_batch_corrected.png",
       plot = vc, width = 12, height = 15, dpi = 300)

