
#load packages

library(tidyverse)
library(edgeR)
library(limma)
library(pheatmap)
library(here)

####createdirectories-------------------

Indir <- here("rawdata/")
outdir <- here("results/beV_RNAseq_monocytes_limma_01")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)



####loaddata------------------------

counts <- read.csv(
  here("rawdata", "RNAseq_raw_counts 1.csv"),
  row.names = 1
)


#####grouping by gene names and take average of genes when there are mutiple---

df_mean <- counts |>
  group_by(Gene.name) |>
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
counts <- df_mean |> na.omit()
counts <- as.data.frame(counts)
rownames(counts) <- counts$Gene.name
counts$Gene.name <- NULL


#create and load metadata---------
## Sample column names (replace with colnames(your_data)[1:9])
samples <- colnames(counts)

## Build metadata by splitting on "."
metadata <- data.frame(
  sample    = samples,
  donor     = sub("\\..*$", "", samples),       # text before '.'
  treatment = sub("^.*\\.", "", samples)        # text after '.'
)

saveRDS(
  metadata,
  file = file.path(outdir, "metadata.rds")
)



#For correlation matrix you need a datamatrix with only number hence removing the first column

corMT <- cor(counts, method = "spearman")
diag(corMT) <- NA   # makes diagonal NAs 

png(file.path(outdir, "heatmapcorMT.png"),
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
ggsave(
  filename = file.path(outdir, "mds_plot.png"),
  plot = mds_plot,
  width = 8,
  height = 6,
  dpi = 300
)



# Relevel treatment so "mock" is baseline

metadata <- metadata |>                      
  mutate(treatment = fct_relevel(treatment, "mock"))

# Create design matrix
mm <- model.matrix(~ donor + treatment, data = metadata)

# Set row names to sample IDs so they appear on y-axis
# Replace 'sample' with the actual column name in metadata containing sample names
rownames(mm) <- metadata$sample

# Save heatmap as PNG

png(file.path(outdir, "design_matrix.png"),
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
dataVoom <- voom(dge, design = mm, plot = FALSE)
show(dataVoom)

#defining a statistical model for DGE
fit <- lmFit(dataVoom, design = mm)
fit <- eBayes(fit) #uses bayesian moderation 
show(coef(fit))

limmaRes <- list() # start an empty list
for(coefx in colnames(coef(fit))){ # run a loop for each coefficient/name coefx chosen by us
  print(coefx)
  # topTable returns the statistics of our genes. We then store the result of each coefficient in a list.
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


saveRDS(
  limmaRes,
  file = file.path(outdir, "limmaRes.rds")
)





# data visulaization ------------------------------------------------------

# ---- Volcano plot (adj.P.Val) ----
volcano_adj <- limmaRes |>
  ggplot(aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(alpha = 0.3) +
  facet_grid(cols = vars(coef))

ggsave(
  filename = file.path(outdir, "volcano_adjPval.png"),
  plot = volcano_adj,
  width = 10,
  height = 6,
  dpi = 300
)

# ---- Volcano hex plot ----
volcano_hex <- limmaRes |>
  ggplot(aes(x = logFC, y = -log10(P.Value))) + 
  geom_hex() +
  facet_grid(cols = vars(coef))

ggsave(
  filename = file.path(outdir, "volcano_hex.png"),
  plot = volcano_hex,
  width = 10,
  height = 6,
  dpi = 300
)

# ---- P-value distribution ----
pval_hist <- limmaRes |>
  ggplot(aes(x = P.Value, fill = factor(floor(AveExpr)))) +
  geom_histogram() +
  facet_grid(cols = vars(coef))

ggsave(
  filename = file.path(outdir, "pval_histogram.png"),
  plot = pval_hist,
  width = 10,
  height = 6,
  dpi = 300
)


# ---- Extract significant genes ----
limmaResSig <- limmaRes |> 
  filter(adj.P.Val < 0.05)

count(limmaResSig, coef)

saveRDS(
  limmaResSig,
  file = file.path(outdir, "limmaRes.rds")
)


#-------------OMVvsHp and OMVvsmock enriched and absent in HP vs mock-----------------

################removing batch effect--------------

design2 <- model.matrix(~ treatment, data = metadata)
expr_corrected <- removeBatchEffect(dataVoom$E, batch = metadata$donor, 
                                    design = design2)


saveRDS(
  expr_corrected,
  file = file.path(outdir, "expr_corrected.rds")
)

########identiyingBEVspecificgenes--------------------------------

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

saveRDS(
  genes_OMV_specific_up,
  file = file.path(outdir, "limmaRes.rds")
)

saveRDS(
  genes_OMV_specific_down,
  file = file.path(outdir, "limmaRes.rds")
)



######plottingexpressionforDEGs------------------------


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
    low = "#0571b0",
    mid = "#F7F7F7",
    high = "#ca0020",
    midpoint = 0,
    name = "Scaled expression"
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

ggsave(
  filename = file.path(outdir, "Expression_of_OMV_specific_upregulated_DEGs_batchcorrected.png"),
  plot = d,
  width = 12,
  height = 15,
  dpi = 300
)


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
  low = "#0571b0",
  mid = "#F7F7F7",
  high = "#ca0020",
  midpoint = 0,
  name = "Scaled expression"
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

ggsave(
  filename = file.path(outdir, "Expression_of_OMV_specific_downregulated_DEGs_batchcorrected.png"),
  plot = d,
  width = 12,
  height = 15,
  dpi = 300
)





