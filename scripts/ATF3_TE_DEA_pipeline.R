############################################################
## ATF3 KO vs Control — Transposable Element (TE) Analysis
############################################################

## ---------------------------
## 1. Load packages
## ---------------------------
library(dplyr)
library(tidyr)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(forcats)

## ---------------------------
## 2. Paths and basic settings
## ---------------------------
# Folder that contains the *.cntTable files (ATF3 experiment)
base_dir  <- "D:/LAB/Walaa_RNA_seq/atf3"

# Output directory for DEA results and plots
output_dir <- file.path(base_dir, "DEA")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## ---------------------------
## 3. Read all *.cntTable files and build count matrix
## ---------------------------

# Each .cntTable is assumed to have:
#   column 1: feature (gene or TE name)
#   column 2: raw counts
files <- list.files(base_dir, pattern = "\\.cntTable$", full.names = TRUE)
files

read_cnt <- function(path) {
  tab <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  data.frame(
    feature = tab[[1]],
    count   = tab[[2]],
    stringsAsFactors = FALSE
  )
}

cnt_list <- lapply(files, read_cnt)

# Sample names from file names (remove .cntTable)
sample_names <- basename(files) %>% sub("\\.cntTable$", "", .)
names(cnt_list) <- sample_names

# Rename "count" column in each table to the sample name
cnt_list_renamed <- Map(function(tab, samp){
  names(tab)[2] <- samp
  tab
}, cnt_list, sample_names)

# Merge all tables by "feature"
cnt_wide <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE),
                   cnt_list_renamed)

# Turn NA to 0
cnt_wide[is.na(cnt_wide)] <- 0

# Convert to matrix with rownames = feature
cnt_wide_df <- as.data.frame(cnt_wide)
rownames(cnt_wide_df) <- cnt_wide_df$feature

count_mat <- as.matrix(cnt_wide_df[ , -1, drop = FALSE])

dim(count_mat)
head(rownames(count_mat))
colnames(count_mat)

## ---------------------------
## 4. Split genes vs TEs
## ---------------------------
# Here we assume gene IDs start with "ENSG"
is_gene <- grepl("^ENSG", rownames(count_mat))

gene_counts <- count_mat[is_gene, , drop = FALSE]
TE_counts   <- count_mat[!is_gene, , drop = FALSE]

dim(gene_counts)
dim(TE_counts)
head(rownames(TE_counts))

## ---------------------------
## 5. Define sample groups (ATF3_KO vs Control)
## ---------------------------
# IMPORTANT: change these to your actual sample names if needed
group <- rep(NA_character_, ncol(count_mat))
names(group) <- colnames(count_mat)

# ATF3 knockout samples
group[c("105_ATF3_8_1_S10_R1_001",
        "105_ATF3_8_2_S11_R1_001",
        "105_ATF3_8_3_S12_R1_001")] <- "atf3_ko"

# Control samples
group[c("105_CONTROL_2_S37_R1_001",
        "105_CONTROL_3_S38_R1_001",
        "105_CONTROL_1_S24_R1_001")] <- "controls"

group
group <- factor(group, levels = c("controls", "atf3_ko"))  # controls as reference
group

## ---------------------------
## 6. edgeR object for TEs + filtering + normalization
## ---------------------------
y_TE <- DGEList(counts = TE_counts, group = group)
y_TE

# Filter lowly expressed TEs
keep_TE <- filterByExpr(y_TE, group = group)
sum(keep_TE)   # number of retained TEs

y_TE <- y_TE[keep_TE, , keep.lib.sizes = FALSE]

# TMM normalization
y_TE <- calcNormFactors(y_TE)
logCPM <- cpm(y_TE, log = TRUE)

# Quick QC: TE counts range
summary(as.numeric(TE_counts))

## ---------------------------
## 7. PCA on normalized TE counts
## ---------------------------
pca <- prcomp(t(logCPM), scale. = TRUE)   # samples as rows

var_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100
round(var_explained[1:6], 2)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = pca$x[,3],
  PC4 = pca$x[,4],
  group = group
)

# Optional: extract some ID from sample name (here just dummy)
pca_df$sample_id <- pca_df$Sample

p12 <- ggplot(pca_df, aes(PC1, PC2, color = group, label = sample_id)) +
  geom_point(size = 6) +
  geom_text(vjust = -0.6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of TE expression (PC1 vs PC2)",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  )

p34 <- ggplot(pca_df, aes(PC3, PC4, color = group, label = sample_id)) +
  geom_point(size = 6) +
  geom_text(vjust = -0.6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of TE expression (PC3 vs PC4)",
    x = paste0("PC3 (", round(var_explained[3], 1), "%)"),
    y = paste0("PC4 (", round(var_explained[4], 1), "%)")
  )

p_pca <- p12 | p34
p_pca

ggsave(file.path(output_dir, "TE_PCA_ATF3KO_vs_Control.png"),
       plot = p_pca, width = 10, height = 5, dpi = 300)

## ---------------------------
## 8. Differential expression for TEs (edgeR GLM)
## ---------------------------
design <- model.matrix(~ group)
design

y_TE <- estimateDisp(y_TE, design)
plotBCV(y_TE)

fit_TE <- glmFit(y_TE, design)
# coef = 2 corresponds to "groupatf3_ko" (ATF3_KO vs controls)
lrt_TE <- glmLRT(fit_TE, coef = 2)

res_TE <- topTags(lrt_TE, n = Inf)$table
head(res_TE)
sum(res_TE$FDR < 0.05)

# Add TE name as a column
res_TE$TE <- rownames(res_TE)

# Save DEA table
write.csv(res_TE,
          file = file.path(output_dir, "TE_DEA_results_ATF3KO_vs_Control.csv"),
          row.names = TRUE, quote = FALSE)

## ---------------------------
## 9. Volcano plot (nice, labeled)
## ---------------------------
res_TE$significant <- res_TE$FDR < 0.05

res_TE$direction <- ifelse(
  res_TE$FDR < 0.05 & res_TE$logFC > 0,  "Up in ATF3_KO",
  ifelse(res_TE$FDR < 0.05 & res_TE$logFC < 0, "Down in ATF3_KO", "Not significant")
)

volcano_df <- res_TE %>%
  mutate(
    neg_logFDR = -log10(FDR)
  )

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = neg_logFDR)) +
  geom_point(aes(color = direction),
             alpha = 0.8, size = 1.8) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_color_manual(
    values = c(
      "Up in ATF3_KO"   = "#d73027",  # red
      "Down in ATF3_KO" = "#4575b4",  # blue
      "Not significant" = "grey80"
    )
  ) +
  labs(
    title = "Differential Expression — Transposable Elements",
    x = "Log2 Fold Change (ATF3_KO / Control)",
    y = expression(-log[10]("FDR")),
    color = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.position = "top",
    legend.margin = margin(b = 5),
    legend.box.margin = margin(b = -5)
  )

# Label top significant TEs
top_labels <- volcano_df %>%
  filter(significant) %>%
  arrange(desc(neg_logFDR)) %>%
  slice_head(n = 10)   # top 10 – change if you want

p_volcano_labeled <- p_volcano +
  geom_text_repel(
    data = top_labels,
    aes(label = TE),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.4,
    point.padding = 0.2,
    min.segment.length = 0
  )

p_volcano_labeled

ggsave(file.path(output_dir, "TE_volcano_ATF3KO_vs_Control_labeled.png"),
       plot = p_volcano_labeled,
       width = 6, height = 5, dpi = 300)

## ---------------------------
## 10. Parse TE into subfamily / family / class
## ---------------------------
# TE names look like: "MER39B:ERV1:LTR" or "L1M2a:L1:LINE"
res_TE$split     <- strsplit(res_TE$TE, ":", fixed = TRUE)
res_TE$subfamily <- sapply(res_TE$split, `[`, 1)
res_TE$family    <- sapply(res_TE$split, `[`, 2)
res_TE$class     <- sapply(res_TE$split, `[`, 3)

table(res_TE$family[res_TE$FDR < 0.05])
table(res_TE$class[res_TE$FDR < 0.05])

## ---------------------------
## 11. “GSEA-style” TE family enrichment (Wilcoxon rank test)
## ---------------------------
# Rank statistic: sign(logFC) * -log10(PValue)
rank_stats <- res_TE %>%
  mutate(rank_score = sign(logFC) * (-log10(PValue))) %>%
  select(TE, class, family, subfamily, rank_score)

# Wilcoxon per family
gsea_results <- rank_stats %>%
  group_by(family) %>%
  summarise(
    n = n(),
    stat = wilcox.test(rank_score, alternative = "two.sided")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(stat, method = "BH")
  ) %>%
  arrange(FDR)

head(gsea_results, 10)

# Also per TE class (LINE, LTR, SINE, etc.)
gsea_class <- rank_stats %>%
  group_by(class) %>%
  summarise(
    n = n(),
    stat = wilcox.test(rank_score, alternative = "two.sided")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(stat, method = "BH")
  ) %>%
  arrange(FDR)

gsea_class

# Save GSEA-style results
write.csv(gsea_results,
          file = file.path(output_dir, "TE_family_enrichment_ATF3KO_vs_Control.csv"),
          row.names = FALSE, quote = FALSE)

write.csv(gsea_class,
          file = file.path(output_dir, "TE_class_enrichment_ATF3KO_vs_Control.csv"),
          row.names = FALSE, quote = FALSE)

## ---------------------------
## 12. Plot top enriched TE families
## ---------------------------



# compute direction for each family
family_dirs <- rank_stats %>%
  group_by(family) %>%
  summarize(
    median_rank = median(rank_score),
    .groups = "drop"
  )

gsea_results2 <- gsea_results %>%
  left_join(family_dirs, by = "family") %>%
  mutate(
    direction = case_when(
      FDR < 0.05 & median_rank > 0 ~ "Enriched UP in ATF3_KO",
      FDR < 0.05 & median_rank < 0 ~ "Enriched DOWN in ATF3_KO",
      TRUE ~ "Not significant"
    )
  )

gsea_results2

library(dplyr)
library(ggplot2)
library(forcats)

# take top 15 families or change n = 20 etc.
plot_df <- gsea_results2 %>%
  slice_min(FDR, n = 15) %>%
  mutate(
    family = fct_reorder(family, -log10(FDR)),
    dir_color = case_when(
      direction == "Enriched UP in ATF3_KO"   ~ "UP",
      direction == "Enriched DOWN in ATF3_KO" ~ "DOWN",
      TRUE                                    ~ "NS"
    )
  )

# choose colors
cols <- c(
  "UP"   = "#d73027",   # red
  "DOWN" = "#4575b4",   # blue
  "NS"   = "grey80"     # light gray
)

p_enrich <- ggplot(plot_df, aes(x = family, y = -log10(FDR), fill = dir_color)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "black", linewidth = 0.6) +
  coord_flip() +
  scale_fill_manual(values = cols, name = "Direction",
                    labels = c(
                      "UP"   = "Enriched in ATF3_KO",
                      "DOWN" = "Depleted in ATF3_KO",
                      "NS"   = "Not significant"
                    )) +
  labs(
    x = "TE family",
    y = expression(-log[10]("FDR")),
    title = "TE Family Enrichment in ATF3_KO (Wilcoxon rank test)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    axis.title.y  = element_text(margin = margin(r = 10)),
    axis.title.x  = element_text(margin = margin(t = 10)),
    legend.position = "top"
  )

p_enrich
ggsave(file.path(output_dir, "TE_family_enrichment_direction_ATF3KO.png"),
       plot = p_enrich, width = 7, height = 5, dpi = 300)
