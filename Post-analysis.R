# ============================================
# CADRES POST-ANALYSIS PIPELINE (C>U vs A>I MUTATIONS)
# Expanded version with user-set input file, synchronized categories,
# and "Uncategorised" for missing locations
# ============================================

# ---- User-defined Input File ----
# Option 1: interactive (opens file explorer)
input_file <- file.choose()
# Option 2: manual path (uncomment and modify below)
# input_file <- "C:/Users/83935/Desktop/JoVE_Manuscript/293T_V6_rMATS-DVR_Result_annotated.txt"

# ---- Set Working Directory Automatically ----
setwd(dirname(input_file))
cat("\n📂 Working directory set to:", getwd(), "\n")

# ---- Load Required Packages ----
required_pkgs <- c("ggplot2", "dplyr", "readr", "stringr", "ggrepel", "forcats")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ---- Read Input Data ----
data <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE)
cat("\n✅ Data loaded from:", input_file, "\n")
cat("Number of rows:", nrow(data), " | Number of columns:", ncol(data), "\n\n")

# ============================================
# DATA PREPARATION
# ============================================

# ---- Compute Sequence Depth ----
data$Depth <- sapply(
  seq_len(nrow(data)),
  function(i) {
    alt_vals <- as.numeric(str_split(data$DMSO_Alt[i], ",")[[1]])
    ref_vals <- as.numeric(str_split(data$DMSO_Ref[i], ",")[[1]])
    mean(alt_vals + ref_vals, na.rm = TRUE)
  }
)

# ---- Handle FDR and Assign Mutation Labels ----
data$FDR[data$FDR == 0] <- 1e-5
data$log10FDR <- -log10(data$FDR)

data <- data %>%
  mutate(
    MutationCategory = case_when(
      Type == "C>T" ~ "C>U",
      Type == "A>G" ~ "A>I",
      TRUE ~ "Other"
    ),
    ColorGroup = case_when(
      FDR < 0.05 & MutationCategory == "C>U" ~ "C>U",
      FDR < 0.05 & MutationCategory == "A>I" ~ "A>I",
      FDR < 0.05 & MutationCategory == "Other" ~ "Other",
      TRUE ~ "Not Significant"
    ),
    Alt_allele_fraction_diff_flipped = -Alt_allele_fraction_diff
  )

# ============================================
# VOLCANO PLOT
# ============================================

top20 <- data %>%
  arrange(desc(log10FDR)) %>%
  slice_head(n = 20)

p_volcano <- ggplot(
  data,
  aes(x = Alt_allele_fraction_diff_flipped,
      y = log10FDR,
      color = ColorGroup,
      fill = ColorGroup,
      size = Depth)
) +
  geom_point(shape = 21, alpha = 0.8, stroke = 0.2, color = "grey30") +
  geom_text_repel(
    data = top20,
    aes(label = Genename),
    size = 5,
    fontface = "bold",
    color = "black",
    max.overlaps = 100
  ) +
  scale_fill_manual(values = c("C>U" = "#FF4D4D",
                               "A>I" = "#FFD300",
                               "Other" = "#0080FF",
                               "Not Significant" = "gray80")) +
  scale_color_manual(values = c("C>U" = "#FF4D4D",
                                "A>I" = "#FFD300",
                                "Other" = "#0080FF",
                                "Not Significant" = "gray80")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
  labs(
    title = "Volcano Plot of RNA Editing Differences (C>U vs A>I)",
    x = expression(-Delta ~ "Alt Allele Fraction (DOX - DMSO)"),
    y = expression(-log[10](FDR))
  ) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

ggsave("Plot_Volcano_CU_AI.png", p_volcano, width = 9, height = 6, dpi = 300)
cat("💾 Saved: Plot_Volcano_CU_AI.png\n")

# ============================================
# LOCATION COMPARISON (C>U vs A>I)
# ============================================

data_sig_05 <- data %>% filter(FDR < 0.05)

# ---- Replace missing/empty Location values with "Uncategorised" ----
data_sig_05$Location[data_sig_05$Location == "" | is.na(data_sig_05$Location)] <- "Uncategorised"

# ---- Ensure both mutation types share all possible Location categories ----
all_locations <- unique(data_sig_05$Location)

location_summary <- data_sig_05 %>%
  filter(MutationCategory %in% c("C>U", "A>I")) %>%
  group_by(MutationCategory, Location) %>%
  summarise(n = n(), .groups = "drop")

# ---- Create a complete grid for all MutationCategory-Location pairs ----
complete_table <- expand.grid(
  MutationCategory = c("C>U", "A>I"),
  Location = all_locations
)

# ---- Merge and fill missing combinations with zero ----
location_summary <- full_join(location_summary, complete_table,
                              by = c("MutationCategory", "Location"))
location_summary$n[is.na(location_summary$n)] <- 0

# ---- Compute percentage per mutation type ----
location_summary <- location_summary %>%
  group_by(MutationCategory) %>%
  mutate(Percentage = 100 * n / sum(n, na.rm = TRUE)) %>%
  ungroup()

# ---- Generate plot ----
p_location <- ggplot(location_summary,
                     aes(
                       x = fct_reorder(Location, -Percentage),
                       y = Percentage,
                       fill = MutationCategory
                     )) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Distribution of Mutation Locations (C>U vs A>I, FDR < 0.05)",
    x = "Genomic Location",
    y = "Percentage of Sites",
    fill = "Mutation Type"
  ) +
  scale_fill_manual(values = c("C>U" = "#FF4D4D",
                               "A>I" = "#FFD300")) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Plot_Location_Comparison_CU_AI.png", p_location, width = 9, height = 6, dpi = 300)
cat("💾 Saved: Plot_Location_Comparison_CU_AI.png (missing → Uncategorised)\n")

# ============================================
# EDITING SHIFT PER MUTATION TYPE
# ============================================

filtered_data <- data %>%
  filter(FDR < 0.05) %>%
  filter(!is.na(Alt_allele_fraction_diff_flipped)) %>%
  filter(MutationCategory %in% c("C>U", "A>I"))

p_shift <- ggplot(filtered_data,
                  aes(x = MutationCategory,
                      y = Alt_allele_fraction_diff_flipped,
                      fill = MutationCategory)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.4) +
  scale_fill_manual(values = c("C>U" = "#FF4D4D", "A>I" = "#FFD300")) +
  labs(
    title = "Editing Shift per Mutation Type (FDR < 0.05)",
    x = "Mutation Type",
    y = expression(-Delta ~ "Alt Allele Fraction (DOX - DMSO)")
  ) +
  theme_classic(base_size = 18)

ggsave("Plot_EditingShift_CU_AI.png", p_shift, width = 9, height = 6, dpi = 300)
cat("💾 Saved: Plot_EditingShift_CU_AI.png\n")

# ============================================
# TOP GENES WITH MULTIPLE EDITING SITES
# ============================================

gene_summary <- data %>%
  filter(FDR < 0.05, Genename != "") %>%
  group_by(Genename) %>%
  summarise(
    n_sites = n(),
    mean_diff = mean(Alt_allele_fraction_diff, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE)
  ) %>%
  arrange(desc(n_sites))

top_genes <- head(gene_summary, 20)

p_genes <- ggplot(top_genes,
                  aes(x = reorder(Genename, n_sites),
                      y = n_sites,
                      fill = mean_diff)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "#377EB8",
                       mid = "white",
                       high = "#E41A1C",
                       midpoint = 0) +
  labs(
    title = "Top Genes with Multiple Editing Sites (FDR < 0.05)",
    x = "Gene",
    y = "Number of Sites",
    fill = "Mean Δ Fraction"
  ) +
  theme_classic(base_size = 18)

ggsave("Plot_TopGenes_CU_AI.png", p_genes, width = 9, height = 6, dpi = 300)
cat("💾 Saved: Plot_TopGenes_CU_AI.png\n")

# ============================================
# EXPORT SIGNIFICANT GENES
# ============================================

sig_genes <- unique(data$Genename[data$FDR < 0.05 & data$Genename != ""])
write.table(sig_genes,
            file = "significant_genes_FDR0.05.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

cat("🧬 Wrote significant genes to: significant_genes_FDR0.05.txt\n")

# ============================================
# EDITING SITE DISTRIBUTION BY LOCATION
# ============================================

distribution <- data %>%
  mutate(Location = ifelse(Location == "" | is.na(Location), "Uncategorised", Location)) %>%
  filter(FDR < 0.05) %>%
  count(Location) %>%
  mutate(Percent = 100 * n / sum(n))

p_dist <- ggplot(distribution,
                 aes(x = fct_reorder(Location, Percent),
                     y = Percent,
                     fill = Location)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "Editing Site Distribution by Genomic Location (FDR < 0.05)",
    x = "Genomic Region",
    y = "Percentage of Sites"
  ) +
  theme_classic(base_size = 18)

ggsave("Plot_LocationDistribution_AllSig.png", p_dist, width = 9, height = 6, dpi = 300)
cat("💾 Saved: Plot_LocationDistribution_AllSig.png\n")

# ============================================
# DOX vs DMSO CORRELATION
# ============================================

get_mean <- function(x) {
  vals <- strsplit(as.character(x), ",")[[1]]
  mean(as.numeric(vals), na.rm = TRUE)
}

data_rep <- data %>%
  mutate(DMSO_mean = sapply(DMSO_Alt_allele_fraction, get_mean)) %>%
  mutate(DOX_mean  = sapply(DOX_Alt_allele_fraction, get_mean)) %>%
  mutate(
    DMSO_depth = sapply(DMSO_Ref, function(x)
      mean(as.numeric(strsplit(as.character(x), ",")[[1]]), na.rm = TRUE)) +
      sapply(DMSO_Alt, function(x)
        mean(as.numeric(strsplit(as.character(x), ",")[[1]]), na.rm = TRUE)),
    DOX_depth = sapply(DOX_Ref, function(x)
      mean(as.numeric(strsplit(as.character(x), ",")[[1]]), na.rm = TRUE)) +
      sapply(DOX_Alt, function(x)
        mean(as.numeric(strsplit(as.character(x), ",")[[1]]), na.rm = TRUE))
  ) %>%
  mutate(Avg_depth = (DMSO_depth + DOX_depth) / 2)

cor_test_res <- cor.test(data_rep$DMSO_mean, data_rep$DOX_mean, method = "pearson")

cat("\n===== Editing Correlation Summary =====\n")
cat("Pearson's r :", round(cor_test_res$estimate, 3), "\n")
cat("P-value     :", formatC(cor_test_res$p.value, format = "e", digits = 2), "\n")
cat("N sites     :", nrow(data_rep), "\n")

p_corr <- ggplot(data_rep,
                 aes(x = DMSO_mean,
                     y = DOX_mean,
                     color = MutationCategory,
                     alpha = Avg_depth)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("C>U" = "#FF4D4D",
                                "A>I" = "#FFD300",
                                "Other" = "#0080FF")) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = sprintf("DOX vs DMSO Editing Correlation (r = %.2f, p = %.2e)",
                    cor_test_res$estimate, cor_test_res$p.value),
    x = "Mean DMSO Alt Allele Fraction",
    y = "Mean DOX Alt Allele Fraction",
    color = "Mutation Type",
    alpha = "Mean Read Depth"
  ) +
  theme_classic(base_size = 18)

ggsave("Plot_Correlation_DMSO_DOX.png", p_corr, width = 9, height = 6, dpi = 300)
cat("💾 Saved: Plot_Correlation_DMSO_DOX.png\n")

cat("\n🎉 All analyses completed successfully with missing locations classified as 'Uncategorised'!\n")

