# ==============================================================================
# Taxonomic Analysis: Omni vs Vegan Diet Microbiome
# Analyses:
#   A. Taxonomic abundance (bar plots: Kraken2 vs Bracken)
#   B. Alpha diversity (Chao1, Berger-Parker, Shannon)
#   C. Beta diversity (Bray-Curtis + PERMANOVA)
#   D. Differential abundance (ALDEx2)
# ==============================================================================
#LOAD PACKAGES -----------------------------------------------------------------
library(dplyr)
library(phyloseq)
library(readr)
library(ALDEx2)
library(conflicted)
library(tidyverse)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggrepel)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(purrr::reduce)

options(error = NULL)

#SET PATHS ---------------------------------------------------------------------
data_dir <- "/Users/ailishm/Desktop/Kraken2_Ailish"
metadata_path <- "/Users/ailishm/Desktop/Kraken2_Ailish/metadata.csv"

#LOAD METADATA -----------------------------------------------------------------
metadata <- read_csv(metadata_path) %>%
  column_to_rownames("SampleID")
metadata$Diet <- factor(metadata$Diet, levels = c("Omni", "Vegan"))
srr_to_sample <- setNames(rownames(metadata), metadata$SRR)

#Check
cat("Metadata loaded:", nrow(metadata), "samples\n")

#COLOUR PALETTES ----------------------------------------------------------------
diet_colors <- c("Omni" = "orchid", "Vegan" = "lightgreen")
top_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(20),
  paste0("sp",1:20))      

#READ FUNCTIONS -----------------------------------------------------------
read_bracken_file <- function(filepath, sample_name) {
  df <- read_tsv(filepath, col_types = cols(), show_col_types = FALSE)
  df <- df[df$taxonomy_lvl == "S", ]
  out <- data.frame(
    name = df$name,
    taxonomy_id = as.character(df$taxonomy_id),
    stringsAsFactors = FALSE)
  out[[sample_name]] <- df$new_est_reads
  return(out)
}

read_kraken_file <- function(filepath, sample_name) {
  df <- read_tsv(filepath, col_types = cols(), show_col_types = FALSE,
                 col_names = c("pct","clade_reads","direct_reads",
                               "rank","taxid","name"))
  df$name <- trimws(df$name)
  df <- df[df$rank == "S", ]
  out <- data.frame(
    name = df$name,
    taxonomy_id = as.character(df$taxid),
    stringsAsFactors = FALSE)
  out[[sample_name]] <- df$direct_reads
  return(out)
}

#BUILD ABUNDANCE TABLES --------------------------------------------------------

## BRACKEN
bracken_files  <- list.files(data_dir, pattern = "\\.bracken$", 
                             full.names = TRUE)
bracken_list   <- lapply(bracken_files, function(f) {
  srr <- sub("\\.bracken$", "", basename(f))
  read_bracken_file(f, srr_to_sample[[srr]])
})
bracken_merged <- purrr::reduce(bracken_list,
                                function(a, b) merge(a, b,
                                                     by = c("name","taxonomy_id"), 
                                                     all = TRUE))
bracken_merged[is.na(bracken_merged)] <- 0

## KRAKEN (species level only)
kraken_files  <- list.files(data_dir, pattern = "\\.kraken\\.report$", 
                            full.names = TRUE)
kraken_list   <- lapply(kraken_files, function(f) {
  srr <- sub("\\.kraken\\.report$", "", basename(f))
  read_kraken_file(f, srr_to_sample[[srr]])
})
kraken_merged <- purrr::reduce(kraken_list,
                               function(a, b) merge(a, b,
                                                    by = c("name","taxonomy_id"), 
                                                    all = TRUE))
kraken_merged[is.na(kraken_merged)] <- 0

#Check species detected
cat("Bracken species detected:", nrow(bracken_merged), "\n")
cat("Kraken  species detected:", nrow(kraken_merged),  "\n")

#BUILD PHYLOSEQ OBJECTS --------------------------------------------------------
build_phyloseq <- function(count_df, meta) {
  # OTU matrix: rows = taxa, cols = samples
  sample_cols <- rownames(meta)
  otu_mat <- as.matrix(count_df[, sample_cols])
  rownames(otu_mat) <- count_df$taxonomy_id
  
  # Taxonomy matrix
  tax_mat <- matrix(count_df$name, ncol = 1,
                    dimnames = list(count_df$taxonomy_id, "Species"))
  
  phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_mat),
    sample_data(meta)
  )
}

ps_bracken <- build_phyloseq(bracken_merged, metadata)
ps_kraken  <- build_phyloseq(kraken_merged,  metadata)

#Filter to keep taxa with >10 reads total
ps_bracken <- prune_taxa(taxa_sums(ps_bracken) > 10, ps_bracken)
ps_kraken  <- prune_taxa(taxa_sums(ps_kraken)  > 10, ps_kraken)

#Check species detected after filtering
cat("Bracken species after filtering:", ntaxa(ps_bracken), "\n")
cat("Kraken  species after filtering:", ntaxa(ps_kraken),  "\n")

# ==============================================================================
# SECTION A1: RELATIVE ABUNDANCE - SPECIES LEVEL (Kraken2 vs Bracken)
# ==============================================================================

make_barplot <- function(ps, title_label) {
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  #identify top 40 taxa by mean relative abundance across all samples
  top40_ids <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:40]
  #get full OTU table (taxa x samples) and species name lookup
  otu_full  <- as.data.frame(otu_table(ps_rel))         
  tax_vec   <- as.vector(tax_table(ps_rel)[, "Species"])
  names(tax_vec) <- taxa_names(ps_rel)
  #label taxa: top 40 keep their name, rest become "Other"
  otu_full$Species <- ifelse(rownames(otu_full) %in% top40_ids,
                             tax_vec[rownames(otu_full)],
                             "Other")
  #collapse all "Other" rows by summing per sample
  sample_cols <- rownames(metadata)
  df_agg <- otu_full %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(across(all_of(sample_cols), sum), .groups = "drop")
  #pivot to long
  df_long <- df_agg %>%
    pivot_longer(cols      = all_of(sample_cols),
                 names_to  = "Sample",
                 values_to = "Abundance") %>%
    mutate(Sample = factor(Sample, levels = rownames(metadata)),
           Diet   = factor(metadata[as.character(Sample), "Diet"],
                           levels = c("Omni", "Vegan")))
  #colours: top 20 get distinct colours, Other gets grey
  sp_names    <- setdiff(unique(df_long$Species), "Other")
  sp_colors   <- setNames(
    colorRampPalette(brewer.pal(12, "Paired"))(length(sp_names)),
    sp_names
  )
  sp_colors["Other"] <- "grey75"
  #order legend: top 40 alphabetically, then Other at bottom
  legend_order <- c(sort(sp_names), "Other")
  df_long$Species <- factor(df_long$Species, levels = legend_order)
  
  ggplot(df_long, aes(x = Sample, y = Abundance, fill = Species)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.15) +
    scale_fill_manual(values = sp_colors) +
    facet_wrap(~Diet, scales = "free_x", nrow = 1) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      legend.text      = element_text(size = 6.5),
      legend.key.size  = unit(0.35, "cm"),
      plot.title       = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "grey90")
    ) +
    labs(title = title_label,
         x = "Sample", y = "Relative Abundance", fill = "Species")
}

bar_kraken  <- make_barplot(ps_kraken,  "Kraken2 — Top 40 Species")
bar_bracken <- make_barplot(ps_bracken, "Bracken  — Top 40 Species")

combined_bar <- bar_kraken / bar_bracken +
  plot_annotation(
    title    = "Taxonomic Abundance: Kraken2 vs Bracken",
    subtitle = "Omni (n=3) vs Vegan (n=3) | Top 40 species by relative abundance",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

print(combined_bar)
ggsave(file.path(data_dir, "A1_abundance_barplot.png"),
       combined_bar, width = 16, height = 12)
cat("Saved: A1_abundance_barplot.png\n")

# ==============================================================================
# SECTION A2: RELATIVE ABUNDANCE — GENUS LEVEL (Kraken2 vs Bracken)
# ==============================================================================

make_genus_barplot <- function(count_df, source_label, title_label) {
  sample_cols <- rownames(metadata)
  #for kraken: pull genus-level rows directly from raw report files
  if (source_label == "Kraken2") {
    genus_list <- lapply(kraken_files, function(f) {
      srr         <- sub("\\.kraken\\.report$", "", basename(f))
      sample_name <- srr_to_sample[[srr]]
      df <- read_tsv(f, col_types = cols(), show_col_types = FALSE,
                     col_names = c("pct","clade_reads","direct_reads",
                                   "rank","taxid","name"))
      df$name <- trimws(df$name)
      df <- df[df$rank == "G", ]
      out <- data.frame(Genus = df$name, stringsAsFactors = FALSE)
      out[[sample_name]] <- df$direct_reads
      return(out)
    })
    genus_merged <- purrr::reduce(genus_list,
                                  function(a, b) merge(a, b, by = "Genus", all = TRUE))
    genus_merged[is.na(genus_merged)] <- 0
    
  } else {
    #for bracken: derive genus by taking first two words of species name
    genus_merged <- count_df
    genus_merged$Genus <- sapply(strsplit(count_df$name, " "), `[`, 1)
    genus_merged <- genus_merged %>%
      dplyr::select(-name, -taxonomy_id) %>%
      dplyr::group_by(Genus) %>%
      dplyr::summarise(across(all_of(sample_cols), sum), .groups = "drop") %>%
      as.data.frame()
  }
  #relative abundance per sample
  num_cols <- sample_cols
  totals   <- colSums(genus_merged[, num_cols])
  rel_mat  <- sweep(genus_merged[, num_cols], 2, totals, "/")
  genus_merged[, num_cols] <- rel_mat
  #top 20 genera by mean relative abundance
  genus_merged$mean_ab <- rowMeans(genus_merged[, num_cols])
  top20_genera <- genus_merged$Genus[order(genus_merged$mean_ab,
                                           decreasing = TRUE)][1:20]
  #label: top 20 keep name, rest = "Other"
  genus_merged$Label <- ifelse(genus_merged$Genus %in% top20_genera,
                               genus_merged$Genus, "Other")
  #aggregate Other
  df_agg <- genus_merged %>%
    dplyr::select(-Genus, -mean_ab) %>%
    dplyr::group_by(Label) %>%
    dplyr::summarise(across(all_of(sample_cols), sum), .groups = "drop")
  #long
  df_long <- df_agg %>%
    pivot_longer(cols      = all_of(sample_cols),
                 names_to  = "Sample",
                 values_to = "Abundance") %>%
    mutate(Sample = factor(Sample, levels = rownames(metadata)),
           Diet   = factor(metadata[as.character(Sample), "Diet"],
                           levels = c("Omni", "Vegan")))
  #colours
  genus_names   <- setdiff(unique(df_long$Label), "Other")
  genus_colors  <- setNames(
    colorRampPalette(brewer.pal(12, "Paired"))(length(genus_names)),
    genus_names
  )
  genus_colors["Other"] <- "grey75"
  legend_order  <- c(sort(genus_names), "Other")
  df_long$Label <- factor(df_long$Label, levels = legend_order)
  
  ggplot(df_long, aes(x = Sample, y = Abundance, fill = Label)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.15) +
    scale_fill_manual(values = genus_colors) +
    facet_wrap(~Diet, scales = "free_x", nrow = 1) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      legend.text      = element_text(size = 6.5),
      legend.key.size  = unit(0.35, "cm"),
      plot.title       = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "grey90")
    ) +
    labs(title = title_label,
         x = "Sample", y = "Relative Abundance", fill = "Genus")
}

genus_kraken  <- make_genus_barplot(kraken_merged,  "Kraken2",
                                    "Kraken2 — Top 20 Genera")
genus_bracken <- make_genus_barplot(bracken_merged, "Bracken",
                                    "Bracken  — Top 20 Genera")

combined_genus <- genus_kraken / genus_bracken +
  plot_annotation(
    title    = "Genus-Level Relative Abundance: Kraken2 vs Bracken",
    subtitle = "Omni (n=3) vs Vegan (n=3) | Top 20 genera, remainder as Other",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    ))

print(combined_genus)
ggsave(file.path(data_dir, "A2_genus_abundance_barplot.png"),
       combined_genus, width = 14, height = 12)
cat("Saved: A2_genus_abundance_barplot.png\n")

# ==============================================================================
# SECTION A3: KRAKEN2 vs BRACKEN COMPARISON
# ==============================================================================

#get relative abundance table from a count df 
get_rel_abund <- function(count_df) {
  sample_cols <- rownames(metadata)
  mat         <- as.matrix(count_df[, sample_cols])
  rownames(mat) <- count_df$name
  # Convert to relative abundance
  sweep(mat, 2, colSums(mat), "/")
}

#SPECIES-LEVEL COMPARISON ------------------------------------------------------
rel_bracken_sp <- get_rel_abund(bracken_merged)
rel_kraken_sp  <- get_rel_abund(kraken_merged)

#Find species present in both
shared_sp <- intersect(rownames(rel_bracken_sp), rownames(rel_kraken_sp))
cat("Species in Bracken:", nrow(rel_bracken_sp), "\n")
cat("Species in Kraken: ", nrow(rel_kraken_sp),  "\n")
cat("Shared species:    ", length(shared_sp),     "\n")

#Build comparison dataframe (mean relative abundance across all samples)
sp_compare <- data.frame(
  Species       = shared_sp,
  Bracken_mean  = rowMeans(rel_bracken_sp[shared_sp, ]),
  Kraken_mean   = rowMeans(rel_kraken_sp[shared_sp, ])
) %>%
  mutate(Difference = Bracken_mean - Kraken_mean,
         Avg        = (Bracken_mean + Kraken_mean) / 2)

#Scatter plot: Kraken vs Bracken mean relative abundance
sp_scatter <- ggplot(sp_compare, aes(x = Kraken_mean, y = Bracken_mean)) +
  geom_point(alpha = 0.4, size = 1.5, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "red", linewidth = 0.7) +
  geom_text_repel(
    data  = sp_compare %>% arrange(desc(abs(Difference))) %>% slice_head(n = 15),
    aes(label = Species), size = 2.5, max.overlaps = 15, color = "grey30"
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  labs(title    = "Species-Level: Kraken2 vs Bracken",
       subtitle = "Each point = one species | dashed line = perfect agreement",
       x        = "Kraken2 Mean Relative Abundance",
       y        = "Bracken Mean Relative Abundance")

#Bland-Altman style plot: average vs difference
sp_ba <- ggplot(sp_compare, aes(x = Avg, y = Difference)) +
  geom_point(alpha = 0.4, size = 1.5, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_text_repel(
    data  = sp_compare %>% arrange(desc(abs(Difference))) %>% slice_head(n = 15),
    aes(label = Species), size = 2.5, max.overlaps = 15, color = "grey30"
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  labs(title    = "Species-Level: Bland-Altman (Bracken - Kraken2)",
       subtitle = "Points above 0 = higher in Bracken | Top 15 differences labelled",
       x        = "Mean Relative Abundance (average of both tools)",
       y        = "Difference (Bracken - Kraken2)")

combined_sp_compare <- sp_scatter + sp_ba +
  plot_annotation(
    title    = "Kraken2 vs Bracken: Species-Level Agreement",
    subtitle = "Mean relative abundance across all 6 samples",
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 10, color = "grey40"))
  )

print(combined_sp_compare)
ggsave(file.path(data_dir, "A3a_species_kraken_vs_bracken.png"),
       combined_sp_compare, width = 14, height = 6)
cat("Saved: A3a_species_kraken_vs_bracken.png\n")

#Pearson correlation
sp_cor <- cor.test(sp_compare$Kraken_mean, sp_compare$Bracken_mean,
                   method = "pearson")
cat(sprintf("\nSpecies-level Pearson r = %.4f, p = %.4e\n",
            sp_cor$estimate, sp_cor$p.value))

#GENUS-LEVEL COMPARISON --------------------------------------------------------

#genus-level relative abundance for Bracken
get_bracken_genus <- function() {
  sample_cols   <- rownames(metadata)
  tmp           <- bracken_merged
  tmp$Genus     <- sapply(strsplit(tmp$name, " "), `[`, 1)
  genus_df      <- tmp %>%
    dplyr::select(-name, -taxonomy_id) %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarise(across(all_of(sample_cols), sum), .groups = "drop") %>%
    as.data.frame()
  mat           <- as.matrix(genus_df[, sample_cols])
  rownames(mat) <- genus_df$Genus
  sweep(mat, 2, colSums(mat), "/")
}

#get Kraken genus-level relative abundance
get_kraken_genus <- function() {
  sample_cols <- rownames(metadata)
  genus_list  <- lapply(kraken_files, function(f) {
    srr         <- sub("\\.kraken\\.report$", "", basename(f))
    sample_name <- srr_to_sample[[srr]]
    df          <- read_tsv(f, col_types = cols(), show_col_types = FALSE,
                            col_names = c("pct","clade_reads","direct_reads",
                                          "rank","taxid","name"))
    df$name     <- trimws(df$name)
    df          <- df[df$rank == "G", ]
    out         <- data.frame(Genus = df$name, stringsAsFactors = FALSE)
    out[[sample_name]] <- df$direct_reads
    return(out)
  })
  genus_df    <- purrr::reduce(genus_list,
                               function(a, b) merge(a, b, by = "Genus", all = TRUE))
  genus_df[is.na(genus_df)] <- 0
  mat           <- as.matrix(genus_df[, sample_cols])
  rownames(mat) <- genus_df$Genus
  sweep(mat, 2, colSums(mat), "/")
}

rel_bracken_g <- get_bracken_genus()
rel_kraken_g  <- get_kraken_genus()

shared_g <- intersect(rownames(rel_bracken_g), rownames(rel_kraken_g))
cat("Genera in Bracken:", nrow(rel_bracken_g), "\n")
cat("Genera in Kraken: ", nrow(rel_kraken_g),  "\n")
cat("Shared genera:    ", length(shared_g),     "\n")

g_compare <- data.frame(
  Genus        = shared_g,
  Bracken_mean = rowMeans(rel_bracken_g[shared_g, ]),
  Kraken_mean  = rowMeans(rel_kraken_g[shared_g, ])
) %>%
  mutate(Difference = Bracken_mean - Kraken_mean,
         Avg        = (Bracken_mean + Kraken_mean) / 2)

#scatter
g_scatter <- ggplot(g_compare, aes(x = Kraken_mean, y = Bracken_mean)) +
  geom_point(alpha = 0.5, size = 2, color = "darkorchid") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "red", linewidth = 0.7) +
  geom_text_repel(
    data  = g_compare %>% arrange(desc(abs(Difference))) %>% slice_head(n = 15),
    aes(label = Genus), size = 2.8, max.overlaps = 15, color = "grey30"
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  labs(title    = "Genus-Level: Kraken2 vs Bracken",
       subtitle = "Each point = one genus | dashed line = perfect agreement",
       x        = "Kraken2 Mean Relative Abundance",
       y        = "Bracken Mean Relative Abundance")

#Bland-Altman
g_ba <- ggplot(g_compare, aes(x = Avg, y = Difference)) +
  geom_point(alpha = 0.5, size = 2, color = "darkorchid") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_text_repel(
    data  = g_compare %>% arrange(desc(abs(Difference))) %>% slice_head(n = 15),
    aes(label = Genus), size = 2.8, max.overlaps = 15, color = "grey30"
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  labs(title    = "Genus-Level: Bland-Altman (Bracken - Kraken2)",
       subtitle = "Points above 0 = higher in Bracken | Top 15 differences labelled",
       x        = "Mean Relative Abundance (average of both tools)",
       y        = "Difference (Bracken - Kraken2)")

combined_g_compare <- g_scatter + g_ba +
  plot_annotation(
    title    = "Kraken2 vs Bracken: Genus-Level Agreement",
    subtitle = "Mean relative abundance across all 6 samples",
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 10, color = "grey40"))
  )

print(combined_g_compare)
ggsave(file.path(data_dir, "A3b_genus_kraken_vs_bracken.png"),
       combined_g_compare, width = 14, height = 6)
cat("Saved: A3b_genus_kraken_vs_bracken.png\n")

#Pearson correlation
g_cor <- cor.test(g_compare$Kraken_mean, g_compare$Bracken_mean,
                  method = "pearson")
cat(sprintf("Genus-level Pearson r = %.4f, p = %.4e\n",
            g_cor$estimate, g_cor$p.value))
# ==============================================================================
# SECTION B: ALPHA DIVERSITY: Chao1, Shannon, Berger-Parker
# ==============================================================================

#Chao1 and Shannon via phyloseq
alpha_std <- estimate_richness(ps_bracken, measures = c("Chao1", "Shannon"))
alpha_std$SampleID <- rownames(alpha_std)

#Berger-Parker: max species reads / total reads per sample
otu_raw <- as.data.frame(otu_table(ps_bracken))  # taxa x samples
bp_vals <- apply(otu_raw, 2, function(x) max(x) / sum(x))
bp_df <- data.frame(SampleID = names(bp_vals), BergerParker = as.numeric(bp_vals))

#Combine
alpha_df <- alpha_std %>%
  left_join(bp_df, by = "SampleID") %>%
  left_join(rownames_to_column(metadata, "SampleID"), by = "SampleID")

cat("\n--- Alpha Diversity ---\n")
print(alpha_df[, c("SampleID","Diet","Chao1","Shannon","BergerParker")])

#Plot
alpha_long <- alpha_df %>%
  pivot_longer(cols      = c(Chao1, Shannon, BergerParker),
               names_to  = "Measure",
               values_to = "Value") %>%
  mutate(Measure = dplyr::recode(Measure,
                          "Chao1" = "Chao1\n(Richness)",
                          "Shannon" = "Shannon\n(Evenness)",
                          "BergerParker" = "Berger-Parker\n(Dominance)"
  ))

alpha_plot <- ggplot(alpha_long,
                     aes(x = Diet, y = Value, color = Diet, fill = Diet)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, linewidth = 0.7) +
  geom_jitter(width = 0.08, size = 3.5, alpha = 0.9) +
  scale_color_manual(values = diet_colors) +
  scale_fill_manual(values  = diet_colors) +
  facet_wrap(~Measure, scales = "free_y", nrow = 1) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 10),
    axis.title.x     = element_blank(),
    plot.title       = element_text(face = "bold", size = 13)
  ) +
  labs(title    = "Alpha Diversity: Omni vs Vegan",
       subtitle = "Bracken species-level counts | n=3 per group",
       y        = "Value")

print(alpha_plot)
ggsave(file.path(data_dir, "B_alpha_diversity.png"),
       alpha_plot, width = 11, height = 5)
cat("Saved: B_alpha_diversity.png\n")

# ==============================================================================
# SECTION C: BETA DIVERSITY: Bray-Curtis dissimilarity + PERMANOVA
# ==============================================================================

ps_rel  <- transform_sample_counts(ps_bracken, function(x) x / sum(x))
bc_dist <- phyloseq::distance(ps_rel, method = "bray")

cat("\nBray-Curtis Distance Matrix:\n")
print(round(as.matrix(bc_dist), 3))

#PCoA --------------------------------------------------------------------------
ord     <- ordinate(ps_rel, method = "PCoA", distance = "bray")
eig     <- ord$values$Eigenvalues
pct_var <- round(100 * eig / sum(eig[eig > 0]), 1)

pcoa_df <- data.frame(
  PCoA1    = ord$vectors[, 1],
  PCoA2    = ord$vectors[, 2],
  SampleID = rownames(ord$vectors)
) %>% left_join(rownames_to_column(metadata, "SampleID"), by = "SampleID")

pcoa_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2,
                                 color = Diet, label = SampleID)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(size = 3.2, show.legend = FALSE) +
  scale_color_manual(values = diet_colors) +
  theme_bw(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold")) +
  labs(title    = "Bray-Curtis PCoA",
       subtitle = "Omni vs Vegan | Bracken species-level",
       x        = paste0("PCoA1 (", pct_var[1], "% variance)"),
       y        = paste0("PCoA2 (", pct_var[2], "% variance)"),
       color    = "Diet")

print(pcoa_plot)
ggsave(file.path(data_dir, "C1_beta_diversity_pcoa.png"),
       pcoa_plot, width = 8, height = 6)
cat("Saved: C1_beta_diversity_pcoa.png\n")

#NMDS --------------------------------------------------------------------------
set.seed(42)
meta_df <- rownames_to_column(metadata, "SampleID")
ord_nmds <- ordinate(ps_rel, method = "NMDS", distance = "bray")
cat(sprintf("\nNMDS stress: %.4f", ord_nmds$stress),
    ifelse(ord_nmds$stress < 0.1, " (excellent)\n",
           ifelse(ord_nmds$stress < 0.2, " (good)\n", " (fair — interpret carefully)\n")))

nmds_df <- data.frame(
  NMDS1    = ord_nmds$points[, 1],
  NMDS2    = ord_nmds$points[, 2],
  SampleID = rownames(ord_nmds$points)
) %>% left_join(meta_df, by = "SampleID")

nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2,
                                 color = Diet, label = SampleID)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(size = 3.2, show.legend = FALSE) +
  scale_color_manual(values = diet_colors) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("Stress = ", round(ord_nmds$stress, 3)),
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "grey40") +
  theme_bw(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold")) +
  labs(title    = "Bray-Curtis NMDS",
       subtitle = "Omni vs Vegan | Bracken species-level",
       x = "NMDS1", y = "NMDS2",
       color = "Diet")

print(nmds_plot)
ggsave(file.path(data_dir, "C2_beta_diversity_nmds.png"),
       nmds_plot, width = 8, height = 6)
cat("Saved: C2_beta_diversity_nmds.png\n")

#PERMANOVA
samp_df   <- data.frame(sample_data(ps_rel))
permanova <- adonis2(bc_dist ~ Diet, data = samp_df, permutations = 999)
cat("\n--- PERMANOVA Results ---\n")
print(permanova)
write_csv(as.data.frame(permanova) %>% rownames_to_column("Term"),
          file.path(data_dir, "C_permanova_results.csv"))
cat("Saved: C_permanova_results.csv\n")

# ================================================================
# SECTION D: DIFFERENTIAL ABUNDANCE — ALDEx2
# Vegan vs Omni, species level (Bracken counts)
# ================================================================

# ALDEx2 needs integer counts
otu_int    <- round(as.data.frame(otu_table(ps_bracken)))
conditions <- as.character(sample_data(ps_bracken)$Diet)
cat("\nALDEx2 sample order and conditions:\n")
print(data.frame(Sample = colnames(otu_int), Diet = conditions))

# Run ALDEx2
set.seed(42)
cat("\nRunning ALDEx2 (128 Monte Carlo instances) — please wait...\n")
aldex_out <- aldex(
  reads      = otu_int,
  conditions = conditions,
  mc.samples = 128,
  test       = "t",
  effect     = TRUE,
  denom      = "iqlr"    # robust for sparse microbiome data
)

# Annotate with species names and sort
aldex_results <- aldex_out %>%
  rownames_to_column("taxonomy_id") %>%
  left_join(bracken_merged %>% select(taxonomy_id, name), by = "taxonomy_id") %>%
  arrange(wi.eBH)

cat("\n--- Top 20 Differentially Abundant Species (ALDEx2, Vegan vs Omni) ---\n")
print(head(aldex_results %>%
  select(name, we.ep, we.eBH, wi.ep, wi.eBH, effect), 20))

write_csv(aldex_results, file.path(data_dir, "D_aldex2_results.csv"))
cat("Saved: D_aldex2_results.csv\n")

# -- D1. MW plot (between-group vs within-group difference) --
mw_plot <- ggplot(aldex_results,
                  aes(x = diff.win, y = diff.btw, color = wi.eBH < 0.05)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#E63946"),
                     labels = c("Not significant", "BH-adjusted p < 0.05")) +
  geom_text_repel(
    data  = aldex_results %>% filter(wi.eBH < 0.05),
    aes(label = name), size = 3, max.overlaps = 15, color = "#E63946"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", size = 13),
        legend.title = element_text(face = "bold")) +
  labs(title    = "ALDEx2: Differential Abundance (Vegan vs Omni)",
       subtitle = "MW plot — positive diff.btw = higher in Vegan",
       x        = "Dispersion (within-group difference, CLR)",
       y        = "Difference between groups (CLR)",
       color    = "Significance")

print(mw_plot)
ggsave(file.path(data_dir, "D_aldex2_MW_plot.png"),
       mw_plot, width = 10, height = 7)
cat("Saved: D_aldex2_MW_plot.png\n")

# -- D2. Effect size bar chart (top 20 species by |effect|) --
top_effect <- aldex_results %>%
  arrange(desc(abs(effect))) %>%
  slice_head(n = 20) %>%
  mutate(name      = str_trunc(name, 40),
         Direction = ifelse(effect > 0, "Higher in Vegan", "Higher in Omni"))

effect_plot <- ggplot(top_effect,
                      aes(x = effect, y = reorder(name, effect),
                          fill = Direction)) +
  geom_col(color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.6) +
  scale_fill_manual(values = c("Higher in Vegan" = "#2A9D8F",
                               "Higher in Omni"  = "#E76F51")) +
  theme_bw(base_size = 11) +
  theme(axis.text.y  = element_text(size = 8),
        plot.title   = element_text(face = "bold", size = 13),
        legend.title = element_text(face = "bold")) +
  labs(title    = "ALDEx2: Top 20 Species by Effect Size",
       subtitle = "CLR effect size | Vegan vs Omni",
       x        = "Effect Size (CLR)",
       y        = NULL,
       fill     = "Direction")

print(effect_plot)
ggsave(file.path(data_dir, "D_aldex2_effect_size.png"),
       effect_plot, width = 10, height = 8)
cat("Saved: D_aldex2_effect_size.png\n")


# ================================================================
# SUMMARY
# ================================================================
cat("\n========================================\n")
cat("  Analysis Complete! Output files:\n")
cat("========================================\n")
cat("A1_abundance_barplot.png\n")
cat("A2_genus_abundance_barplot.png\n")
cat("A3a_species_kraken_vs_bracken.png\n")
cat("A3b_genus_kraken_vs_bracken.png\n")
cat("B_alpha_diversity.png\n")
cat("C1_beta_diversity_pcoa.png\n")
cat("C2_beta_diversity_nmds.png\n")
cat("C_permanova_results.csv\n")
cat("D_aldex2_results.csv\n")
cat("D_aldex2_effect_size.png\n")
cat("D_aldex2_MW_plot.png\n")
cat("All saved to:", data_dir, "\n")

