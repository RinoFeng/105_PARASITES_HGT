
# Revised plotting script for phylogeny-retained HGT candidates
# Input: revision/main_phylogeny_retained_HGT_8985.tsv
# Purpose: regenerate key manuscript/supplement figures using the
#          8,985 gene-tree curated HGT candidates.
#
# Notes on object definitions:
#   - all_filtered_en in the old workflow refers to the unintegrated dual-pipeline pool
#     (nucleotide-derived + protein-derived records) and is NOT recreated here.
#   - integrate_data in this script refers only to the 8,985 phylogeny-retained,
#     integrated HGT candidates used for the revised main analyses.
#   - retained_integrated_en is an explicit alias for integrate_data.
#   - data_type_overlap / Discovery_mode indicates the original discovery mode
#     (Only Protein, Overlap [Protein-ref], Only Nucleotide), not an evidence tier.

# --------------------------- -1 修改源文件----------------------------
library(readr)
library(dplyr)
library(stringr)


# 1. 文件路径


retained_file <- "revision/main_phylogeny_retained_HGT_8985.tsv"

# 如果你从项目根目录运行，用这个
overlap_file <- "../all_HGT.prot.one2one.blastx.iden60_cov80.acc"

# 如果你是在 revision/ 文件夹内部运行，则改成：
# overlap_file <- "../all_HGT.prot.one2one.blastx.iden60_cov80.acc"

out_file <- "revision/main_phylogeny_retained_HGT_8985.tsv"
backup_file <- "revision/main_phylogeny_retained_HGT_8985.before_reannotate_data_type_overlap.tsv"


# 2. 读取数据


hgt <- read_tsv(retained_file, show_col_types = FALSE)

overlap_prot <- read.table(
  overlap_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

overlap_acc <- overlap_prot[[1]] %>%
  as.character() %>%
  str_trim() %>%
  unique()


# 3. 基础检查


if (!"Accession" %in% colnames(hgt)) {
  stop("ERROR: retained HGT table does not contain column: Accession")
}

if (!"Data_type" %in% colnames(hgt)) {
  stop("ERROR: retained HGT table does not contain column: Data_type")
}

cat("Total retained HGTs:", nrow(hgt), "\n")
cat("Unique retained Accession:", n_distinct(hgt$Accession), "\n")
cat("Overlap protein-ref accession loaded:", length(overlap_acc), "\n")
cat("Overlap accessions found in retained table:", sum(hgt$Accession %in% overlap_acc), "\n\n")

cat("Original data_type_overlap distribution:\n")
if ("data_type_overlap" %in% colnames(hgt)) {
  print(table(hgt$data_type_overlap, useNA = "ifany"))
} else {
  cat("No existing data_type_overlap column found. It will be created.\n")
}


# 4. 备份原文件


write_tsv(hgt, backup_file)
cat("\nBackup written to:", backup_file, "\n")


# 5. 重新注释 data_type_overlap


hgt_fixed <- hgt %>%
  mutate(
    Accession = as.character(Accession),
    Data_type = as.character(Data_type),
    
    data_type_overlap = case_when(
      Accession %in% overlap_acc ~ "Overlap (Protein-ref)",
      
      str_detect(Data_type, regex("^Nucleotide$", ignore_case = TRUE)) ~ "Only Nucleotide",
      
      str_detect(Data_type, regex("^Protein$", ignore_case = TRUE)) ~ "Only Protein",
      
      TRUE ~ "Unclassified"
    ),
    
    data_type_overlap = factor(
      data_type_overlap,
      levels = c("Only Protein", "Overlap (Protein-ref)", "Only Nucleotide", "Unclassified")
    )
  )

cat("\nCorrected data_type_overlap distribution:\n")
print(table(hgt_fixed$data_type_overlap, useNA = "ifany"))

cat("\nCross-tabulation: Data_type × corrected data_type_overlap\n")
print(table(hgt_fixed$Data_type, hgt_fixed$data_type_overlap, useNA = "ifany"))


# 6. 检查是否有 Unclassified


unclassified_n <- sum(hgt_fixed$data_type_overlap == "Unclassified", na.rm = TRUE)

if (unclassified_n > 0) {
  warning(paste0("There are ", unclassified_n, " Unclassified records. Please check Data_type values."))
  print(
    hgt_fixed %>%
      filter(data_type_overlap == "Unclassified") %>%
      count(Data_type, sort = TRUE)
  )
}


# 7. 覆盖写回源文件


hgt_fixed <- hgt_fixed %>%
  mutate(data_type_overlap = as.character(data_type_overlap))

write_tsv(hgt_fixed, out_file)

cat("\nCorrected file written to:", out_file, "\n")
cat("Done.\n")

# --------------------------- 0. Packages ---------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(ggsci)
  library(scales)
  library(svglite)
  library(stringr)
  library(forcats)
  library(ggpubr)
  library(ggh4x)
})

# --------------------------- 1. Paths ------------------------------
# Run this script from the project root. The retained HGT table is expected
# under revision/ according to the revised manuscript workflow.
hgt_file <- "revision/main_phylogeny_retained_HGT_8985.tsv"
if (!file.exists(hgt_file)) {
  stop("Cannot find retained HGT table: ", hgt_file,
       "\nPlease put main_phylogeny_retained_HGT_8985.tsv under revision/ or edit hgt_file.")
}

outdir <- "revision/plots_retained_HGT"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --------------------------- 2. Themes -----------------------------
big_text_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  legend.position = "top",
  strip.text = element_text(face = "bold", size = 18),
  panel.grid = element_blank(),
  axis.text = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title = element_text(size = 20, face = "bold"),
  plot.subtitle = element_text(hjust = 0.5, size = 18, face = "bold"),
  legend.text = element_text(size = 15),
  legend.title = element_text(size = 16, face = "bold")
)

common_theme <- theme_minimal(base_size = 20) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

clade_colors <- c(
  "Microbial eukaryotes" = "#2E86AB",
  "Metazoan parasites" = "#F18F01"
)

life_colors <- c(
  "Intracellular" = "#DC0000B2",
  "Extracellular" = "#3C5488B2",
  "Free-living" = "#00A087B2"
)

discovery_colors <- c(
  "Only Protein" = "#2E5077",
  "Overlap (Protein-ref)" = "#F18F01",
  "Only Nucleotide" = "#4DA1A9"
)

hgt_category_colors <- c(
  "Donor-placed HGT candidates" = "#8B1E3F",
  "Donor-leaning HGT candidates" = "#F18F01",
  "Lineage-retained ancestral HGT candidates" = "#2E86AB"
)

# --------------------------- 3. Load and clean ----------------------
integrate_data <- readr::read_tsv(hgt_file, show_col_types = FALSE) %>%
  mutate(
    manuscript_HGT_category = factor(
      manuscript_HGT_category,
      levels = c(
        "Donor-placed HGT candidates",
        "Donor-leaning HGT candidates",
        "Lineage-retained ancestral HGT candidates"
      )
    ),
    manuscript_HGT_category_short = case_when(
      manuscript_HGT_category == "Donor-placed HGT candidates" ~ "Donor-placed",
      manuscript_HGT_category == "Donor-leaning HGT candidates" ~ "Donor-leaning",
      manuscript_HGT_category == "Lineage-retained ancestral HGT candidates" ~ "Lineage-retained ancestral",
      TRUE ~ as.character(manuscript_HGT_category)
    ),
    manuscript_HGT_category_short = factor(
      manuscript_HGT_category_short,
      levels = c("Donor-placed", "Donor-leaning", "Lineage-retained ancestral")
    ),
    Discovery_mode = case_when(
      data_type_overlap %in% c("Only Protein", "Protein-only", "Protein only") ~ "Only Protein",
      data_type_overlap %in% c("Only Nucleotide", "Nucleotide-only", "Nucleotide only") ~ "Only Nucleotide",
      str_detect(as.character(data_type_overlap), "Overlap") ~ "Overlap (Protein-ref)",
      TRUE ~ as.character(data_type_overlap)
    ),
    Discovery_mode = factor(
      Discovery_mode,
      levels = c("Only Protein", "Overlap (Protein-ref)", "Only Nucleotide")
    ),
    # Keep the old column name as a plotting alias so downstream code remains unchanged.
    data_type_overlap = Discovery_mode,
    Data_type = factor(Data_type, levels = c("Protein", "Nucleotide")),
    Parasite_clade = case_when(
      Parasite_clade %in% c("Protozoa", "Microbial eukaryotes") ~ "Microbial eukaryotes",
      Parasite_clade %in% c("Metazoa", "Metazoan parasites") ~ "Metazoan parasites",
      TRUE ~ as.character(Parasite_clade)
    ),
    Parasite_clade = factor(Parasite_clade, levels = c("Microbial eukaryotes", "Metazoan parasites")),
    Parasitism_type = factor(Parasitism_type, levels = c("Intracellular", "Extracellular", "Free-living")),
    Length = as.numeric(Length),
    HGT_coverage = as.numeric(HGT_coverage),
    HGT_identity = as.numeric(HGT_identity),
    GC_Content = as.numeric(GC_Content),
    Contig_Length = as.numeric(Contig_Length),
    dist_to_edge = as.numeric(dist_to_edge),
    hgt_contig_fraction = as.numeric(hgt_contig_fraction)
  )

# Explicit retained integrated object for revised analyses.
retained_integrated_en <- integrate_data

# Do NOT overwrite old all_filtered_en here. In the original workflow, all_filtered_en
# represented the unintegrated nucleotide-derived + protein-derived candidate pool.
# If you need to compare against the old discovery pool, load it separately as
# old_all_filtered_en or initial_all_filtered_en.

# Main subsets used repeatedly.
donor_side_data <- integrate_data %>%
  filter(manuscript_HGT_category %in% c("Donor-placed HGT candidates", "Donor-leaning HGT candidates"))

strict_donor_placed_data <- integrate_data %>%
  filter(manuscript_HGT_category == "Donor-placed HGT candidates")

# --------------------------- 4. Helper functions -------------------
parse_donor_taxonomy <- function(df) {
  df %>%
    mutate(
      Donor_Kingdom = case_when(
        str_detect(Target_lineage, "Bacteria") ~ "Bacteria",
        str_detect(Target_lineage, "Viruses") ~ "Viruses",
        str_detect(Target_lineage, "Fungi") ~ "Fungi",
        str_detect(Target_lineage, "Viridiplantae|Streptophyta|Chlorophyta") ~ "Plants",
        str_detect(Target_lineage, "Metazoa") ~ "Metazoa",
        str_detect(Target_lineage, "Archaea") ~ "Archaea",
        str_detect(Target_lineage, "Eukaryota") ~ "Other Eukaryotes",
        TRUE ~ "Unknown"
      ),
      lineage_split = str_split(Target_lineage, ";"),
      Donor_Phylum = map_chr(lineage_split, function(x) {
        x <- trimws(x)
        x <- x[!x %in% c("", "organisms", "root", "cellular organisms", "Eukaryota", "Bacteria", "Viruses", "Archaea")]
        x <- x[!str_detect(x, "group$|environmental samples|unclassified")]
        if (length(x) == 0) return("Unclassified")
        out <- x[1]
        ifelse(nchar(out) > 30, "Others", out)
      })
    ) %>%
    select(-lineage_split)
}

save_plot <- function(plot, filename, width, height) {
  ggsave(file.path(outdir, paste0(filename, ".svg")), plot, width = width, height = height, device = svglite, fix_text_size = F)
}

# COG macro group definition retained from the previous script.
cog_groups <- list(
  Unknown = c("S"),
  `Transport and Metabolism` = c("E", "Q", "G", "I", "C", "P", "H", "F"),
  `Cellular Processes and Signalling` = c("O", "M", "T", "U", "V", "D", "N", "Y", "Z"),
  `Genetic Information Processing` = c("J", "L", "K", "A", "B")
)

add_cog_group <- function(df) {
  df %>%
    mutate(
      COG_category = ifelse(is.na(COG_category) | COG_category == "-", NA_character_, COG_category),
      COG_Group = case_when(
        COG_category %in% cog_groups[["Unknown"]] ~ "Unknown",
        COG_category %in% cog_groups[["Transport and Metabolism"]] ~ "Transport and Metabolism",
        COG_category %in% cog_groups[["Cellular Processes and Signalling"]] ~ "Cellular Processes and Signalling",
        COG_category %in% cog_groups[["Genetic Information Processing"]] ~ "Genetic Information Processing",
        TRUE ~ NA_character_
      ),
      COG_Group = factor(
        COG_Group,
        levels = c("Unknown", "Transport and Metabolism", "Cellular Processes and Signalling", "Genetic Information Processing")
      )
    )
}

# --------------------------- 5. Summary tables ----------------------
main_category_summary <- integrate_data %>%
  count(manuscript_HGT_category, name = "Count") %>%
  mutate(Percent = round(Count / sum(Count) * 100, 2))
write_tsv(main_category_summary, file.path(outdir, "Table_main_HGT_category_summary.tsv"))

source_by_category <- integrate_data %>%
  count(manuscript_HGT_category, data_type_overlap, name = "Count") %>%
  group_by(manuscript_HGT_category) %>%
  mutate(Percent_in_category = round(Count / sum(Count) * 100, 2)) %>%
  ungroup()
write_tsv(source_by_category, file.path(outdir, "Table_discovery_mode_by_HGT_category.tsv"))

contig_by_category <- integrate_data %>%
  count(manuscript_HGT_category, contig_context_tier, name = "Count") %>%
  group_by(manuscript_HGT_category) %>%
  mutate(Percent_in_category = round(Count / sum(Count) * 100, 2)) %>%
  ungroup()
write_tsv(contig_by_category, file.path(outdir, "Table_contig_tier_by_HGT_category.tsv"))

species_by_category <- integrate_data %>%
  count(Specie, Parasite_clade, manuscript_HGT_category, name = "Count") %>%
  pivot_wider(names_from = manuscript_HGT_category, values_from = Count, values_fill = 0)
write_tsv(species_by_category, file.path(outdir, "Table_species_by_main_HGT_category.tsv"))

# --------------------------- 6. Figure 1: retained HGT categories ----
p_category <- main_category_summary %>%
  mutate(Category = fct_reorder(manuscript_HGT_category, Count)) %>%
  ggplot(aes(x = Category, y = Count, fill = manuscript_HGT_category)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_text(aes(label = comma(Count)), hjust = -0.1, size = 5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = comma) +
  scale_fill_manual(values = hgt_category_colors) +
  labs(x = NULL, y = "Number of HGT candidates", title = "Phylogeny-retained HGT candidates") +
  theme_bw() +
  big_text_theme +
  theme(legend.position = "none", axis.text.y = element_text(size = 14, face = "bold"))

p_category
save_plot(p_category, "Figure1A_main_HGT_category_summary", 9, 4.5)

# Discovery mode composition within the three retained HGT patterns.
p_discovery_by_category <- source_by_category %>%
  ggplot(aes(x = manuscript_HGT_category, y = Count, fill = data_type_overlap)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = discovery_colors, name = "Discovery mode") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Number of HGT candidates", title = "Discovery modes within retained HGT categories") +
  theme_bw() +
  big_text_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 13))

p_discovery_by_category
save_plot(p_discovery_by_category, "Figure1B_discovery_mode_by_HGT_category", 11, 6)

# --------------------------- 7. Figure 2: HGT abundance and features --
# Quantity comparison: revised to use the 8,985 retained set, not protein/nucleotide as final groups.
hgt_metric_clade <- integrate_data %>%
  filter(!is.na(Parasite_clade)) %>%
  group_by(Parasite_clade) %>%
  summarise(
    total_hgt = n(),
    species_count = n_distinct(Specie),
    avg_hgt_per_species = round(total_hgt / species_count, 1),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(total_hgt, avg_hgt_per_species), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         total_hgt = "Total HGT candidates",
                         avg_hgt_per_species = "Avg. HGT per species"))
write_tsv(hgt_metric_clade, file.path(outdir, "Table_HGT_abundance_by_clade.tsv"))

p_abundance_clade <- ggplot(hgt_metric_clade, aes(x = Parasite_clade, y = Value, fill = Parasite_clade)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(aes(label = comma(Value)), vjust = -0.3, size = 5) +
  facet_wrap(vars(Metric), scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = comma) +
  scale_fill_manual(values = clade_colors) +
  labs(x = "Parasite group", y = "Number of HGT candidates", title = "HGT abundance after phylogenetic curation") +
  theme_bw() +
  big_text_theme +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 16))
p_abundance_clade

# Feature distributions by clade.
numeric_data_clade <- integrate_data %>%
  filter(!is.na(Parasite_clade)) %>%
  select(Parasite_clade, Length, HGT_coverage, HGT_identity, GC_Content) %>%
  pivot_longer(cols = c(Length, HGT_coverage, HGT_identity, GC_Content), names_to = "Indicator", values_to = "Value") %>%
  filter(!is.na(Value), is.finite(Value)) %>%
  mutate(Indicator = factor(
    Indicator,
    levels = c("Length", "HGT_coverage", "HGT_identity", "GC_Content"),
    labels = c("Sequence length", "HGT coverage", "HGT identity", "GC content (%)")
  ))

p_features_clade <- ggplot(numeric_data_clade, aes(x = Parasite_clade, y = Value, fill = Parasite_clade)) +
  geom_boxplot(alpha = 0.7, linewidth = 0.3, outlier.size = 0.5, outlier.alpha = 0.25) +
  facet_wrap(vars(Indicator), scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = FALSE, size = 5, label.x = 1.5) +
  scale_fill_manual(values = clade_colors) +
  labs(x = "Parasite group", y = "Value", title = "Sequence features of retained HGT candidates") +
  theme_bw() +
  big_text_theme +
  theme(legend.position = "none", axis.text.x = element_text(size = 13, color = "black"))
p_features_clade

final_features_clade <- p_abundance_clade / p_features_clade + plot_layout(heights = c(1, 1.4))
final_features_clade
save_plot(final_features_clade, "Figure2_HGT_abundance_and_features_by_clade", 13, 9)

# --------------------------- 8. Figure S: lifestyle features --------
hgt_metric_life <- integrate_data %>%
  filter(!is.na(Parasitism_type)) %>%
  group_by(Parasitism_type) %>%
  summarise(
    total_hgt = n(),
    species_count = n_distinct(Specie),
    avg_hgt_per_species = round(total_hgt / species_count, 1),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(total_hgt, avg_hgt_per_species), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         total_hgt = "Total HGT candidates",
                         avg_hgt_per_species = "Avg. HGT per species"))
write_tsv(hgt_metric_life, file.path(outdir, "Table_HGT_abundance_by_lifestyle.tsv"))

p_abundance_life <- ggplot(hgt_metric_life, aes(x = Parasitism_type, y = Value, fill = Parasitism_type)) +
  geom_col(alpha = 0.9, width = 0.65) +
  geom_text(aes(label = comma(Value)), vjust = -0.3, size = 5) +
  facet_wrap(vars(Metric), scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = comma) +
  scale_fill_manual(values = life_colors) +
  labs(x = "Lifestyle", y = "Number of HGT candidates", title = "HGT abundance by lifestyle") +
  theme_bw() +
  big_text_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
p_abundance_life

numeric_data_life <- integrate_data %>%
  filter(!is.na(Parasitism_type)) %>%
  select(Parasitism_type, Length, HGT_coverage, HGT_identity, GC_Content) %>%
  pivot_longer(cols = c(Length, HGT_coverage, HGT_identity, GC_Content), names_to = "Indicator", values_to = "Value") %>%
  filter(!is.na(Value), is.finite(Value)) %>%
  mutate(Indicator = factor(
    Indicator,
    levels = c("Length", "HGT_coverage", "HGT_identity", "GC_Content"),
    labels = c("Sequence length", "HGT coverage", "HGT identity", "GC content (%)")
  ))

life_groups <- levels(droplevels(integrate_data$Parasitism_type))
my_comparisons <- combn(life_groups[!is.na(life_groups)], 2, simplify = FALSE)

p_features_life <- ggplot(numeric_data_life, aes(x = Parasitism_type, y = Value, fill = Parasitism_type)) +
  geom_boxplot(alpha = 0.7, linewidth = 0.3, outlier.size = 0.5, outlier.alpha = 0.25) +
  facet_wrap(vars(Indicator), scales = "free_y", ncol = 4) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", hide.ns = FALSE, size = 4) +
  scale_fill_manual(values = life_colors) +
  labs(x = "Lifestyle", y = "Value", title = "Sequence features by lifestyle") +
  theme_bw() +
  big_text_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1, size = 13))
p_features_life

final_lifestyle <- p_abundance_life / p_features_life + plot_layout(heights = c(1, 1.4))
final_lifestyle
save_plot(final_lifestyle, "FigureS_HGT_abundance_and_features_by_lifestyle", 14, 10)

# --------------------------- 9. Figure: donor landscape by species ----
plot_hgt_with_anno <- function(data, title_text = "") {
  plot_data <- parse_donor_taxonomy(data)
  species_order <- plot_data %>%
    count(Specie, name = "n") %>%
    arrange(n) %>%
    pull(Specie)
  plot_data <- plot_data %>% mutate(Specie = factor(Specie, levels = species_order))
  species_anno <- plot_data %>%
    select(Specie, Parasitism_type, Parasite_clade) %>%
    distinct(Specie, .keep_all = TRUE) %>%
    tidyr::drop_na(Specie)
  bar_df <- plot_data %>% count(Specie, Donor_Kingdom, name = "n")
  p_main <- ggplot(bar_df, aes(x = Specie, y = n, fill = Donor_Kingdom)) +
    geom_col(width = 0.8) +
    coord_flip() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = comma) +
    scale_fill_npg(name = "Donor source") +
    labs(x = "Parasite species", y = "Number of retained HGT candidates", title = title_text) +
    theme_bw() +
    big_text_theme +
    theme(axis.text.y = element_text(face = "italic", size = 9, color = "black"),
          axis.text.x = element_text(size = 11, color = "black"),
          panel.grid.major.y = element_blank())
  p_anno1 <- ggplot(species_anno, aes(x = "Lifestyle", y = Specie, fill = Parasitism_type)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_manual(values = life_colors, name = "Lifestyle") +
    scale_y_discrete(drop = FALSE) +
    theme_minimal() +
    big_text_theme +
    theme(axis.title = element_blank(), axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), panel.grid = element_blank())
  p_anno2 <- ggplot(species_anno, aes(x = "Group", y = Specie, fill = Parasite_clade)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_manual(values = clade_colors, name = "Parasite group") +
    scale_y_discrete(drop = FALSE) +
    theme_minimal() +
    big_text_theme +
    theme(axis.title = element_blank(), axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), panel.grid = element_blank())
  p_main + p_anno1 + p_anno2 +
    plot_layout(widths = c(8, 0.45, 0.45), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "vertical")
}

fig_landscape_all <- plot_hgt_with_anno(integrate_data, "Retained HGT candidates across parasite species")
print(fig_landscape_all)
save_plot(fig_landscape_all, "FigureS_species_landscape_retained_HGT", 11, 17)

fig_landscape_donor_side <- plot_hgt_with_anno(donor_side_data, "Donor-side retained HGT candidates across parasite species")
print(fig_landscape_donor_side)
save_plot(fig_landscape_donor_side, "FigureS_species_landscape_donor_side_HGT", 11, 17)

# --------------------------- 10. Donor composition plots ------------
# Bar plot by lifestyle and donor kingdom. The main donor-source interpretation should
# primarily use donor-side candidates; all retained candidates are plotted as supplement.
plot_donor_by_group <- function(data, group_col, title_text, proportional = TRUE) {
  plot_df <- parse_donor_taxonomy(data) %>%
    filter(!is.na(.data[[group_col]])) %>%
    count(.data[[group_col]], Donor_Kingdom, name = "Count")
  if (proportional) {
    plot_df <- plot_df %>% group_by(.data[[group_col]]) %>% mutate(Value = Count / sum(Count)) %>% ungroup()
    y_lab <- "Proportion of HGT candidates"
    scale_y <- scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.02)))
  } else {
    plot_df <- plot_df %>% mutate(Value = Count)
    y_lab <- "Number of HGT candidates"
    scale_y <- scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05)))
  }
  ggplot(plot_df, aes(x = .data[[group_col]], y = Value, fill = Donor_Kingdom)) +
    geom_col(width = 0.75, color = "black", linewidth = 0.15) +
    scale_fill_npg(name = "Donor source") +
    scale_y +
    labs(x = NULL, y = y_lab, title = title_text) +
    theme_bw() +
    big_text_theme +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

p_donor_life_donor_side <- plot_donor_by_group(donor_side_data, "Parasitism_type", "Donor composition by lifestyle: donor-side HGTs", TRUE)
p_donor_clade_donor_side <- plot_donor_by_group(donor_side_data, "Parasite_clade", "Donor composition by parasite group: donor-side HGTs", TRUE)
final_donor_comp <- p_donor_life_donor_side / p_donor_clade_donor_side + plot_layout(heights = c(1, 1))
final_donor_comp
save_plot(final_donor_comp, "Figure3_donor_composition_donor_side_HGT", 11, 10)

p_donor_life_all <- plot_donor_by_group(integrate_data, "Parasitism_type", "Donor affinity by lifestyle: all retained HGTs", TRUE)
p_donor_clade_all <- plot_donor_by_group(integrate_data, "Parasite_clade", "Donor affinity by parasite group: all retained HGTs", TRUE)
p_donor_life_all / p_donor_clade_all 
save_plot(p_donor_life_all / p_donor_clade_all, "FigureS_donor_composition_all_retained_HGT", 11, 10)

# Sunburst-style donor plots, retained from previous figure style.
prepare_sunburst_data <- function(hgt_df) {
  hgt_merged <- parse_donor_taxonomy(hgt_df)
  total_n <- nrow(hgt_merged)
  df_phylum <- hgt_merged %>%
    count(Donor_Kingdom, Donor_Phylum, name = "n") %>%
    group_by(Donor_Kingdom) %>%
    arrange(desc(n), .by_group = TRUE) %>%
    mutate(Phylum_Grouped = ifelse(n / total_n < 0.005 | row_number() > 6,
                                   paste0("Other ", Donor_Kingdom), Donor_Phylum)) %>%
    group_by(Donor_Kingdom, Phylum_Grouped) %>%
    summarise(Count = sum(n), .groups = "drop") %>%
    rename(Phylum = Phylum_Grouped)
  inner_ring <- df_phylum %>%
    group_by(Donor_Kingdom) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    mutate(Fraction = Count / sum(Count), ymax = cumsum(Fraction), ymin = c(0, head(ymax, -1)),
           Label = paste0(Donor_Kingdom, "\n(", percent(Fraction, 0.1), ")"))
  outer_ring <- df_phylum %>%
    mutate(Global_Fraction = Count / sum(Count), ymax = cumsum(Global_Fraction), ymin = c(0, head(ymax, -1)),
           Label = ifelse(Global_Fraction > 0.01, paste0(Phylum, "\n(", Count, ")"), ""))
  list(inner = inner_ring, outer = outer_ring)
}

draw_sunburst <- function(data_list, title_text) {
  ggplot() +
    geom_rect(data = data_list$inner, aes(xmin = 0, xmax = 2, ymin = ymin, ymax = ymax, fill = Donor_Kingdom), color = "white", linewidth = 1) +
    geom_rect(data = data_list$outer, aes(xmin = 2.05, xmax = 3.5, ymin = ymin, ymax = ymax, fill = Donor_Kingdom), color = "white", alpha = 0.6, linewidth = 0.5) +
    geom_text(data = data_list$inner %>% filter(Fraction > 0.02), aes(x = 1, y = (ymin + ymax) / 2, label = Label), color = "white", fontface = "bold", size = 3.2) +
    geom_text(data = data_list$outer %>% filter(Label != ""), aes(x = 2.8, y = (ymin + ymax) / 2, label = Label), color = "black", size = 2.3, check_overlap = FALSE) +
    scale_fill_npg() +
    coord_polar(theta = "y") +
    xlim(c(0, 3.5)) +
    labs(title = title_text) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), legend.position = "none")
}

p_sun_all <- draw_sunburst(prepare_sunburst_data(integrate_data), "All retained HGT candidates")
p_sun_all
p_sun_donor <- draw_sunburst(prepare_sunburst_data(donor_side_data), "Donor-side retained HGT candidates")
p_sun_donor
save_plot(p_sun_all + p_sun_donor, "FigureS_donor_sunburst_retained_HGT", 14, 7)

# --------------------------- 11. COG functional plots ----------------
processed_cog_data <- integrate_data %>%
  add_cog_group() %>%
  filter(!is.na(COG_Group), str_detect(COG_category, "^[A-Z]$"))

cog_global_order <- processed_cog_data %>%
  count(COG_Group, COG_category, name = "Global_Count") %>%
  arrange(desc(Global_Count))

# A. Function distribution stratified by discovery mode.
stats_all <- processed_cog_data %>%
  count(COG_Group, COG_category, data_type_overlap, name = "Count") %>%
  mutate(Group = "All retained")
stats_lifestyle <- processed_cog_data %>%
  count(COG_Group, COG_category, data_type_overlap, Parasitism_type, name = "Count") %>%
  rename(Group = Parasitism_type)

plot_data_cog <- bind_rows(stats_all, stats_lifestyle) %>%
  left_join(cog_global_order, by = c("COG_Group", "COG_category")) %>%
  mutate(
    COG_category = fct_reorder(COG_category, Global_Count, .desc = TRUE),
    Group = fct_relevel(as.factor(Group), "All retained", "Intracellular", "Extracellular", "Free-living"),
    data_type_overlap = fct_relevel(data_type_overlap, "Only Protein", "Overlap (Protein-ref)", "Only Nucleotide")
  )

p_cog_top <- plot_data_cog %>%
  filter(Group == "All retained") %>%
  ggplot(aes(x = COG_category, y = Count, fill = data_type_overlap)) +
  geom_col(color = "black", linewidth = 0.15, width = 0.8) +
  scale_fill_manual(values = discovery_colors, name = "Discovery mode") +
  facet_grid(~ COG_Group, scales = "free_x", space = "free_x", switch = "x", labeller = labeller(COG_Group = label_wrap_gen(width = 15))) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "No. of HGT candidates", title = "Functional distribution of retained HGT candidates") +
  theme_bw() +
  big_text_theme +
  theme(axis.ticks.x = element_blank(), strip.placement = "outside", strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 10), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

p_cog_bottom <- plot_data_cog %>%
  filter(Group != "All retained") %>%
  ggplot(aes(x = COG_category, y = Count, fill = data_type_overlap)) +
  geom_col(color = "black", linewidth = 0.15, width = 0.8) +
  scale_fill_manual(values = discovery_colors, name = "Discovery mode") +
  facet_grid(Group ~ COG_Group, scales = "free", space = "free_x", labeller = labeller(COG_Group = label_wrap_gen(width = 15))) +
  scale_y_continuous(labels = comma) +
  labs(x = "COG categories", y = "No. of HGT candidates", title = "Functional distribution by lifestyle") +
  theme_bw() +
  big_text_theme +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.spacing.y = unit(0.6, "cm"))

final_cog_discovery <- p_cog_top / p_cog_bottom + plot_layout(heights = c(1, 2))
print(final_cog_discovery)
save_plot(final_cog_discovery, "Figure4_COG_lifestyle_by_discovery_mode_retained_HGT", 16, 12)

# B. Function distribution stratified by the three HGT retention patterns.
p_cog_category <- processed_cog_data %>%
  count(COG_Group, COG_category, manuscript_HGT_category, name = "Count") %>%
  left_join(cog_global_order, by = c("COG_Group", "COG_category")) %>%
  mutate(COG_category = fct_reorder(COG_category, Global_Count, .desc = TRUE)) %>%
  ggplot(aes(x = COG_category, y = Count, fill = manuscript_HGT_category)) +
  geom_col(color = "black", linewidth = 0.15, width = 0.8) +
  scale_fill_manual(values = hgt_category_colors, name = "HGT pattern") +
  facet_grid(~ COG_Group, scales = "free_x", space = "free_x", switch = "x", labeller = labeller(COG_Group = label_wrap_gen(width = 15))) +
  scale_y_continuous(labels = comma) +
  labs(x = "COG categories", y = "No. of HGT candidates", title = "COG distribution by HGT retention pattern") +
  theme_bw() +
  big_text_theme +
  theme(axis.text.x = element_text(size = 13, face = "bold"), strip.placement = "outside", strip.background = element_rect(fill = "grey95", color = NA),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
print(p_cog_category)
save_plot(p_cog_category, "FigureS_COG_by_HGT_retention_pattern", 15, 8)

# C. Heatmap: species-resolved COG category counts using retained HGTs.
# C. ComplexHeatmap: species-resolved COG functional profiles using retained HGTs
# This replaces the simplified ggplot geom_tile heatmap.

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

library(RColorBrewer)

if (!require("ComplexHeatmap", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap", update = FALSE)
}
library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)


# 1. Read taxonomy information

tax_raw <- read.delim(
  "../105_parasite_taxonomy_info.txt",
  header = FALSE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE
)

get_phylum_group <- function(lineage) {
  if (grepl("Apicomplexa", lineage)) return("Apicomplexa")
  if (grepl("Kinetoplastea|Kinetoplastida", lineage)) return("Kinetoplastida")
  if (grepl("Nematoda", lineage)) return("Nematoda")
  if (grepl("Platyhelminthes", lineage)) return("Platyhelminthes")
  if (grepl("Arthropoda", lineage)) return("Arthropoda")
  return("Other microbial eukaryotes")
}

tax_map <- tax_raw %>%
  mutate(
    Phylum_Group = sapply(V2, get_phylum_group),
    Species_Name_Full = str_trim(sub(".*;", "", V3)),
    Genus = word(Species_Name_Full, 1)
  ) %>%
  select(Genus, Phylum_Group) %>%
  distinct()


# 2. Prepare retained HGT COG data

if (!exists("processed_cog_data")) {
  stop("请先生成 processed_cog_data。该对象应来自 8,985 条 retained HGT，并包含 Specie/Name、COG_category、COG_Group 等列。")
}

# 兼容旧表中的 Name 和新版表中的 Specie
clean_data <- processed_cog_data %>%
  mutate(
    Species_Name = case_when(
      "Specie" %in% colnames(.) ~ as.character(Specie),
      TRUE ~ NA_character_
    ),
    Species_Name = ifelse(is.na(Species_Name), "Unknown_species", Species_Name),
    Genus = word(Species_Name, 1)
  ) %>%
  left_join(tax_map, by = "Genus") %>%
  mutate(
    Phylum_Group = case_when(
      !is.na(Phylum_Group) ~ Phylum_Group,
      "Parasite_clade" %in% colnames(.) & Parasite_clade == "Metazoa" ~ "Metazoan parasites",
      TRUE ~ "Other microbial eukaryotes"
    )
  ) %>%
  filter(!is.na(COG_category), COG_category != "") %>%
  filter(str_detect(COG_category, "^[A-Z]$"))

# 如果你想去掉 Unknown COG S，可以保留这一行；
# 如果想展示 S 类 unknown function，就注释掉这一行。
clean_data <- clean_data %>%
  filter(COG_category != "S")


# 3. Count COG categories and normalize within each species

tog_counts <- clean_data %>%
  group_by(Species_Name, Phylum_Group, Parasitism_type, COG_category, COG_Group) %>%
  summarise(count = n(), .groups = "drop")

tog_norm <- tog_counts %>%
  complete(
    nesting(Species_Name, Phylum_Group, Parasitism_type),
    nesting(COG_category, COG_Group),
    fill = list(count = 0)
  ) %>%
  group_by(Species_Name) %>%
  mutate(
    species_total = sum(count),
    percent = ifelse(species_total > 0, count / species_total * 100, 0)
  ) %>%
  ungroup()


# 4. Build heatmap matrix
# rows = COG categories; columns = species
# values = within-species percentage

mat_raw <- tog_norm %>%
  select(Species_Name, COG_category, percent) %>%
  pivot_wider(names_from = Species_Name, values_from = percent) %>%
  column_to_rownames("COG_category") %>%
  as.matrix()

final_matrix <- mat_raw
final_matrix[is.na(final_matrix)] <- 0

# Remove empty rows or columns if any
final_matrix <- final_matrix[rowSums(final_matrix) > 0, colSums(final_matrix) > 0, drop = FALSE]


# 5. Column and row annotations

species_meta <- clean_data %>%
  select(Species_Name, Phylum_Group, Parasitism_type) %>%
  distinct(Species_Name, .keep_all = TRUE)

mat_col_annotation <- data.frame(Species_Name = colnames(final_matrix)) %>%
  left_join(species_meta, by = "Species_Name") %>%
  mutate(
    Phylum_Group = ifelse(is.na(Phylum_Group), "Unknown", Phylum_Group),
    Parasitism_type = ifelse(is.na(Parasitism_type), "Unknown", Parasitism_type)
  ) %>%
  as.data.frame()

cog_meta <- clean_data %>%
  select(COG_category, COG_Group) %>%
  distinct()

mat_row_annotation <- data.frame(COG_category = rownames(final_matrix)) %>%
  left_join(cog_meta, by = "COG_category") %>%
  mutate(
    COG_Group = ifelse(is.na(COG_Group), "Unknown", as.character(COG_Group))
  ) %>%
  as.data.frame()

# Species total retained HGT counts for bottom barplot
total_hgt_counts <- tog_counts %>%
  group_by(Species_Name) %>%
  summarise(Total = sum(count), .groups = "drop")

mat_bottom_data <- data.frame(Species_Name = colnames(final_matrix)) %>%
  left_join(total_hgt_counts, by = "Species_Name") %>%
  mutate(Total = ifelse(is.na(Total), 0, Total))


# 6. Color settings

positive_values <- final_matrix[final_matrix > 0]
if (length(positive_values) == 0) {
  max_limit <- 1
} else {
  max_limit <- as.numeric(quantile(positive_values, 0.98, na.rm = TRUE))
}

col_heatmap <- colorRamp2(
  c(0, max_limit / 4, max_limit),
  c("white", "#FD8D3C", "#800026")
)

func_levels <- unique(mat_row_annotation$COG_Group)
func_colors <- structure(
  colorRampPalette(brewer.pal(min(12, max(3, length(func_levels))), "Paired"))(length(func_levels)),
  names = func_levels
)

phylum_colors <- c(
  "Apicomplexa" = "#E64B35",
  "Kinetoplastida" = "#4DBBD5",
  "Nematoda" = "#00A087",
  "Platyhelminthes" = "#3C5488",
  "Arthropoda" = "#F39B7F",
  "Other microbial eukaryotes" = "#8491B4",
  "Metazoan parasites" = "#7E6148",
  "Unknown" = "grey80"
)

# Add colors dynamically for missing phylum groups
missing_phylum <- setdiff(unique(mat_col_annotation$Phylum_Group), names(phylum_colors))
if (length(missing_phylum) > 0) {
  phylum_colors <- c(
    phylum_colors,
    structure(rep("grey70", length(missing_phylum)), names = missing_phylum)
  )
}

parasitism_levels <- unique(mat_col_annotation$Parasitism_type)
parasitism_colors <- structure(
  brewer.pal(max(3, length(parasitism_levels)), "Set2")[seq_along(parasitism_levels)],
  names = parasitism_levels
)
if ("Unknown" %in% names(parasitism_colors)) {
  parasitism_colors["Unknown"] <- "grey80"
}

all_colors <- list(
  Phylum = phylum_colors,
  Parasitism = parasitism_colors,
  Func_Class = func_colors
)


# 7. Draw ComplexHeatmap

ht <- Heatmap(
  final_matrix,
  name = "Percent(%)",
  col = col_heatmap,
  
  top_annotation = HeatmapAnnotation(
    Phylum = mat_col_annotation$Phylum_Group,
    Parasitism = mat_col_annotation$Parasitism_type,
    col = all_colors,
    annotation_legend_param = list(
      Phylum = list(
        title_gp = gpar(fontsize = 18, fontface = "bold"),
        labels_gp = gpar(fontsize = 16)
      ),
      Parasitism = list(
        title = "Parasitism type",
        title_gp = gpar(fontsize = 18, fontface = "bold"),
        labels_gp = gpar(fontsize = 16)
      )
    )
  ),
  
  bottom_annotation = HeatmapAnnotation(
    "Total retained HGTs" = anno_barplot(
      mat_bottom_data$Total,
      height = unit(3, "cm"),
      gp = gpar(fill = "#8491B4", col = "transparent"),
      axis_param = list(
        direction = "normal",
        gp = gpar(fontsize = 12)
      )
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 16, fontface = "bold")
  ),
  
  left_annotation = rowAnnotation(
    Func_Class = mat_row_annotation$COG_Group,
    col = all_colors,
    annotation_legend_param = list(
      Func_Class = list(
        title_gp = gpar(fontsize = 18, fontface = "bold"),
        labels_gp = gpar(fontsize = 16)
      )
    )
  ),
  
  row_split = mat_row_annotation$COG_Group,
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  row_gap = unit(2, "mm"),
  border = TRUE,
  rect_gp = gpar(col = "grey95", lwd = 0.2),
  
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 16, fontface = "italic"),
  row_title_gp = gpar(fontsize = 18, fontface = "bold"),
  column_title_gp = gpar(fontsize = 18, fontface = "bold"),
  
  heatmap_legend_param = list(
    title = "HGT percentage\n(%)",
    title_gp = gpar(fontsize = 18, fontface = "bold"),
    labels_gp = gpar(fontsize = 16),
    legend_height = unit(4, "cm")
  ),
  
  column_names_max_height = unit(8, "cm"),
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "right"
)


# 8. Save figure and source table

dir.create("revision/plots_retained_HGT", recursive = TRUE, showWarnings = FALSE)

svglite(
  "revision/plots_retained_HGT/FigureS_COG_heatmap_species_retained_HGT_ComplexHeatmap.svg",
  width = 16,
  height = 8,
  fix_text_size = F
)
draw(
  ht,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

pdf(
  "revision/plots_retained_HGT/FigureS_COG_heatmap_species_retained_HGT_ComplexHeatmap.pdf",
  width = 16,
  height = 8
)
draw(
  ht,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

write.table(
  clean_data,
  "revision/plots_retained_HGT/HGT_COG_function_table_retained.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  tog_norm,
  "revision/plots_retained_HGT/HGT_COG_species_percent_matrix_long.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

draw(
  ht,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
#heatmap_df <- processed_cog_data %>%
#  count(Specie, Parasite_clade, COG_category, name = "Count") %>%
#  group_by(Specie) %>%
#  mutate(Total = sum(Count), Proportion = Count / Total) %>%
#  ungroup()
#
#p_cog_heat <- heatmap_df %>%
#  mutate(Specie = fct_reorder(Specie, Total, .desc = FALSE)) %>%
#  ggplot(aes(x = COG_category, y = Specie, fill = Proportion)) +
#  geom_tile(color = "white", linewidth = 0.1) +
#  facet_grid(Parasite_clade ~ ., scales = "free_y", space = "free_y") +
#  scale_fill_gradient(low = "white", high = "#B2182B", labels = percent, name = "Within-species\nproportion") +
#  labs(x = "COG category", y = NULL, title = "Species-resolved functional profiles of retained HGT candidates") +
#  theme_bw() +
#  big_text_theme +
#  theme(axis.text.y = element_text(face = "italic", size = 7), axis.text.x = element_text(face = "bold", size = 12), legend.position = "right")
#
#print(p_cog_heat)
#save_plot(p_cog_heat, "FigureS_COG_heatmap_species_retained_HGT", 10, 18)

# --------------------------- 12. Cryptosporidium panels -------------
crypto_data <- integrate_data %>% filter(str_detect(Specie, "Cryptosporidium"))
if (nrow(crypto_data) > 0) {
  p_crypto_donor <- draw_sunburst(prepare_sunburst_data(crypto_data), "Cryptosporidium spp. retained HGT donor affinity")
  save_plot(p_crypto_donor, "FigureS_Cryptosporidium_donor_sunburst_retained_HGT", 7, 7)

  crypto_cog <- crypto_data %>% add_cog_group() %>% filter(!is.na(COG_Group), str_detect(COG_category, "^[A-Z]$"))
  if (nrow(crypto_cog) > 0) {
    crypto_order <- crypto_cog %>% count(COG_Group, COG_category, name = "Global_Count") %>% arrange(desc(Global_Count))
    p_crypto_cog <- crypto_cog %>%
      count(COG_Group, COG_category, data_type_overlap, name = "Count") %>%
      left_join(crypto_order, by = c("COG_Group", "COG_category")) %>%
      mutate(COG_category = fct_reorder(COG_category, Global_Count, .desc = TRUE)) %>%
      ggplot(aes(x = COG_category, y = Count, fill = data_type_overlap)) +
      geom_col(position = "stack", color = "black", linewidth = 0.15, width = 0.8) +
      scale_fill_manual(values = discovery_colors, name = "Discovery mode") +
      facet_grid(~ COG_Group, scales = "free_x", space = "free_x", switch = "x", labeller = labeller(COG_Group = label_wrap_gen(width = 15))) +
      labs(x = "COG categories", y = "No. of HGT candidates", title = "Cryptosporidium spp. HGT functional distribution") +
      theme_bw() +
      big_text_theme +
      theme(axis.text.x = element_text(face = "bold", size = 12), strip.placement = "outside",
            strip.background = element_rect(fill = "grey95", color = NA), panel.grid.major.x = element_blank())
    save_plot(p_crypto_cog, "FigureS_Cryptosporidium_COG_retained_HGT", 12, 6)
  }
}

# --------------------------- 13. Contig-context QC -------------------
p_contig_tier <- integrate_data %>%
  count(contig_context_tier, manuscript_HGT_category, name = "Count") %>%
  mutate(contig_context_tier = fct_reorder(contig_context_tier, Count, .fun = sum)) %>%
  ggplot(aes(x = contig_context_tier, y = Count, fill = manuscript_HGT_category)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.15) +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = hgt_category_colors, name = "HGT pattern") +
  labs(x = NULL, y = "Number of retained HGT candidates", title = "Contig-context support of retained HGT candidates") +
  theme_bw() +
  big_text_theme +
  theme(axis.text.y = element_text(size = 13, face = "bold"))

p_contig_tier
save_plot(p_contig_tier, "FigureS_contig_context_tier_retained_HGT", 10, 6)

p_contig_length <- ggplot(integrate_data, aes(x = Contig_Length, fill = contig_context_tier)) +
  geom_histogram(bins = 60, alpha = 0.8) +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(x = "Contig length (bp, log10)", y = "Number of HGT candidates", fill = "Contig-context tier",
       title = "Contig length distribution of retained HGT candidates") +
  theme_bw() +
  big_text_theme
p_contig_length
save_plot(p_contig_length, "FigureS_contig_length_distribution_retained_HGT", 10, 6)


# --------------------------- 14. Downstream analyses hook ------------
# For downstream external tables such as Ka/Ks, insertion time, RNA expression, or synteny,
# always filter the external table by retained HGT identifiers before reusing the previous
# plotting blocks. The preferred keys are hgt_id and Accession.
retained_hgt_ids <- integrate_data %>% distinct(hgt_id, Accession, manuscript_HGT_category, data_type_overlap)
write_tsv(retained_hgt_ids, file.path(outdir, "retained_HGT_id_mapping_for_downstream_filters.tsv"))

# Example pattern:
# kaks <- read_tsv("your_kaks_table.tsv") %>%
#   semi_join(retained_hgt_ids, by = "Accession") %>%
#   left_join(retained_hgt_ids, by = "Accession")
# Then reuse the old Ka/Ks plot code with this filtered object.

message("Done. Revised retained-HGT plots and summary tables were written to: ", outdir)

# --------------------------- 15. Ka/Ks with MRCA-------------------------------

#
# Ka/Ks vs MRCA-based retention age
# using final 8,985 retained HGTs
# No internal HGT evidence grading


library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(svglite)

dir.create("revision/plots_retained_HGT", recursive = TRUE, showWarnings = FALSE)

if (!exists("big_text_theme")) {
  big_text_theme <- theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 18),
    panel.grid = element_blank(),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16, face = "bold")
  )
}

# 1. Read final retained HGT dataset
integrate_data <- read_tsv(
  "revision/main_phylogeny_retained_HGT_8985.tsv",
  show_col_types = FALSE
) %>%
  mutate(
    Accession = as.character(Accession),
    
    # Discovery_mode only indicates the original discovery route.
    # It is NOT an evidence grade.
    Discovery_mode = case_when(
      "data_type_overlap" %in% colnames(.) &
        data_type_overlap %in% c("Only Protein", "Protein-only", "Protein only") ~ "Only Protein",
      
      "data_type_overlap" %in% colnames(.) &
        data_type_overlap %in% c("Only Nucleotide", "Nucleotide-only", "Nucleotide only") ~ "Only Nucleotide",
      
      "data_type_overlap" %in% colnames(.) &
        str_detect(as.character(data_type_overlap), "Overlap") ~ "Overlap (Protein-ref)",
      
      TRUE ~ "Unclassified"
    ),
    
    Discovery_mode = factor(
      Discovery_mode,
      levels = c("Only Protein", "Overlap (Protein-ref)", "Only Nucleotide", "Unclassified")
    )
  )

retained_ids <- integrate_data %>%
  select(Accession, Discovery_mode, everything())

cat("Final retained HGT candidates:", n_distinct(retained_ids$Accession), "\n")


# 2. Read MRCA age and Ka/Ks data
time_data <- read.csv(
  "../Final_HGT_Robust_TBLASTN.csv",
  stringsAsFactors = FALSE
) %>%
  mutate(
    Gene_ID = as.character(Gene_ID),
    Age_MYA = as.numeric(Age_MYA)
  )

kaks_data <- read.delim(
  "../all_HGT.prot.kaks.tsv",
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Sequence = as.character(Sequence),
    `Ka/Ks` = as.numeric(`Ka/Ks`),
    Ka = as.numeric(Ka),
    Ks = as.numeric(Ks)
  )


# 3. Keep only final retained HGTs and merge age + Ka/Ks
retained_time <- retained_ids %>%
  inner_join(time_data, by = c("Accession" = "Gene_ID"))
retained_kaks <- retained_ids %>%
  inner_join(kaks_data, by = c("Accession" = "Sequence"))

merged_data <- retained_time %>%
  inner_join(kaks_data, by = c("Accession" = "Sequence")) %>%
  mutate(
    Age_MYA = as.numeric(Age_MYA),
    `Ka/Ks` = as.numeric(`Ka/Ks`),
    Ka = as.numeric(Ka),
    Ks = as.numeric(Ks)
  ) %>%
  filter(
    !is.na(Age_MYA),
    !is.na(`Ka/Ks`),
    is.finite(Age_MYA),
    is.finite(`Ka/Ks`),
    Age_MYA > 0,
    abs(Age_MYA - 1448.548) > 1e-6,
    `Ka/Ks` < 5
  ) %>%
  mutate(
    Retention_age_class = case_when(
      Age_MYA > 541 ~ "Precambrian",
      Age_MYA > 252 ~ "Paleozoic",
      Age_MYA > 66 ~ "Mesozoic",
      Age_MYA > 0 ~ "Cenozoic",
      TRUE ~ "Tip / unresolved"
    ),
    Retention_age_class = factor(
      Retention_age_class,
      levels = c("Cenozoic", "Mesozoic", "Paleozoic", "Precambrian")
    )
  )

cat("Retained HGTs with MRCA age:", n_distinct(retained_time$Accession), "\n")
cat("Retained HGTs with both MRCA age and Ka/Ks:", n_distinct(merged_data$Accession), "\n")
cat("Rows used for Ka/Ks plot:", nrow(merged_data), "\n\n")

cat("Discovery mode distribution in plotted data:\n")
print(table(merged_data$Discovery_mode, useNA = "ifany"))

write_tsv(
  merged_data,
  "revision/plots_retained_HGT/Table_final_retained_HGT_MRCA_KaKs_merged.tsv"
)


# 4. Scatter plot
p_scatter <- ggplot(
  merged_data,
  aes(x = Age_MYA, y = `Ka/Ks`)
) +
  geom_point(
    aes(),
    alpha = 0.55,
    color = "grey25",
    shape = 21,
    size = 2
  ) +
  geom_smooth(
    method = "loess",
    color = "#E74C3C",
    fill = "#FADBD8",
    linewidth = 1.2,
    se = TRUE
  ) +
  stat_cor(
    method = "spearman",
    label.x.npc = "right",
    label.y.npc = "top",
    hjust = 1,
    size = 6,
    color = "black",
    fontface = "bold"
  ) +
  theme_classic() +
  labs(
    x = "MRCA-based minimum retention age (MYA)",
    y = "Selective pressure (Ka/Ks)"
  ) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  ) +
  big_text_theme +
  theme(
    legend.position = "top",
    axis.text = element_text(color = "black", face = "bold")
  )

print(p_scatter)

svglite(
  "revision/plots_retained_HGT/Figure6_final_retained_HGT_MRCA_age_KaKs_scatter.svg",
  width = 9,
  height = 7,
  fix_text_size = F
)
print(p_scatter)
dev.off()




# 5. Optional: Ka/Ks by MRCA-based retention age class
p_box_age <- merged_data %>
  %
  filter(!is.na(Retention_age_class)) %>%
  ggplot(aes(x = Retention_age_class, y = `Ka/Ks`, fill = Retention_age_class)) +
  geom_boxplot(
    alpha = 0.75,
    outlier.size = 0.6,
    outlier.alpha = 0.35,
    linewidth = 0.35
  ) +
  geom_jitter(
    width = 0.18,
    alpha = 0.25,
    size = 0.7,
    color = "grey20"
  ) +
  stat_compare_means(
    method = "kruskal.test",
    label.y.npc = "top",
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "Cenozoic" = "#91BFDB",
      "Mesozoic" = "#FFFFBF",
      "Paleozoic" = "#FC8D59",
      "Precambrian" = "#D73027"
    ),
    name = "MRCA-based\nretention age"
  ) +
  theme_bw() +
  labs(
    x = "MRCA-based retention age class",
    y = "Selective pressure (Ka/Ks)"
  ) +
  big_text_theme +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold")
  )

print(p_box_age)

svglite(
  "revision/plots_retained_HGT/Figure6_final_retained_HGT_KaKs_by_MRCA_age_class.svg",
  width = 8,
  height = 7,
  fix_text_size = F
)
print(p_box_age)
dev.off()

#

# --------------------------- 16. Donor taxonomy PcoA and niches---------------
# Parasite niche donor distribution using retained HGTs and HGT_target_lineage.txt

library(tidyverse)
library(vegan)
library(patchwork)
library(ggsci)
library(svglite)
library(scales)

dir.create("revision/plots_retained_HGT", recursive = TRUE, showWarnings = FALSE)

if (!exists("big_text_theme")) {
  big_text_theme <- theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 18),
    panel.grid = element_blank(),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16, face = "bold")
  )
}

parse_short_lineage <- function(short_lineage) {
  if (is.na(short_lineage) || short_lineage == "") {
    x <- rep(NA_character_, 7)
  } else {
    x <- strsplit(short_lineage, ";", fixed = TRUE)[[1]]
    x <- str_squish(x)
    if (length(x) < 7) {
      x <- c(x, rep(NA_character_, 7 - length(x)))
    }
    if (length(x) > 7) {
      x <- x[1:7]
    }
  }
  
  tibble(
    Domain = x[1],
    Phylum = x[2],
    Class = x[3],
    Order = x[4],
    Family = x[5],
    Genus = x[6],
    Species_donor = x[7]
  )
}

integrate_data <- readr::read_tsv(
  "revision/main_phylogeny_retained_HGT_8985.tsv",
  show_col_types = FALSE
) %>%
  mutate(
    Accession = as.character(Accession),
    Data_type = as.character(Data_type),
    data_type_overlap = as.character(data_type_overlap)
  )

target_lineage <- read.delim(
  "HGT_target_lineage.txt",
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  fill = TRUE,
  colClasses = "character"
)

colnames(target_lineage)[1:3] <- c(
  "Target_taxid_from_lineage_file",
  "Target_full_lineage",
  "Target_short_lineage"
)

if (nrow(target_lineage) != nrow(integrate_data)) {
  stop(
    paste0(
      "Row number mismatch: integrate_data has ",
      nrow(integrate_data),
      " rows, but HGT_target_lineage.txt has ",
      nrow(target_lineage),
      " rows. Please check whether both files are in the same order."
    )
  )
}

donor_taxonomy <- map_dfr(
  target_lineage$Target_short_lineage,
  parse_short_lineage
)

hgt_data <- bind_cols(
  integrate_data,
  target_lineage %>%
    select(Target_taxid_from_lineage_file, Target_full_lineage, Target_short_lineage),
  donor_taxonomy
) %>%
  mutate(
    across(
      c(Domain, Phylum, Class, Order, Family, Genus, Species_donor),
      ~ifelse(is.na(.x) | .x == "", "Unclassified", .x)
    )
  )

cat("Total retained HGTs:", nrow(hgt_data), "\n")
cat("Donor phylum distribution before top-N grouping:\n")
print(hgt_data %>% count(Phylum, sort = TRUE) %>% head(20))

hgt_data <- hgt_data %>%
  mutate(
    Parasite_clade_display = case_when(
      Parasite_clade %in% c("Protozoa", "Microbial eukaryotes", "microbial eukaryotes") ~ "Microbial eukaryotes",
      Parasite_clade %in% c("Metazoa", "Metazoan parasites", "metazoan parasites") ~ "Metazoan parasites",
      TRUE ~ as.character(Parasite_clade)
    ),
    Parasitism_type = case_when(
      str_detect(as.character(Parasitism_type), "胞内寄生|Intracellular") ~ "Intracellular",
      str_detect(as.character(Parasitism_type), "胞外寄生|Extracellular") ~ "Extracellular",
      str_detect(as.character(Parasitism_type), "自由生活|Free-living") ~ "Free-living",
      TRUE ~ "Unclassified"
    ),
    Niche = paste(Parasitism_type, Parasite_clade_display, sep = " ")
  ) %>%
  filter(Parasitism_type != "Unclassified") %>%
  filter(!is.na(Specie), Specie != "")

top_phyla <- hgt_data %>%
  count(Phylum, sort = TRUE) %>%
  slice_max(order_by = n, n = 8, with_ties = FALSE) %>%
  pull(Phylum)

phylum_levels <- c(
  setdiff(top_phyla, c("Unclassified", "Others")),
  intersect("Unclassified", top_phyla),
  "Others"
)

niche_order <- hgt_data %>%
  count(Niche, sort = TRUE) %>%
  pull(Niche)

hgt_data <- hgt_data %>%
  mutate(
    Phylum_group = ifelse(Phylum %in% top_phyla, Phylum, "Others"),
    Phylum_group = factor(Phylum_group, levels = phylum_levels),
    Niche = factor(Niche, levels = niche_order)
  )

bar_data <- hgt_data %>%
  count(Niche, Phylum_group, name = "Count") %>%
  group_by(Niche) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

fig_niche_bar <- ggplot(bar_data, aes(x = Niche, y = Proportion, fill = Phylum_group)) +
  geom_bar(
    stat = "identity",
    position = "fill",
    width = 0.7,
    color = "black",
    linewidth = 0.2
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_d3(name = "Donor phylum") +
  labs(
    x = "Parasite ecological niche",
    y = "Relative abundance of donor phyla"
  ) +
  theme_classic() +
  big_text_theme +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

fig_niche_bar

otu_table <- hgt_data %>%
  count(Specie, Phylum, name = "Count") %>%
  pivot_wider(names_from = Phylum, values_from = Count, values_fill = 0) %>%
  column_to_rownames(var = "Specie")

otu_table <- otu_table[rowSums(otu_table) > 0, , drop = FALSE]

otu_table_rel <- vegan::decostand(otu_table, method = "total")

metadata <- hgt_data %>%
  select(Specie, Niche, Parasite_clade_display, Parasitism_type) %>%
  distinct(Specie, .keep_all = TRUE) %>%
  column_to_rownames(var = "Specie")

metadata <- metadata[rownames(otu_table_rel), , drop = FALSE]
metadata$Niche <- factor(metadata$Niche)

bray_dist <- vegan::vegdist(otu_table_rel, method = "bray")

pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)

pcoa_points <- as.data.frame(pcoa_res$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points$Specie <- rownames(pcoa_points)

pcoa_points <- pcoa_points %>%
  left_join(
    metadata %>%
      rownames_to_column("Specie"),
    by = "Specie"
  )

positive_eig <- pcoa_res$eig[pcoa_res$eig > 0]
pc1_explained <- round(pcoa_res$eig[1] / sum(positive_eig) * 100, 1)
pc2_explained <- round(pcoa_res$eig[2] / sum(positive_eig) * 100, 1)

adonis_res <- vegan::adonis2(bray_dist ~ Niche, data = metadata, permutations = 999)

p_value <- adonis_res$`Pr(>F)`[1]
r_squared <- round(adonis_res$R2[1], 3)

p_label <- ifelse(
  p_value < 0.001,
  "p < 0.001",
  paste0("p = ", signif(p_value, 3))
)

stat_label <- sprintf("PERMANOVA:\nR² = %s\n%s", r_squared, p_label)

ellipse_data <- pcoa_points %>%
  group_by(Niche) %>%
  summarise(n = n(), .groups = "drop")

niche_for_ellipse <- ellipse_data %>%
  filter(n >= 3) %>%
  pull(Niche)

fig_niche_pcoa_base <- ggplot(
  pcoa_points,
  aes(x = PCoA1, y = PCoA2, color = Niche, shape = Niche)
) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_nejm(name = "Ecological niche") +
  scale_shape_discrete(name = "Ecological niche") +
  labs(
    x = paste0("PCoA1 (", pc1_explained, "%)"),
    y = paste0("PCoA2 (", pc2_explained, "%)")
  ) +
  theme_bw() +
  big_text_theme +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  annotate(
    "text",
    x = min(pcoa_points$PCoA1, na.rm = TRUE),
    y = min(pcoa_points$PCoA2, na.rm = TRUE),
    label = stat_label,
    hjust = 0,
    vjust = 0,
    size = 5,
    fontface = "italic"
  )

if (length(niche_for_ellipse) > 0) {
  fig_niche_pcoa <- fig_niche_pcoa_base +
    stat_ellipse(
      data = pcoa_points %>% filter(Niche %in% niche_for_ellipse),
      aes(fill = Niche),
      geom = "polygon",
      alpha = 0.1,
      show.legend = FALSE
    ) +
    scale_fill_nejm()
} else {
  fig_niche_pcoa <- fig_niche_pcoa_base
}

fig_niche_pcoa

composite_fig <- (fig_niche_bar + fig_niche_pcoa) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(t = 0.5, r = 2, b = 0.5, l = 2, unit = "cm")
    )
  )

print(composite_fig)

svglite(
  "revision/plots_retained_HGT/FigureS_retained_HGT_donor_phylum_by_parasite_niche.svg",
  width = 15,
  height = 7,
  fix_text_size = F
)
print(composite_fig)
dev.off()


write_tsv(
  hgt_data,
  "revision/plots_retained_HGT/Table_retained_HGT_with_HGT_target_lineage_taxonomy.tsv"
)

write_tsv(
  bar_data,
  "revision/plots_retained_HGT/Table_retained_HGT_donor_phylum_by_niche.tsv"
)

pcoa_export <- pcoa_points %>%
  mutate(
    PERMANOVA_R2 = r_squared,
    PERMANOVA_p = p_value
  )

write_tsv(
  pcoa_export,
  "revision/plots_retained_HGT/Table_retained_HGT_donor_phylum_PCoA_points.tsv"
)

# --------------------------- 17. MRCA time analysis----------------------------
## HGT MRCA-based retention age and COG functional composition

library(tidyverse)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(scales)

dir.create("revision/plots_retained_HGT", recursive = TRUE, showWarnings = FALSE)

if (!exists("big_text_theme")) {
  big_text_theme <- theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 18),
    panel.grid = element_blank(),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16, face = "bold")
  )
}

## 1. Read final retained HGT dataset and MRCA-based age table

retained_hgt <- readr::read_tsv(
  "revision/main_phylogeny_retained_HGT_8985.tsv",
  show_col_types = FALSE
) %>%
  mutate(
    Accession = as.character(Accession),
    COG_category = as.character(COG_category),
    Parasite_clade = as.character(Parasite_clade)
  )

## Read MRCA-based retention age tables for protein- and nucleotide-represented HGTs

time_data_prot <- read.csv(
  "../Final_HGT_Robust_TBLASTN.csv",
  stringsAsFactors = FALSE
) %>%
  mutate(
    Gene_ID = as.character(Gene_ID),
    Age_MYA = as.numeric(Age_MYA),
    MRCA_source = "TBLASTN",
    MRCA_Data_type = "Protein"
  )

time_data_nucl <- read.csv(
  "../Final_HGT_Robust_TBLASTX.csv",
  stringsAsFactors = FALSE
) %>%
  mutate(
    Gene_ID = as.character(Gene_ID),
    Age_MYA = as.numeric(Age_MYA),
    MRCA_source = "TBLASTX",
    MRCA_Data_type = "Nucleotide"
  )

time_data <- bind_rows(
  time_data_prot,
  time_data_nucl
) %>%
  select(Gene_ID, Age_MYA, MRCA_source, MRCA_Data_type, everything())

## 2. COG functional mapping

cog_map <- c(
  "J" = "Translation",
  "A" = "RNA processing",
  "K" = "Transcription",
  "L" = "Replication",
  "B" = "Chromatin structure",
  "D" = "Cell cycle",
  "Y" = "Nuclear structure",
  "V" = "Defense mechanisms",
  "T" = "Signal transduction",
  "M" = "Cell wall/membrane",
  "N" = "Cell motility",
  "Z" = "Cytoskeleton",
  "W" = "Extracellular structures",
  "U" = "Intracellular trafficking",
  "O" = "Posttranslational modification",
  "C" = "Energy production",
  "G" = "Carbohydrate transport",
  "E" = "Amino acid transport",
  "F" = "Nucleotide transport",
  "H" = "Coenzyme transport",
  "I" = "Lipid transport",
  "P" = "Inorganic ion transport",
  "Q" = "Secondary metabolites",
  "R" = "General function prediction",
  "S" = "Function unknown"
)

era_levels <- c(
  "Precambrian\n(>541 MYA)",
  "Paleozoic\n(252-541 MYA)",
  "Mesozoic\n(66-252 MYA)",
  "Cenozoic\n(0-66 MYA)"
)

standardize_parasite_group <- function(x) {
  x <- as.character(x)
  
  case_when(
    str_detect(x, regex("protozoa|microbial eukaryote|microbial_eukaryote", ignore_case = TRUE)) ~ "Protozoa",
    str_detect(x, regex("metazoa|metazoan", ignore_case = TRUE)) ~ "Metazoa",
    TRUE ~ "Other"
  )
}

## 3. Merge retained HGTs with MRCA-based age table

## Merge final retained HGTs with matched MRCA-based retention ages

analysis_base <- retained_hgt %>%
  mutate(
    Data_type = as.character(Data_type),
    Data_type_for_MRCA = case_when(
      Data_type == "Protein" ~ "Protein",
      Data_type == "Nucleotide" ~ "Nucleotide",
      TRUE ~ Data_type
    )
  ) %>%
  inner_join(
    time_data,
    by = c(
      "Accession" = "Gene_ID",
      "Data_type_for_MRCA" = "MRCA_Data_type"
    )
  ) %>%
  mutate(
    Age_MYA = as.numeric(Age_MYA),
    Parasite_group = standardize_parasite_group(Parasite_clade),
    MRCA_retention_age_class = case_when(
      Age_MYA <= 66 ~ "Cenozoic\n(0-66 MYA)",
      Age_MYA > 66 & Age_MYA <= 252 ~ "Mesozoic\n(66-252 MYA)",
      Age_MYA > 252 & Age_MYA <= 541 ~ "Paleozoic\n(252-541 MYA)",
      Age_MYA > 541 ~ "Precambrian\n(>541 MYA)",
      TRUE ~ NA_character_
    ),
    MRCA_retention_age_class = factor(
      MRCA_retention_age_class,
      levels = era_levels
    )
  ) %>%
  filter(
    !is.na(Age_MYA),
    is.finite(Age_MYA),
    Age_MYA > 0,
    abs(Age_MYA - 1448.548) > 1e-6,
    !is.na(MRCA_retention_age_class)
  )

cat("Final retained HGTs:", n_distinct(retained_hgt$Accession), "\n")
cat("Protein MRCA records:", n_distinct(time_data_prot$Gene_ID), "\n")
cat("Nucleotide MRCA records:", n_distinct(time_data_nucl$Gene_ID), "\n")
cat("Retained HGTs with matched MRCA age:", n_distinct(analysis_base$Accession), "\n\n")

cat("Matched MRCA source distribution:\n")
print(table(analysis_base$MRCA_source, useNA = "ifany"))

cat("\nMatched Data_type distribution:\n")
print(table(analysis_base$Data_type, useNA = "ifany"))

cat("\nMRCA retention age class distribution:\n")
print(table(analysis_base$MRCA_retention_age_class, useNA = "ifany"))

write_tsv(
  analysis_base,
  "revision/plots_retained_HGT/Table_retained_HGT_with_MRCA_retention_age.tsv"
)

## 4. Function for alluvial plot and residual heatmap

run_cog_age_analysis <- function(input_df, group_label, output_prefix, alluvial_width = 10, residual_width = 9) {
  
  message("Running COG retention-age analysis for: ", group_label)
  
  input_df <- input_df %>%
    filter(!is.na(COG_category), COG_category != "-", COG_category != "", COG_category != "NA")
  
  plot_data <- input_df %>%
    mutate(
      COG_category_clean = str_replace_all(COG_category, "[^A-Z]", "")
    ) %>%
    filter(COG_category_clean != "") %>%
    mutate(COG_category_split = str_split(COG_category_clean, "")) %>%
    unnest(COG_category_split) %>%
    filter(COG_category_split %in% names(cog_map)) %>%
    mutate(
      COG_Full = cog_map[COG_category_split]
    ) %>%
    count(MRCA_retention_age_class, COG_Full, name = "Count") %>%
    group_by(MRCA_retention_age_class) %>%
    mutate(
      Percentage = Count / sum(Count) * 100
    ) %>%
    ungroup()
  
  if (nrow(plot_data) == 0) {
    warning("No valid COG data for: ", group_label)
    return(NULL)
  }
  
  plot_data$MRCA_retention_age_class <- factor(
    plot_data$MRCA_retention_age_class,
    levels = era_levels
  )
  
  write_tsv(
    plot_data,
    paste0("revision/plots_retained_HGT/", output_prefix, "_COG_by_MRCA_retention_age_plot_data.tsv")
  )
  
  p_alluvial <- ggplot(
    plot_data,
    aes(
      x = MRCA_retention_age_class,
      y = Percentage,
      fill = COG_Full,
      stratum = COG_Full,
      alluvium = COG_Full
    )
  ) +
    geom_alluvium(
      alpha = 0.4,
      width = 1 / 4,
      color = "white",
      linewidth = 0.1
    ) +
    geom_stratum(
      width = 1 / 4,
      color = "white",
      linewidth = 0.5
    ) +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal() +
    labs(
      x = "MRCA-based minimum retention age class",
      y = "Relative proportion (%)",
      fill = "COG functional category",
      title = group_label
    ) +
    big_text_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold", color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "right"
    )
  
  print(p_alluvial)
  
  svglite(
    paste0("revision/plots_retained_HGT/", output_prefix, "_COG_alluvial_MRCA_retention_age.svg"),
    width = alluvial_width,
    height = 8,
    fix_text_size = F
  )
  print(p_alluvial)
  dev.off()
  
  stat_data <- plot_data %>%
    select(MRCA_retention_age_class, COG_Full, Count) %>%
    pivot_wider(
      names_from = COG_Full,
      values_from = Count,
      values_fill = 0
    ) %>%
    arrange(MRCA_retention_age_class)
  
  stat_matrix <- as.matrix(stat_data[, -1, drop = FALSE])
  rownames(stat_matrix) <- stat_data$MRCA_retention_age_class
  
  if (nrow(stat_matrix) < 2 || ncol(stat_matrix) < 2) {
    warning("Not enough rows or columns for chi-square residual analysis: ", group_label)
    return(list(
      plot_data = plot_data,
      alluvial = p_alluvial,
      residual_heatmap = NULL,
      top_hits = NULL
    ))
  }
  
  chisq_res <- suppressWarnings(chisq.test(stat_matrix))
  
  residuals_df <- as.data.frame(chisq_res$residuals)
  residuals_df$MRCA_retention_age_class <- rownames(residuals_df)
  
  sig_changes <- residuals_df %>%
    pivot_longer(
      cols = -MRCA_retention_age_class,
      names_to = "COG_Category",
      values_to = "Residual"
    ) %>%
    mutate(
      Status = case_when(
        Residual > 1.96 ~ "Enriched",
        Residual < -1.96 ~ "Depleted",
        TRUE ~ "Stable"
      ),
      MRCA_retention_age_class = factor(
        MRCA_retention_age_class,
        levels = era_levels
      )
    )
  
  top_hits <- sig_changes %>%
    filter(Status != "Stable") %>%
    arrange(desc(abs(Residual)))
  
  write_tsv(
    sig_changes,
    paste0("revision/plots_retained_HGT/", output_prefix, "_COG_chisq_standardized_residuals.tsv")
  )
  
  write_tsv(
    top_hits,
    paste0("revision/plots_retained_HGT/", output_prefix, "_COG_significant_residual_hits.tsv")
  )
  
  p_residual <- ggplot(
    sig_changes,
    aes(x = MRCA_retention_age_class, y = COG_Category, fill = Residual)
  ) +
    geom_tile(color = "white") +
    geom_text(
      aes(
        label = round(Residual, 2),
        fontface = ifelse(abs(Residual) > 1.96, "bold", "plain")
      ),
      size = 4.5,
      color = "black"
    ) +
    scale_fill_gradient2(
      low = "#377EB8",
      mid = "white",
      high = "#E41A1C",
      midpoint = 0
    ) +
    theme_minimal() +
    labs(
      x = "MRCA-based minimum retention age class",
      y = "COG functional category",
      fill = "Standardized\nresidual",
      title = group_label
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 16, face = "bold"),
      axis.text.y = element_text(color = "black", size = 14, face = "bold"),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 15, face = "bold")
    ) +
    big_text_theme
  
  print(p_residual)
  
  svglite(
    paste0("revision/plots_retained_HGT/", output_prefix, "_COG_residuals_MRCA_retention_age.svg"),
    width = residual_width,
    height = 8,
    fix_text_size = F
  )
  print(p_residual)
  dev.off()

  
  return(list(
    plot_data = plot_data,
    alluvial = p_alluvial,
    residual_heatmap = p_residual,
    top_hits = top_hits,
    chisq = chisq_res
  ))
}

## 5. Run three analyses: all retained HGTs, protozoa, and metazoa

res_all <- run_cog_age_analysis(
  input_df = analysis_base,
  group_label = "All final retained HGT candidates",
  output_prefix = "All_retained_HGT",
  alluvial_width = 11,
  residual_width = 9
)

res_protozoa <- run_cog_age_analysis(
  input_df = analysis_base %>% filter(Parasite_group == "Protozoa"),
  group_label = "Protozoa / microbial eukaryote retained HGT candidates",
  output_prefix = "Protozoa_retained_HGT",
  alluvial_width = 11,
  residual_width = 9
)

res_metazoa <- run_cog_age_analysis(
  input_df = analysis_base %>% filter(Parasite_group == "Metazoa"),
  group_label = "Metazoa retained HGT candidates",
  output_prefix = "Metazoa_retained_HGT",
  alluvial_width = 9,
  residual_width = 9
)

## 6. Optional overview table for manuscript reporting

summary_by_group_age <- analysis_base %>%
  mutate(
    Analysis_group = case_when(
      Parasite_group == "Protozoa" ~ "Protozoa / microbial eukaryotes",
      Parasite_group == "Metazoa" ~ "Metazoa",
      TRUE ~ "Other"
    )
  ) %>%
  count(Analysis_group, MRCA_retention_age_class, name = "HGT_count") %>%
  group_by(Analysis_group) %>%
  mutate(Proportion = HGT_count / sum(HGT_count)) %>%
  ungroup()

summary_all_age <- analysis_base %>%
  count(MRCA_retention_age_class, name = "HGT_count") %>%
  mutate(
    Analysis_group = "All final retained HGT candidates",
    Proportion = HGT_count / sum(HGT_count)
  ) %>%
  select(Analysis_group, MRCA_retention_age_class, HGT_count, Proportion)

summary_age_final <- bind_rows(summary_all_age, summary_by_group_age)

write_tsv(
  summary_age_final,
  "revision/plots_retained_HGT/Table_retained_HGT_MRCA_retention_age_summary_by_group.tsv"
)

print(summary_age_final)


## MRCA-based retention age distribution for final retained HGTs

library(tidyverse)
library(ggplot2)
library(patchwork)
library(svglite)
library(grid)

dir.create("revision/plots_retained_HGT", recursive = TRUE, showWarnings = FALSE)

if (!exists("big_text_theme")) {
  big_text_theme <- theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 18),
    panel.grid = element_blank(),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 16),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16, face = "bold")
  )
}

retained_file <- "revision/main_phylogeny_retained_HGT_8985.tsv"
prot_mrca_file <- "../Final_HGT_Robust_TBLASTN.csv"
nucl_mrca_file <- "../Final_HGT_Robust_TBLASTX.csv"

retained_hgt <- readr::read_tsv(
  retained_file,
  show_col_types = FALSE
) %>%
  mutate(
    Accession = as.character(Accession),
    Data_type = as.character(Data_type),
    Parasite_clade = as.character(Parasite_clade)
  )

time_data_prot <- read.csv(
  prot_mrca_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Gene_ID = as.character(Gene_ID),
    Age_MYA = as.numeric(Age_MYA),
    MRCA_source = "TBLASTN",
    MRCA_Data_type = "Protein"
  )

time_data_nucl <- read.csv(
  nucl_mrca_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Gene_ID = as.character(Gene_ID),
    Age_MYA = as.numeric(Age_MYA),
    MRCA_source = "TBLASTX",
    MRCA_Data_type = "Nucleotide"
  )

time_data <- bind_rows(
  time_data_prot,
  time_data_nucl
) %>%
  group_by(Gene_ID, MRCA_Data_type) %>%
  slice(1) %>%
  ungroup()

standardize_parasite_group <- function(x) {
  x <- as.character(x)

  case_when(
    str_detect(x, regex("protozoa|microbial eukaryote|microbial_eukaryote", ignore_case = TRUE)) ~ "Protozoa",
    str_detect(x, regex("metazoa|metazoan", ignore_case = TRUE)) ~ "Metazoa",
    TRUE ~ "Other"
  )
}

analysis_base <- retained_hgt %>%
  mutate(
    Data_type_for_MRCA = case_when(
      str_detect(Data_type, regex("^Protein$", ignore_case = TRUE)) ~ "Protein",
      str_detect(Data_type, regex("^Nucleotide$", ignore_case = TRUE)) ~ "Nucleotide",
      TRUE ~ Data_type
    )
  ) %>%
  left_join(
    time_data,
    by = c(
      "Accession" = "Gene_ID",
      "Data_type_for_MRCA" = "MRCA_Data_type"
    )
  ) %>%
  mutate(
    Age_MYA = as.numeric(Age_MYA),
    Parasite_group = standardize_parasite_group(Parasite_clade),
    MRCA_retention_age_class = case_when(
      Age_MYA <= 66 ~ "Cenozoic\n(0-66 MYA)",
      Age_MYA > 66 & Age_MYA <= 252 ~ "Mesozoic\n(66-252 MYA)",
      Age_MYA > 252 & Age_MYA <= 541 ~ "Paleozoic\n(252-541 MYA)",
      Age_MYA > 541 ~ "Precambrian\n(>541 MYA)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Age_MYA),
    is.finite(Age_MYA),
    Age_MYA > 0,
    Age_MYA <= 1500,
    abs(Age_MYA - 1448.548) > 1e-6
  )

cat("Final retained HGTs:", n_distinct(retained_hgt$Accession), "\n")
cat("Retained HGTs with matched MRCA age:", n_distinct(analysis_base$Accession), "\n\n")

cat("MRCA source distribution:\n")
print(table(analysis_base$MRCA_source, useNA = "ifany"))

cat("\nData_type distribution after MRCA matching:\n")
print(table(analysis_base$Data_type, useNA = "ifany"))

cat("\nParasite group distribution after MRCA matching:\n")
print(table(analysis_base$Parasite_group, useNA = "ifany"))

cat("\nMRCA-based retention age class distribution:\n")
print(table(analysis_base$MRCA_retention_age_class, useNA = "ifany"))

write_tsv(
  analysis_base,
  "revision/plots_retained_HGT/Table_final_retained_HGT_MRCA_distribution_input.tsv"
)

parasite_eras <- data.frame(
  Era = c(
    "Cenozoic\n(0-66 MYA)",
    "Mesozoic\n(66-252 MYA)",
    "Paleozoic\n(252-541 MYA)",
    "Precambrian\n(>541 MYA)"
  ),
  start = c(0, 66, 252, 541),
  end = c(66, 252, 541, 1500),
  color = c("#FAD7A1", "#A9DFBF", "#AED6F1", "#E5E7E9")
)

plot_mrca_density <- function(input_df, group_label, output_prefix, x_max = 1500, width = 16, height = 6) {

  plot_df <- input_df %>%
    filter(
      !is.na(Age_MYA),
      is.finite(Age_MYA),
      Age_MYA > 0,
      Age_MYA <= x_max
    )

  if (nrow(plot_df) < 3) {
    warning("Too few records for density plot: ", group_label)
    return(NULL)
  }

  age_summary <- plot_df %>%
    mutate(
      MRCA_retention_age_class = factor(
        MRCA_retention_age_class,
        levels = parasite_eras$Era
      )
    ) %>%
    count(MRCA_retention_age_class, name = "Count") %>%
    mutate(
      Proportion = Count / sum(Count)
    )

  write_tsv(
    age_summary,
    paste0("revision/plots_retained_HGT/", output_prefix, "_MRCA_retention_age_distribution_summary.tsv")
  )

  p_density <- ggplot() +
    geom_rect(
      data = parasite_eras,
      aes(xmin = start, xmax = end, ymin = 0, ymax = Inf, fill = Era),
      alpha = 0.4
    ) +
    scale_fill_manual(values = setNames(parasite_eras$color, parasite_eras$Era)) +
    geom_density(
      data = plot_df,
      aes(x = Age_MYA),
      fill = "#2C3E50",
      color = "#1A252F",
      alpha = 0.6,
      linewidth = 1,
      adjust = 1
    ) +
    geom_rug(
      data = plot_df,
      aes(x = Age_MYA),
      color = "#E74C3C",
      alpha = 0.45,
      sides = "b",
      length = unit(0.03, "npc")
    ) +
    geom_text(
      data = parasite_eras,
      aes(x = (start + end) / 2, y = Inf, label = Era),
      vjust = 1.5,
      size = 5.5,
      fontface = "bold",
      color = "#34495E"
    ) +
    scale_x_reverse(
      breaks = c(0, 66, 252, 541, 1000, 1500),
      limits = c(x_max, 0),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    theme_classic() +
    labs(
      title = group_label,
      subtitle = "Density shows the distribution of MRCA-based minimum retention age assignments",
      x = "MRCA-based minimum retention age (million years ago)",
      y = "Density of retained HGT candidates",
      fill = "Retention age class"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "grey40", margin = margin(b = 15)),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
    ) +
    big_text_theme +
    theme(
      plot.margin = margin(t = 0.5, r = 3, b = 0.5, l = 0.5, unit = "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

  print(p_density)

  svglite(
    paste0("revision/plots_retained_HGT/", output_prefix, "_MRCA_retention_age_density.svg"),
    width = width,
    height = height,
    fix_text_size = F
  )
  print(p_density)
  dev.off()



  return(p_density)
}

p_density_all <- plot_mrca_density(
  input_df = analysis_base,
  group_label = "All final retained HGT candidates",
  output_prefix = "All_retained_HGT",
  x_max = 1500,
  width = 16,
  height = 6
)

p_density_protozoa <- plot_mrca_density(
  input_df = analysis_base %>% filter(Parasite_group == "Protozoa"),
  group_label = "Protozoa / microbial eukaryote retained HGT candidates",
  output_prefix = "Protozoa_retained_HGT",
  x_max = 1500,
  width = 16,
  height = 6
)

p_density_metazoa <- plot_mrca_density(
  input_df = analysis_base %>% filter(Parasite_group == "Metazoa"),
  group_label = "Metazoa retained HGT candidates",
  output_prefix = "Metazoa_retained_HGT",
  x_max = 1500,
  width = 16,
  height = 6
)

combined_density <- p_density_all / p_density_protozoa / p_density_metazoa +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(t = 0.5, r = 1.5, b = 0.5, l = 1.5, unit = "cm")
    )
  )

print(combined_density)


age_distribution_summary <- analysis_base %>%
  mutate(
    MRCA_retention_age_class = factor(
      MRCA_retention_age_class,
      levels = parasite_eras$Era
    ),
    Analysis_group = case_when(
      Parasite_group == "Protozoa" ~ "Protozoa / microbial eukaryotes",
      Parasite_group == "Metazoa" ~ "Metazoa",
      TRUE ~ "Other"
    )
  ) %>%
  count(Analysis_group, MRCA_retention_age_class, name = "Count") %>%
  group_by(Analysis_group) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

age_distribution_all <- analysis_base %>%
  mutate(
    MRCA_retention_age_class = factor(
      MRCA_retention_age_class,
      levels = parasite_eras$Era
    )
  ) %>%
  count(MRCA_retention_age_class, name = "Count") %>%
  mutate(
    Analysis_group = "All final retained HGT candidates",
    Proportion = Count / sum(Count)
  ) %>%
  select(Analysis_group, MRCA_retention_age_class, Count, Proportion)

age_distribution_summary_final <- bind_rows(
  age_distribution_all,
  age_distribution_summary
)

write_tsv(
  age_distribution_summary_final,
  "revision/plots_retained_HGT/Table_final_retained_HGT_MRCA_retention_age_distribution_All_Protozoa_Metazoa.tsv"
)

print(age_distribution_summary_final)

# --------------------------- 18. GC with MRCA----------------------------------
## GC content difference versus MRCA-based retention age

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(svglite)
library(patchwork)

dir.create("revision/plots_retained_HGT", recursive = TRUE, showWarnings = FALSE)

if (!exists("big_text_theme")) {
  big_text_theme <- theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 18),
    panel.grid = element_blank(),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 16),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16, face = "bold")
  )
}

retained_file <- "revision/main_phylogeny_retained_HGT_8985.tsv"
prot_mrca_file <- "../Final_HGT_Robust_TBLASTN.csv"
nucl_mrca_file <- "../Final_HGT_Robust_TBLASTX.csv"
prot_gc_file <- "../all_HGT.prot.CG.txt"
nucl_gc_file <- "../all_HGT.nucl.CG.txt"

if (!file.exists(retained_file) && file.exists("main_phylogeny_retained_HGT_8985.tsv")) {
  retained_file <- "main_phylogeny_retained_HGT_8985.tsv"
}

if (!file.exists(prot_gc_file) && file.exists("../all_HGT.prot.CG.txt")) {
  prot_gc_file <- "../all_HGT.prot.CG.txt"
}

if (!file.exists(nucl_gc_file) && file.exists("../all_HGT.nucl.CG.txt")) {
  nucl_gc_file <- "../all_HGT.nucl.CG.txt"
}

read_gc_file <- function(file, type_label) {
  df <- read.table(
    file,
    header = FALSE,
    fill = TRUE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  colnames(df)[1:5] <- c("Genome", "Accession", "HGT_CG", "Genome_CG", "HGT_diff")
  
  df %>%
    select(Genome, Accession, HGT_CG, Genome_CG, HGT_diff) %>%
    mutate(
      Accession = as.character(Accession),
      Genome = as.character(Genome),
      HGT_CG = as.numeric(HGT_CG),
      Genome_CG = as.numeric(Genome_CG),
      HGT_diff = as.numeric(HGT_diff),
      GC_Data_type = type_label
    )
}

gc_prot <- read_gc_file(prot_gc_file, "Protein")
gc_nucl <- read_gc_file(nucl_gc_file, "Nucleotide")

gc_data <- bind_rows(gc_prot, gc_nucl) %>%
  filter(
    !is.na(Accession),
    !is.na(HGT_CG),
    !is.na(Genome_CG),
    !is.na(HGT_diff)
  ) %>%
  group_by(Accession, GC_Data_type) %>%
  slice(1) %>%
  ungroup()

if (max(gc_data$HGT_CG, na.rm = TRUE) <= 1.5 && max(gc_data$Genome_CG, na.rm = TRUE) <= 1.5) {
  gc_data <- gc_data %>%
    mutate(
      HGT_CG = HGT_CG * 100,
      Genome_CG = Genome_CG * 100,
      HGT_diff = HGT_diff * 100
    )
}

retained_hgt <- readr::read_tsv(
  retained_file,
  show_col_types = FALSE
) %>%
  mutate(
    Accession = as.character(Accession),
    Data_type = as.character(Data_type),
    Parasite_clade = as.character(Parasite_clade)
  )

time_data_prot <- read.csv(
  prot_mrca_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Gene_ID = as.character(Gene_ID),
    Age_MYA = as.numeric(Age_MYA),
    MRCA_source = "TBLASTN",
    MRCA_Data_type = "Protein"
  )

time_data_nucl <- read.csv(
  nucl_mrca_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Gene_ID = as.character(Gene_ID),
    Age_MYA = as.numeric(Age_MYA),
    MRCA_source = "TBLASTX",
    MRCA_Data_type = "Nucleotide"
  )

time_data <- bind_rows(time_data_prot, time_data_nucl) %>%
  group_by(Gene_ID, MRCA_Data_type) %>%
  slice(1) %>%
  ungroup()

standardize_parasite_group <- function(x) {
  x <- as.character(x)
  
  case_when(
    str_detect(x, regex("protozoa|microbial eukaryote|microbial_eukaryote", ignore_case = TRUE)) ~ "Protozoa",
    str_detect(x, regex("metazoa|metazoan", ignore_case = TRUE)) ~ "Metazoa",
    TRUE ~ "Other"
  )
}

analysis_gc_mrca <- retained_hgt %>%
  mutate(
    Data_type_for_MRCA = case_when(
      str_detect(Data_type, regex("^Protein$", ignore_case = TRUE)) ~ "Protein",
      str_detect(Data_type, regex("^Nucleotide$", ignore_case = TRUE)) ~ "Nucleotide",
      TRUE ~ Data_type
    )
  ) %>%
  left_join(
    gc_data,
    by = c(
      "Accession" = "Accession",
      "Data_type_for_MRCA" = "GC_Data_type"
    )
  ) %>%
  left_join(
    time_data,
    by = c(
      "Accession" = "Gene_ID",
      "Data_type_for_MRCA" = "MRCA_Data_type"
    )
  ) %>%
  mutate(
    Age_MYA = as.numeric(Age_MYA),
    HGT_CG = as.numeric(HGT_CG),
    Genome_CG = as.numeric(Genome_CG),
    HGT_diff = as.numeric(HGT_diff),
    Parasite_group = standardize_parasite_group(Parasite_clade),
    MRCA_retention_age_class = case_when(
      Age_MYA <= 66 ~ "Cenozoic",
      Age_MYA > 66 & Age_MYA <= 252 ~ "Mesozoic",
      Age_MYA > 252 & Age_MYA <= 541 ~ "Paleozoic",
      Age_MYA > 541 ~ "Precambrian",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Age_MYA),
    !is.na(HGT_diff),
    is.finite(Age_MYA),
    is.finite(HGT_diff),
    Age_MYA > 0,
    Age_MYA <= 1500,
    abs(Age_MYA - 1448.548) > 1e-6
  )

cat("Final retained HGTs:", n_distinct(retained_hgt$Accession), "\n")
cat("Retained HGTs with GC and MRCA age:", n_distinct(analysis_gc_mrca$Accession), "\n\n")

cat("Data type distribution:\n")
print(table(analysis_gc_mrca$Data_type, useNA = "ifany"))

cat("\nMRCA source distribution:\n")
print(table(analysis_gc_mrca$MRCA_source, useNA = "ifany"))

cat("\nParasite group distribution:\n")
print(table(analysis_gc_mrca$Parasite_group, useNA = "ifany"))

cat("\nMRCA retention age class distribution:\n")
print(table(analysis_gc_mrca$MRCA_retention_age_class, useNA = "ifany"))

write_tsv(
  analysis_gc_mrca,
  "revision/plots_retained_HGT/Table_final_retained_HGT_GC_MRCA_merged.tsv"
)

plot_gc_mrca <- function(input_df, group_label, output_prefix, width = 9, height = 7) {
  plot_df <- input_df %>%
    filter(
      !is.na(Age_MYA),
      !is.na(HGT_diff),
      is.finite(Age_MYA),
      is.finite(HGT_diff)
    )
  
  if (nrow(plot_df) < 5) {
    warning("Too few records for GC-MRCA fitting plot: ", group_label)
    return(NULL)
  }
  
  p <- ggplot(plot_df, aes(x = Age_MYA, y = HGT_diff)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    geom_point(
      aes(fill = Data_type_for_MRCA),
      alpha = 0.55,
      color = "grey25",
      shape = 21,
      size = 2.2
    ) +
    geom_smooth(
      method = "loess",
      color = "#E74C3C",
      fill = "#FADBD8",
      alpha = 0.35,
      linewidth = 1.2,
      se = TRUE
    ) +
    stat_cor(
      method = "spearman",
      label.x.npc = "right",
      label.y.npc = "top",
      hjust = 1,
      size = 6,
      color = "black",
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c(
        "Protein" = "#2E5077",
        "Nucleotide" = "#4DA1A9"
      ),
      name = "Sequence type"
    ) +
    theme_classic() +
    labs(
      title = group_label,
      x = "MRCA-based minimum retention age (MYA)",
      y = "GC difference (HGT GC - genome GC, %)"
    ) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
    ) +
    big_text_theme +
    theme(
      legend.position = "top",
      axis.text = element_text(color = "black", face = "bold")
    )
  
  print(p)
  
  svglite(
    paste0("revision/plots_retained_HGT/", output_prefix, "_GC_difference_vs_MRCA_retention_age.svg"),
    width = width,
    height = height,
    fix_text_size = F
  )
  print(p)
  dev.off()
  
  return(p)
}

p_gc_all <- plot_gc_mrca(
  input_df = analysis_gc_mrca,
  group_label = "All final retained HGT candidates",
  output_prefix = "All_retained_HGT",
  width = 9,
  height = 7
)


