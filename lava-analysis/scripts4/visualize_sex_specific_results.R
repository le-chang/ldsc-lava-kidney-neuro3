#!/usr/bin/env Rscript

# LAVA Sex-Specific Results Visualization

library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(viridis)

# Load significant results
output_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/results/summary_results"
plot_dir <- file.path(output_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)

bivar_sig <- fread(file.path(output_dir, "bivar_significant_results.tsv"))
summary_data <- fread(file.path(output_dir, "summary_by_sex_and_pair.tsv"))

# ============================================================================
# PLOT 1: Number of Significant Loci by Sex and Phenotype Pair
# ============================================================================

p1 <- ggplot(summary_data, aes(x = pheno_pair, y = n_sig_loci, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = n_sig_loci), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("females" = "#E69F00", "males" = "#56B4E9"),
                    labels = c("Females", "Males")) +
  labs(title = "Significant Local Genetic Correlations by Sex",
       subtitle = paste("Bonferroni-corrected threshold applied"),
       x = "Phenotype Pair",
       y = "Number of Significant Loci",
       fill = "Sex") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(plot_dir, "fig1_nsig_by_sex.pdf"), p1, width = 8, height = 6)
ggsave(file.path(plot_dir, "fig1_nsig_by_sex.png"), p1, width = 8, height = 6, dpi = 300)

# ============================================================================
# PLOT 2: Distribution of Genetic Correlations (rho)
# ============================================================================

p2 <- ggplot(bivar_sig, aes(x = rho, fill = sex)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("females" = "#E69F00", "males" = "#56B4E9"),
                    labels = c("Females", "Males")) +
  facet_wrap(~pheno_pair, scales = "free_y", ncol = 1) +
  labs(title = "Distribution of Significant Local Genetic Correlations",
       x = "Local Genetic Correlation (ρ)",
       y = "Count",
       fill = "Sex") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "lightgray"))

ggsave(file.path(plot_dir, "fig2_rho_distribution.pdf"), p2, width = 10, height = 8)
ggsave(file.path(plot_dir, "fig2_rho_distribution.png"), p2, width = 10, height = 8, dpi = 300)

# ============================================================================
# PLOT 3: Effect Size Comparison (Violin Plot)
# ============================================================================

p3 <- ggplot(bivar_sig, aes(x = pheno_pair, y = rho, fill = sex)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), 
               outlier.shape = NA, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("females" = "#E69F00", "males" = "#56B4E9"),
                    labels = c("Females", "Males")) +
  labs(title = "Effect Size Distribution by Sex and Phenotype Pair",
       x = "Phenotype Pair",
       y = "Local Genetic Correlation (ρ)",
       fill = "Sex") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(plot_dir, "fig3_rho_violin.pdf"), p3, width = 10, height = 6)
ggsave(file.path(plot_dir, "fig3_rho_violin.png"), p3, width = 10, height = 6, dpi = 300)

# ============================================================================
# PLOT 4: P-value Distribution (Manhattan-style)
# ============================================================================

bivar_sig_plot <- bivar_sig %>%
  arrange(pheno_pair, chr, pos) %>%
  group_by(pheno_pair, sex) %>%
  mutate(
    index = row_number(),
    neglog10p = -log10(p)
  )

p4 <- ggplot(bivar_sig_plot, aes(x = index, y = neglog10p, color = as.factor(chr))) +
  geom_point(alpha = 0.7, size = 2) +
  facet_grid(sex ~ pheno_pair, scales = "free_x") +
  scale_color_viridis_d(option = "turbo") +
  labs(title = "Significance of Local Genetic Correlations",
       subtitle = "Significant loci shown by chromosome",
       x = "Locus Index",
       y = "-log10(P-value)",
       color = "Chr") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "lightgray"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(file.path(plot_dir, "fig4_pvalue_manhattan.pdf"), p4, width = 14, height = 8)
ggsave(file.path(plot_dir, "fig4_pvalue_manhattan.png"), p4, width = 14, height = 8, dpi = 300)

# ============================================================================
# PLOT 5: Proportion of Positive vs Negative Correlations
# ============================================================================

direction_data <- summary_data %>%
  select(sex, pheno_pair, n_positive_rg, n_negative_rg) %>%
  pivot_longer(cols = c(n_positive_rg, n_negative_rg),
               names_to = "direction",
               values_to = "count") %>%
  mutate(direction = ifelse(direction == "n_positive_rg", "Positive", "Negative"))

p5 <- ggplot(direction_data, aes(x = pheno_pair, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 3.5) +
  facet_wrap(~sex, ncol = 2) +
  scale_fill_manual(values = c("Positive" = "#00BA38", "Negative" = "#F8766D")) +
  labs(title = "Direction of Significant Local Genetic Correlations",
       x = "Phenotype Pair",
       y = "Number of Loci",
       fill = "Direction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "lightgray"))

ggsave(file.path(plot_dir, "fig5_direction_barplot.pdf"), p5, width = 10, height = 6)
ggsave(file.path(plot_dir, "fig5_direction_barplot.png"), p5, width = 10, height = 6, dpi = 300)

# ============================================================================
# PLOT 6: Sex Comparison Scatter Plot
# ============================================================================

# Reshape data for comparison
if (nrow(bivar_sig) > 0) {
  comparison_data <- bivar_sig %>%
    select(locus, pheno_pair, sex, rho, p) %>%
    pivot_wider(
      names_from = sex,
      values_from = c(rho, p),
      names_sep = "_"
    ) %>%
    filter(!is.na(rho_females) & !is.na(rho_males))
  
  if (nrow(comparison_data) > 0) {
    p6 <- ggplot(comparison_data, aes(x = rho_females, y = rho_males)) +
      geom_point(aes(color = pheno_pair), alpha = 0.7, size = 3) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5) +
      labs(title = "Sex Comparison of Local Genetic Correlations",
           subtitle = "Loci significant in both sexes",
           x = "ρ (Females)",
           y = "ρ (Males)",
           color = "Phenotype Pair") +
      theme_bw() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 14))
    
    ggsave(file.path(plot_dir, "fig6_sex_comparison.pdf"), p6, width = 8, height = 8)
    ggsave(file.path(plot_dir, "fig6_sex_comparison.png"), p6, width = 8, height = 8, dpi = 300)
  }
}

cat("\nAll plots saved to:", plot_dir, "\n")
