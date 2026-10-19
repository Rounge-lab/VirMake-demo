
setwd("~/VirMake-demo") ## Correct path

library(tidyverse)
library(gtsummary)
library(cowplot)
library(paletteer)

source("scripts/utils.R")


# load data ---------------------------------------------------------------

buck_samples <- read_tsv("data/buck/sample_stats-buck-3.tsv", col_types = cols())
buck_genomes <- read_tsv("data/buck/vOTU_stats-buck-3.tsv", col_types = cols())
buck_abundance <- read_tsv("data/buck/rel_abundance_table-buck-3.tsv", col_types = cols())
buck_sample_info <- read_tsv("data/buck/filereport_read_run_PRJEB38681_tsv.txt", col_types = cols())

han_samples <- read_tsv("data/han/sample_stats-han-2.tsv", col_types = cols())
han_genomes <- read_tsv("data/han/vOTU_stats-han-2.tsv", col_types = cols())
han_abundance <- read_tsv("data/han/rel_abundance_table-han-2.tsv", col_types = cols())

moreno_samples <- read_tsv("data/moreno/sample_stats-moreno-5.tsv", col_types = cols())
moreno_genomes <- read_tsv("data/moreno/vOTU_stats-moreno-5.tsv", col_types = cols())
moreno_abundance <- read_tsv("data/moreno/rel_abundance_table-moreno-5.tsv", col_types = cols())

# res_usage <- read_csv2("data/time_res.csv", col_types = cols())
res_usage <- read_csv2("data/time_res-rev.csv", col_types = cols())

res_usage_by_module <- read_delim("data/virmake-res-revision.csv", delim = ";", col_types = cols())

# summaries ---------------------------------------------------------------

## Resource usage plot
res_usage_plot <- 
  res_usage_by_module %>% 
  filter(!str_detect(Process, "Memory|total")) %>% 
  pivot_longer(-Process, names_to = "Dataset", values_to = "CPU hours") %>% 
  mutate(Process = factor(Process, levels = res_usage_by_module$Process[-c(7,8)])) %>% 
  ggplot(aes(x = Dataset, y = `CPU hours`)) +
  geom_col(fill = "#6ba3d6") +
  facet_wrap(~Process, scales = "free_y") +
  theme_minimal() +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 6))


ggsave2(filename = "results/plots/resource_use_plot.pdf", 
        plot = res_usage_plot,
        height = 60, width = 80, units = "mm")

## Summary tables
output_by_sample_table <- 
  moreno_samples %>% 
  mutate(dataset = "Moreno-Gallego") %>% 
  bind_rows(han_samples %>% 
              mutate(dataset = "Li")) %>% 
  bind_rows(buck_samples %>%
              inner_join(buck_sample_info %>%
                           filter(str_detect(scientific_name, "metagenome")) %>%
                           select(sample_id = run_accession), by = "sample_id") %>%
              mutate(dataset = "Buck")) %>%
  mutate(dataset = factor(dataset, levels = c("Moreno-Gallego", "Li", "Buck"))) %>% 
  summarize_samples(by = "dataset")

output_by_genome_summary <-
  moreno_genomes %>% 
  mutate(dataset = "Moreno-Gallego") %>% 
  bind_rows(han_genomes %>% 
              mutate(dataset = "Li")) %>% 
  bind_rows(buck_genomes %>%
              mutate(dataset = "Buck")) %>%
  mutate(dataset = factor(dataset, levels = c("Moreno-Gallego", "Li", "Buck"))) %>% 
  summarize_genomes(by_var = "dataset")

theme_gtsummary_compact()
results_tab <- 
  tbl_stack(list("Resource usage" = summarize_resource_usage(res_usage) %>% modify_header(all_stat_cols() ~ "**{level}**"), 
               "by sample" = output_by_sample_table %>% modify_header(label = ""),
               "by vOTU" = output_by_genome_summary), 
          group_header = c("Resource usage", "Summary of output by sample", "Summary of output by vOTU")) %>%
  modify_footnote(all_stat_cols() ~ "Totals for resource usage; median (Q1,Q3) for sample-level summary; n (%) for vOTU-level summary") %>% 
  as_gt() %>% 
  gt::tab_style(
    style = gt::cell_text(style = "italic"),
    locations = gt::cells_row_groups(groups = everything())
  ) 
results_tab %>% 
  gt::gtsave(filename = "results/tables/results_summary-2.html")
results_tab %>% 
  gt::gtsave(filename = "results/tables/results_summary-2.rtf")
results_tab %>% 
  gt::gtsave(filename = "results/tables/results_summary-2.docx")
reset_gtsummary_theme()


sample_input_output_plot <-
  list(list(moreno_samples, NA), 
     list(han_samples, NA),
     list(buck_samples %>%
            left_join(buck_sample_info %>% select(sample_id = run_accession, scientific_name), by = "sample_id") %>%
            filter(str_detect(scientific_name, "metagenome")))) %>%
  lapply(function(x) plot_input_v_output(x[[1]])) %>% 
  plot_grid(plotlist = ., nrow = 3, labels = "auto", label_size = 8)

ggsave2(filename = "results/plots/sample_input_output_plot-2.pdf", 
        plot = sample_input_output_plot, height = 90, width = 90, units = "mm")

## Summarize by genomes

genomes_summary_plot <-
  list(list(moreno_genomes, "PrettyCols::Autumn"), 
       list(han_genomes, "ggthemes::Classic_Blue_Red_12"),
       list(buck_genomes, "tvthemes::CrystalGems")
       ) %>%
  lapply(function(x) summarize_genomes_plot(x[[1]], x[[2]])) %>% 
  plot_grid(plotlist = ., nrow = 3, labels = "auto", label_size = 8)

ggsave2(filename = "results/plots/genomes_summary_plot-2.pdf", 
        plot = genomes_summary_plot, height = 90, width = 90, units = "mm")


if (!dir.exists("results/tables")) dir.create("results/tables")

output_by_sample_table %>% 
  as_gt() %>% 
  gt::gtsave(filename = "results/tables/summary_by_sample-2.html")

# additional analyses -----------------------------------------------------

# buck_samples %>%
#   left_join(buck_sample_info %>% select(sample_id = run_accession, scientific_name), by = "sample_id") %>% 
#   mutate(virus_detection = case_when(n_genomes > 0 & n_present > 0 ~ "recovered and detected",
#                                      n_genomes > 0 & n_present == 0 ~ "recovered, none detected",
#                                      n_genomes == 0 & n_present > 0 ~ "none recovered, some detected",
#                                      n_genomes == 0 & n_present == 0 ~ "none recovered or detected")) %>%  
#   mutate(input_type = case_when(str_detect(scientific_name, "metagenome") ~ scientific_name,
#                                 TRUE ~ "SAG")) %>% 
#   ggplot(aes(x = total_reads, y = N50, color = virus_detection)) +
#   geom_point() +
#   facet_wrap(~input_type, nrow = 1) +
#   scale_x_log10() +
#   scale_y_log10()

# buck_samples %>%
#   left_join(buck_sample_info %>% select(sample_id = run_accession, scientific_name), by = "sample_id") %>% 
#   mutate(virus_detection = case_when(n_genomes > 0 & n_present > 0 ~ "recovered and detected",
#                                      n_genomes > 0 & n_present == 0 ~ "recovered, none detected",
#                                      n_genomes == 0 & n_present > 0 ~ "none recovered, some detected",
#                                      n_genomes == 0 & n_present == 0 ~ "none recovered or detected")) %>%  
#   mutate(input_type = case_when(str_detect(scientific_name, "metagenome") ~ scientific_name,
#                                 TRUE ~ "SAG")) %>% 
#   # filter(virus_detection %in% "recovered and detected",
#   #        input_type %in% "SAG") %>% 
#   # select(sample_id) %>% 
#   left_join(buck_abundance %>% 
#               filter(vOTU %in% "vOTU06874") %>% 
#               pivot_longer(-vOTU, values_to = "abundance", names_to = "sample_id"), by = "sample_id") %>% 
#   mutate(phix_detected = abundance > 0) %>% 
#   count(phix_detected, virus_detection, input_type)

# buck_samples %>%
#   left_join(buck_sample_info %>% select(sample_id = run_accession, scientific_name), by = "sample_id") %>% 
#   mutate(input_type = case_when(str_detect(scientific_name, "metagenome") ~ scientific_name,
#                               TRUE ~ "SAG")) %>% 
#   filter(!str_detect(scientific_name, "metagenome")) %>% 
#   count(n_present, scientific_name) %>% 
#   pivot_wider(names_from = n_present, values_from = n, values_fill = 0)

# buck_samples %>%
#   left_join(buck_sample_info %>% select(sample_id = run_accession, scientific_name), by = "sample_id") %>% 
#   mutate(input_type = case_when(str_detect(scientific_name, "metagenome") ~ scientific_name,
#                                 TRUE ~ "SAG")) %>% 
#   filter(!str_detect(scientific_name, "metagenome")) %>% 
#   count(n_genomes, scientific_name) %>% 
#   pivot_wider(names_from = n_genomes, values_from = n, values_fill = 0)

# buck_abundance %>%
#   pivot_longer(-vOTU, names_to = "sample_id", values_to = "abundance") %>% 
#   left_join(buck_sample_info %>% select(sample_id = run_accession, scientific_name), by = "sample_id") %>% 
#   mutate(input_type = case_when(str_detect(scientific_name, "metagenome") ~ scientific_name,
#                                 TRUE ~ "SAG")) %>% 
#   filter(!str_detect(scientific_name, "metagenome")) %>% 
#   filter(abundance > 0) %>% 
#   count(vOTU) %>% 
#   left_join(buck_genomes, by = "vOTU")

# buck_sample_with_potential_phix <- 
#   buck_genomes %>% 
#   filter(str_detect(Phylum, "Phix")) %>% 
#   filter(vOTU == "vOTU06874")

# buck_abundance %>%
#   pivot_longer(-vOTU, names_to = "sample_id", values_to = "abundance") %>% 
#   left_join(buck_sample_info %>% select(sample_id = run_accession, scientific_name), by = "sample_id") %>% 
#   mutate(input_type = case_when(str_detect(scientific_name, "metagenome") ~ scientific_name,
#                                 TRUE ~ "SAG")) %>% 
#   filter(!str_detect(scientific_name, "metagenome")) %>% 
#   filter(abundance > 0) %>% 
#   filter(vOTU != "vOTU06874") %>% 
#   left_join(buck_genomes, by = "vOTU") %>% 
#   left_join(buck_sample_info %>% select(sample_id = run_accession, everything())) 

# buck_samples %>% 
#   mutate(potential_isolates = total_reads < 3e7) %>% 
#   filter(potential_isolates) 
