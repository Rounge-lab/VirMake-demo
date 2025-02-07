

reorg_family <- function(dat, n = 10) {
  dat %>% 
    mutate(Family = str_replace(Family, "^n.a.$", "No reference"),
           Family = str_replace(Family, "^Unclassified$", "Reference unclassified"),
           Family = fct_infreq(Family),
           Family = fct_other(Family, keep = levels(Family)[1:min(n, length(levels(Family)))]))
}

summarize_samples <- function(data_by_sample, by_var = NULL) {
  
  summarize_vars <- c(`Viral genomes identified` = "n_genomes",
                      `vOTUs mapped` = "n_present",
                      `Percent reads mapped` = "mapped_percent")
  
  tmp <-
    data_by_sample %>% 
    mutate(mapped_percent = str_remove(mapped_percent, "%") %>% as.double()) %>% 
    mutate(fraction_provirus = n_provirus/n_genomes) %>% 
    rename(all_of(summarize_vars))
  if (!is.null(by_var)) tmp <- tmp %>% rename(by_var = all_of(by_var))
  tmp %>% 
    tbl_summary(include = all_of(names(summarize_vars)), by = by_var)
}


summarize_genomes <- function(data_by_genome, by_var = NULL) {
  
  summarize_vars <- c(`Genomes per vOTU` = "mult_genomes", 
                      `Virus integration` = "any_provirus", 
                      `Genome Family` = "Family")
  
  tmp <-
    data_by_genome %>% 
    reorg_family() %>% 
    mutate(mult_genomes = ifelse(n_genomes > 1, "Multiple genomes", "Single genome")) %>% 
    mutate(any_provirus = ifelse(pangenome_n_provirus > 0, 
                                 "Integrated genome detected", 
                                 "No integration")) %>% 
    rename(all_of(summarize_vars))
  
  if (!is.null(by_var)) {
    tmp %>% 
      rename(by_var = all_of(by_var)) %>% 
      tbl_summary(include = all_of(names(summarize_vars)), by = by_var)
    
  } else {
    tmp %>% 
      tbl_summary(include = all_of(names(summarize_vars)))
  }
}



summarize_resource_usage <- function(res_usage_df) {
  res_usage_df %>% 
    pivot_longer(-Dataset) %>% 
    mutate(name = factor(name, 
                         levels = c("Moreno", "Han", "Buck"), 
                         labels = c("Moreno-Gallego", "Li", "Buck"))) %>% 
    mutate(variable = factor(Dataset, 
                             levels = c("Samples", "Dataset size (Gbp)", 
                                        "Memory needed (GB)", "CPU time (hours)", 
                                        "Wall time (hours)"))) %>% 
    tbl_custom_summary(
      by = name,
      label = list(variable ~ "Totals per dataset"),
      stat_fns = everything() ~ function(data, ...) tibble(rep_val = mean(data$value)),
      statistic = everything() ~ "{rep_val}",
      include = variable
    )
}

summarize_genomes_plot <- function(data_by_genome, pal = "nord::frost") {
  
  n_genomes_plot <-
    data_by_genome %>% 
    reorg_family() %>% 
    group_by(Family) %>% 
    summarize(n_reps = n(),
              n_genomes = sum(n_genomes), .groups = "drop") %>% 
    ggplot(aes(y = Family)) +
    geom_col(aes(x = n_genomes), fill = "gray", width = 0.7) +
    geom_col(aes(x = n_reps, fill = Family)) +
    scale_fill_paletteer_d(palette = pal) +
    theme_light() +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, 3),
                  labels = scales::trans_format("log10",scales::math_format(10^.x))) +
    labs(x = "number of\ngenomes")
  
  genome_len_plot <-
    data_by_genome %>% 
    mutate(pangenome_median_length = pangenome_median_length/1e3) %>% 
    reorg_family() %>% 
    ggplot(aes(x = pangenome_median_length, y = Family, fill = Family)) +
    geom_boxplot(outlier.shape = 20, outlier.size = 0.4) +
    scale_fill_paletteer_d(palette = pal) +
    scale_x_log10() +
    theme_light() +
    theme(axis.text.y = element_blank()) +
    labs(x = "genome\nlength (kb)")
  
  int_genome_plot <-
    data_by_genome %>% 
    reorg_family() %>% 
    group_by(Family) %>% 
    summarize(frac_integrated = sum(pangenome_n_provirus)/sum(n_genomes), .groups = "drop") %>% 
    ggplot(aes(x = frac_integrated, y = Family, color = Family)) +
    geom_point() +
    scale_color_paletteer_d(palette = pal) +
    scale_x_continuous(limits = c(0,1)) +
    theme_light() +
    theme(axis.text.y = element_blank()) +
    labs(x = "fraction\nintegrated")
  
  plots <- list(n_genomes_plot, genome_len_plot, int_genome_plot)
  leg <- get_legend(n_genomes_plot)
  combined_plot <- plot_grid(plotlist = lapply(plots, function(x) x + labs(y = "") + 
                                                 theme(legend.position = "none",
                                                       text = element_text(size = 6))), 
                             nrow = 1, rel_widths = c(1.5,1,1))
  combined_plot
}

plot_input_v_output <- function(data_by_sample) {
  
  data_by_sample %>% 
    pivot_longer(c(N50, total_reads), names_to = "input_var", values_to = "input_val") %>% 
    mutate(input_var = factor(input_var, levels = c("total_reads", "N50"), labels = c("Reads", "Assembly N50"))) %>% 
    pivot_longer(c(n_genomes, n_present), names_to = "output_var", values_to = "output_val") %>% 
    mutate(output_var = factor(output_var, levels = c("n_genomes", "n_present"), labels = c("Genomes recovered", "vOTUs detected"))) %>% 
    ggplot(aes(x = input_val, y = output_val)) +
    geom_point(size = 0.5, shape = 20) +
    scale_x_log10() +
    labs(x = "", y = str_wrap("Genomes recovered/ vOTUs detected", 20)) +
    facet_wrap(~input_var+output_var, scales = "free_x", nrow = 1) +
    theme_light() +
    theme(text = element_text(size = 6))
  
  
}



