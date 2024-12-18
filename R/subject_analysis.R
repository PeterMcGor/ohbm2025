# Load required libraries
library(tidyverse)
library(ggplot2)
library(patchwork)

#' Calculate distribution-based distances with units (mm³)
#' @param data1,data2 Numeric vectors (paired by site)
calculate_unit_distances <- function(data1, data2) {
  emp_wass <- mean(abs(sort(data1) - sort(data2)))
  
  mu1 <- mean(data1)
  mu2 <- mean(data2)
  sigma1 <- sd(data1)
  sigma2 <- sd(data2)
  gauss_wass <- sqrt((mu1 - mu2)^2 + (sigma1 - sigma2)^2)
  
  list(
    empirical_wasserstein = emp_wass,
    gaussian_wasserstein = gauss_wass
  )
}

#' Calculate Jensen-Shannon divergence
#' @param data1,data2 Numeric vectors (paired by site)
calculate_js_divergence <- function(data1, data2) {
  mu1 <- mean(data1)
  mu2 <- mean(data2)
  sigma1 <- sd(data1)
  sigma2 <- sd(data2)
  
  m <- (mu1 + mu2)/2
  s <- sqrt((sigma1^2 + sigma2^2)/4)
  
  0.5 * (
    (log(s/sigma1) + (sigma1^2 + (mu1 - m)^2)/(2 * s^2) - 0.5) +
      (log(s/sigma2) + (sigma2^2 + (mu2 - m)^2)/(2 * s^2) - 0.5)
  )
}

#' Calculate rank-based distances
#' @param data1,data2 Numeric vectors (paired by site)
calculate_rank_distances <- function(data1, data2) {
  n <- length(data1)
  
  r1 <- rank(data1)
  r2 <- rank(data2)
  footrule <- sum(abs(r1 - r2))
  footrule_norm <- footrule / (n * (n-1)/2)
  
  kt <- suppressWarnings(
    cor.test(data1, data2, method="kendall", exact=FALSE)
  )
  tau_dist <- (1 - kt$estimate)/2
  
  calculate_wilcoxon_effect <- function(x, y, paired=TRUE) {
    if(all(x == y)) return(0)
    wx <- wilcox.test(x, y, paired=paired, exact=FALSE, correct=TRUE)
    z <- qnorm(wx$p.value/2)
    return(abs(z)/sqrt(if(paired) n else 2*n))
  }
  
  effect_paired <- calculate_wilcoxon_effect(data1, data2, paired=TRUE)
  effect_unpaired <- calculate_wilcoxon_effect(data1, data2, paired=FALSE)
  
  list(
    footrule = footrule_norm,
    kendall_tau = tau_dist,
    wilcox_effect_paired = effect_paired,
    wilcox_effect_unpaired = effect_unpaired
  )
}

#' Calculate all distances between two samples
#' @param data1,data2 Numeric vectors (must be paired by site)
calculate_all_distances <- function(data1, data2) {
  if(length(data1) != length(data2)) stop("Data vectors must have same length")
  
  unit_dist <- calculate_unit_distances(data1, data2)
  js_dist <- calculate_js_divergence(data1, data2)
  rank_dist <- calculate_rank_distances(data1, data2)
  
  c(unit_dist, js_divergence = js_dist, rank_dist)
}

#' Calculate intra-subject variability metrics
calculate_intra_subject <- function(data) {
  data %>%
    group_by(region, software, sub) %>%
    summarize(
      within_std = sd(volume),
      within_cv = sd(volume) / mean(volume),
      within_mad = mad(volume),
      within_iqr = IQR(volume),
      n_measurements = n(),
      .groups = 'keep'
    ) %>%
    group_by(region, software) %>%
    summarize(
      mean_within_std = mean(within_std, na.rm = TRUE),
      mean_within_cv = mean(within_cv, na.rm = TRUE),
      mean_within_mad = mean(within_mad, na.rm = TRUE),
      mean_within_iqr = mean(within_iqr, na.rm = TRUE),
      .groups = 'keep'
    )
}

#' Calculate pairwise distances between subjects
calculate_subject_distances <- function(data) {
  subjects <- unique(data$sub)
  n_subjects <- length(subjects)
  dists <- list()
  
  for(i in 1:(n_subjects-1)) {
    for(j in (i+1):n_subjects) {
      data_i <- filter(data, sub == subjects[i]) %>% arrange(site)
      data_j <- filter(data, sub == subjects[j]) %>% arrange(site)
      
      if(!all(data_i$site == data_j$site)) {
        stop(paste("Sites don't match for subjects", subjects[i], "and", subjects[j]))
      }
      
      curr_dists <- calculate_all_distances(data_i$volume, data_j$volume)
      
      for(metric in names(curr_dists)) {
        dists[[metric]] <- c(dists[[metric]], curr_dists[[metric]])
      }
    }
  }
  
  lapply(dists, mean)
}

#' Create horizontal barplot
#' @param data Data frame containing the measurements
#' @param metric_col Column name for the metric to plot
#' @param title Plot title
#' @param group_prefix Optional prefix for filtering regions
#' @param sort_by How to sort the regions ("metric" or "alphabetical")
create_horizontal_barplot <- function(data, metric_col, title, group_prefix = NULL, 
                                      sort_by = "metric") {
  data <- subset(data, region != "CSF")
  
  filtered_data <- if (!is.null(group_prefix)) {
    if (is.numeric(group_prefix)) {
      data %>% 
        group_by(region) %>%
        summarise(
          across(all_of(metric_col), max),
          .groups = 'drop'
        ) %>%
        arrange(desc(.data[[metric_col]])) %>%
        slice_head(n = group_prefix) %>%
        pull(region) %>%
        {filter(data, region %in% .)}
    } else if (group_prefix == "other") {
      data %>% 
        filter(!startsWith(region, "lh_") & 
                 !startsWith(region, "rh_") &
                 !startsWith(region, "lhCortex") &
                 !startsWith(region, "rhCortex") &
                 !startsWith(region, "lhCerebral") &
                 !startsWith(region, "rhCerebral"))
    } else {
      data %>% 
        filter(startsWith(region, group_prefix))
    }
  } else {
    data
  }
  
  region_order <- if(sort_by == "metric") {
    filtered_data %>%
      group_by(region) %>%
      summarise(
        metric_value = max(.data[[metric_col]]),
        .groups = 'drop'
      ) %>%
      arrange(metric_value) %>%
      pull(region)
  } else if(sort_by == "alphabetical") {
    sort(unique(filtered_data$region))
  } else {
    unique(filtered_data$region)
  }
  
  ggplot(filtered_data, aes(x = .data[[metric_col]], y = region, fill = software)) +
    geom_col(position = position_dodge(width = 0.9)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    ) +
    scale_y_discrete(limits = region_order) +
    labs(
      title = title,
      x = gsub("_", " ", metric_col),
      y = NULL
    ) +
    scale_fill_manual(values = c(
      "FSv7.4.1" = "#5ea290",
      "FSv8.0.0" = "#0267b1",
      "FSv8.0.0_robust" = "#e9a600",
      "SynthSeg+v2.0" = "#9e2444"
    ))
}

#' Main analysis function
#' @export
analyze_subject_differences <- function(data) {
  intra_subject <- calculate_intra_subject(data)
  inter_subject <- calculate_inter_subject(data)
  
  results <- inner_join(inter_subject, intra_subject, 
                        by = c("region", "software")) %>%
    mutate(
      js_ratio = mean_js_div / mean_within_cv,
      footrule_ratio = mean_footrule / mean_within_cv,
      kendall_ratio = mean_kendall_tau / mean_within_cv,
      wilcox_paired_ratio = mean_wilcox_paired / mean_within_cv,
      wilcox_unpaired_ratio = mean_wilcox_unpaired / mean_within_cv
    )
  
  unit_metrics <- c("mean_empirical_wasserstein", "mean_gaussian_wasserstein")
  dimensionless_metrics <- c("mean_js_div", "mean_footrule", "mean_kendall_tau",
                             "mean_wilcox_paired", "mean_wilcox_unpaired")
  ratio_metrics <- c("js_ratio", "footrule_ratio", "kendall_ratio",
                     "wilcox_paired_ratio", "wilcox_unpaired_ratio")
  
  list(
    results = results,
    unit_metrics = unit_metrics,
    dimensionless_metrics = dimensionless_metrics,
    ratio_metrics = ratio_metrics
  )
}

# Helper functions for plotting in the Shiny app
plot_cv_analysis <- function(analysis_results, region) {
  analysis_results$results %>%
    filter(region == !!region) %>%
    create_horizontal_barplot("mean_within_cv", 
                              paste("Within-subject CV for", region))
}

plot_inter_subject_metrics <- function(analysis_results, region, metric_type) {
  data <- analysis_results$results %>%
    filter(region == !!region)
  
  metrics <- switch(metric_type,
                    "unit" = analysis_results$unit_metrics,
                    "dimensionless" = analysis_results$dimensionless_metrics,
                    "ratio" = analysis_results$ratio_metrics)
  
  plots <- lapply(metrics, function(metric) {
    create_horizontal_barplot(data, metric, 
                              paste(gsub("_", " ", metric), "for", region))
  })
  
  wrap_plots(plots, ncol = 1)
}

plot_summary_statistics <- function(analysis_results, region) {
  data <- analysis_results$results %>%
    filter(region == !!region)
  
  metrics <- c("mean_within_cv", 
               analysis_results$unit_metrics[1],
               analysis_results$dimensionless_metrics[1],
               analysis_results$ratio_metrics[1])
  
  plots <- lapply(metrics, function(metric) {
    create_horizontal_barplot(data, metric, 
                              paste(gsub("_", " ", metric), "for", region))
  })
  
  wrap_plots(plots, ncol = 2)
}