# =================================================================
# Analysis functions for neuroimaging software comparison
# Author: [Your Name]
# =================================================================

#' Prepare Data for Analysis
#' 
#' @param data Input data frame with brain region measurements
#' @param normalize_by_icv Boolean, whether to normalize by intracranial volume
#' @param scale_volumes Boolean, whether to scale the volumes
#' @param metadata_cols Character vector of metadata column names to exclude
#' @return Data frame in long format with processed measurements
#' @export
prepare_data <- function(data, 
                         normalize_by_icv = FALSE, 
                         scale_volumes = FALSE, 
                         metadata_cols = c('id', 'software', 'sub', 'site', 
                                           'contrast', 'ses', 'sequence', 
                                           'Unnamed: 0')) {
  
  # Filter to keep only regions present in both software types
  fs_data <- data[data$software == 'FSv7.4.1', ]
  synth_data <- data[data$software == 'SynthSeg+v2.0', ]
  
  fs_regions <- names(fs_data)[colSums(!is.na(fs_data)) > 0]
  synth_regions <- names(synth_data)[colSums(!is.na(synth_data)) > 0]
  common_regions <- intersect(fs_regions, synth_regions)
  common_regions <- setdiff(common_regions, metadata_cols)
  
  # Keep only common regions and metadata
  keep_cols = intersect(names(data), c(metadata_cols, common_regions))
  processed_data <- data[, keep_cols]
  
  # Normalize by ICV if requested
  if (normalize_by_icv && "EstimatedTotalIntraCranialVol" %in% names(processed_data)) {
    icv_col <- processed_data$EstimatedTotalIntraCranialVol
    region_cols <- setdiff(common_regions, "EstimatedTotalIntraCranialVol")
    processed_data[region_cols] <- lapply(processed_data[region_cols], 
                                          function(x) x / icv_col)
  }
  
  # Scale if requested 
  if (scale_volumes) {
    processed_data[region_cols] <- lapply(processed_data[common_regions], scale)
  }
  
  # Convert to long format
  processed_data %>%
    pivot_longer(
      cols = all_of(common_regions),
      names_to = "region",
      values_to = "volume"
    )
}

#' Analyze Region
#'
#' Performs comprehensive analysis of volumetric data for a specific brain 
#' region, including statistical tests, correlations, and visualizations.
#'
#' @param data Processed data frame in long format
#' @param region_name Name of the region to analyze
#' @param pairwise_test Statistical test function (default: wilcox.test)
#' @param p_adjust_method Method for p-value adjustment
#' @return List containing:
#'         - pairwise_tests: Results of statistical comparisons
#'         - correlations: Correlation analysis results
#'         - mixed_model: Summary of mixed effects model
#'         - final_plot: Combined visualization
#' @export
analyze_region <- function(data, 
                           region_name, 
                           pairwise_test = wilcox.test, 
                           p_adjust_method = "bonferroni") {
  
  # Filter data for the specified region
  region_data <- data %>%
    filter(region == region_name)
  
  # Generate software pairs for comparison
  software_list <- as.character(unique(region_data$software))
  software_pairs <- expand.grid(
    software1 = software_list,
    software2 = software_list,
    stringsAsFactors = FALSE
  ) %>%
    filter(software1 != software2) %>%
    filter(!duplicated(t(apply(., 1, sort))))
  
  # Prepare data for comparisons
  comparable_region_data <- region_data %>%
    select(sub, site, software, volume) %>%
    pivot_wider(
      names_from = software,
      values_from = volume
    )
  
  # Perform statistical tests
  pairwise_results <- perform_pairwise_tests(comparable_region_data, 
                                             software_pairs, 
                                             pairwise_test)
  pairwise_results$p_adjusted <- p.adjust(pairwise_results$p_value, 
                                          method = p_adjust_method)
  
  # Fit mixed effects model
  region_model <- lmer(volume ~ software + (1|sub) + (1|site), 
                       data = region_data)
  
  # Calculate correlations
  correlations <- calculate_correlations(comparable_region_data, software_pairs)
  
  # Create visualizations
  final_plot <- create_visualizations(region_data, comparable_region_data,
                                      pairwise_results, software_pairs,
                                      correlations, region_name)
  
  # Return results
  list(
    pairwise_tests = pairwise_results,
    correlations = correlations,
    mixed_model = summary(region_model),
    final_plot = final_plot
  )
}

# Helper Functions ======================================================

#' Perform Pairwise Statistical Tests
#' @keywords internal
perform_pairwise_tests <- function(comparable_data, software_pairs, pairwise_test) {
  map_dfr(1:nrow(software_pairs), function(i) {
    software1 <- software_pairs$software1[i]
    software2 <- software_pairs$software2[i]
    
    test <- tryCatch({
      pairwise_test(
        comparable_data[[software1]], 
        comparable_data[[software2]], 
        paired = TRUE
      )
    }, warning = function(w) {
      pairwise_test(
        comparable_data[[software1]], 
        comparable_data[[software2]], 
        paired = TRUE,
        exact = FALSE
      )
    })
    
    data.frame(
      software1 = software1,
      software2 = software2,
      comparison = paste(software1, "vs", software2),
      statistic = if("statistic" %in% names(test)) test$statistic else NA,
      p_value = test$p.value
    )
  })
}

#' Calculate Correlations Between Software Pairs
#' @keywords internal
calculate_correlations <- function(comparable_data, software_pairs) {
  map_dfr(1:nrow(software_pairs), function(i) {
    software1 <- software_pairs$software1[i]
    software2 <- software_pairs$software2[i]
    
    cor_test <- cor.test(comparable_data[[software1]], 
                         comparable_data[[software2]], 
                         method = "pearson")
    
    data.frame(
      software1 = software1,
      software2 = software2,
      correlation = cor_test$estimate,
      cor_pvalue = cor_test$p.value
    )
  })
}

#' Create Visualizations
#' @keywords internal
create_visualizations <- function(region_data, comparable_data,
                                  pairwise_results, software_pairs,
                                  correlations, region_name) {
  
  # Setup for plots
  software_levels <- levels(factor(region_data$software))
  max_y <- max(region_data$volume, na.rm = TRUE)
  y_range <- max_y - min(region_data$volume, na.rm = TRUE)
  shape_values <- c(49:57)
  
  # Create significance annotations
  sig_annotations <- create_significance_annotations(pairwise_results, 
                                                     software_levels,
                                                     max_y, y_range)
  
  # Create box plot
  box_plot <- create_box_plot(region_data, sig_annotations, 
                              shape_values, y_range, region_name)
  
  # Create paired plots and Bland-Altman plots
  paired_plots <- create_paired_plots(comparable_data, software_pairs,
                                      correlations, shape_values)
  bland_altman_plots <- create_bland_altman_plots(comparable_data,
                                                  software_pairs,
                                                  shape_values)
  
  # Combine plots
  n_comparisons <- nrow(software_pairs)
  shared_legend <- get_legend(box_plot)
  box_plot <- box_plot + theme(legend.position = "none")
  
  box_plot / 
    shared_legend /
    wrap_plots(paired_plots, ncol = n_comparisons) /
    wrap_plots(bland_altman_plots, ncol = n_comparisons) +
    plot_layout(heights = c(1, 0.2, 
                            n_comparisons/floor(n_comparisons), 
                            n_comparisons/floor(n_comparisons)))
}


#' Create Significance Annotations
#' @keywords internal
create_significance_annotations <- function(pairwise_results, software_levels, max_y, y_range) {
  pairwise_results %>%
    mutate(
      y_position = seq(max_y + 0.1*y_range, 
                       max_y + (0.1 + 0.15*nrow(pairwise_results))*y_range, 
                       length.out = nrow(pairwise_results)),
      significance = sprintf("p = %.2e", p_adjusted),
      xmin = match(software1, software_levels),
      xmax = match(software2, software_levels),
      sig_color = case_when(
        p_adjusted < 0.001 ~ "#f50253",
        p_adjusted < 0.01 ~ "#8000b0",
        p_adjusted < 0.05 ~ "#5501ff",
        TRUE ~ "black"
      )
    )
}

#' Create Box Plot
#' @keywords internal
create_box_plot <- function(region_data, sig_annotations, shape_values, y_range, region_name) {
  # Define jitter based on number of subjects
  gj <- if (length(unique(region_data$sub)) < 10) {
    geom_jitter(aes(color = site, shape = sub), 
                size = 3, 
                position = position_jitter(width = 0.2))
  } else {
    geom_jitter(aes(color = site), 
                size = 3, 
                position = position_jitter(width = 0.2))
  }
  
  # Create base plot
  box_plot <- ggplot(region_data, aes(x = software, y = volume)) +
    geom_boxplot() +
    gj +
    scale_shape_manual(values = shape_values) +
    labs(title = paste("Volume Comparison for", region_name),
         y = "Volume") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "horizontal")
  
  # Add significance bars
  box_plot +
    geom_segment(data = sig_annotations,
                 aes(x = xmin, xend = xmax,
                     y = y_position, yend = y_position),
                 color = sig_annotations$sig_color) +
    geom_segment(data = sig_annotations,
                 aes(x = xmin, xend = xmin,
                     y = y_position - 0.03*y_range, yend = y_position),
                 color = sig_annotations$sig_color) +
    geom_segment(data = sig_annotations,
                 aes(x = xmax, xend = xmax,
                     y = y_position - 0.03*y_range, yend = y_position),
                 color = sig_annotations$sig_color) +
    geom_text(data = sig_annotations,
              aes(x = (xmin + xmax)/2, 
                  y = y_position + 0.05*y_range,
                  label = significance),
              size = 3,
              color = sig_annotations$sig_color)
}

#' Create Paired Plots
#' @keywords internal
create_paired_plots <- function(comparable_data, software_pairs, correlations, shape_values) {
  map(1:nrow(software_pairs), function(i) {
    software1 <- software_pairs$software1[i]
    software2 <- software_pairs$software2[i]
    
    # Get correlation for this pair
    cor_value <- correlations$correlation[i]
    
    # Calculate axis limits
    x_min <- min(comparable_data[[software1]], na.rm = TRUE)
    x_max <- max(comparable_data[[software1]], na.rm = TRUE)
    y_max <- max(comparable_data[[software2]], na.rm = TRUE)
    
    ggplot(comparable_data, 
           aes(x = .data[[software1]], y = .data[[software2]])) +
      geom_point(aes(shape = sub, color = site), size = 3) +
      scale_shape_manual(values = shape_values) +
      geom_abline(intercept = 0, slope = 1, 
                  linetype = "dashed", color = "red") +
      annotate("text",
               x = x_min + (x_max - x_min)*0.3,
               y = y_max,
               label = sprintf("italic(R) == %.3f", cor_value),
               parse = TRUE) +
      theme_minimal() +
      theme(legend.position = "none")
  })
}

#' Create Bland-Altman Plots
#' @keywords internal
create_bland_altman_plots <- function(comparable_data, software_pairs, shape_values) {
  map(1:nrow(software_pairs), function(i) {
    software1 <- software_pairs$software1[i]
    software2 <- software_pairs$software2[i]
    
    # Calculate differences and means
    measurements1 <- comparable_data[[software1]]
    measurements2 <- comparable_data[[software2]]
    mean_measurements <- (measurements1 + measurements2) / 2
    diff_measurements <- measurements1 - measurements2
    
    # Calculate statistics
    mean_diff <- mean(diff_measurements)
    sd_diff <- sd(diff_measurements)
    upper_loa <- mean_diff + 1.96 * sd_diff
    lower_loa <- mean_diff - 1.96 * sd_diff
    
    # Create plot data
    plot_data <- data.frame(
      mean = mean_measurements,
      difference = diff_measurements,
      sub = comparable_data$sub,
      site = comparable_data$site
    )
    
    # Create plot
    ggplot(plot_data, aes(x = mean, y = difference)) +
      geom_point(aes(shape = sub, color = site), size = 3) +
      scale_shape_manual(values = shape_values) +
      geom_hline(yintercept = mean_diff, 
                 linetype = "dashed", color = "red") +
      geom_hline(yintercept = upper_loa, 
                 linetype = "dotted", color = "blue") +
      geom_hline(yintercept = lower_loa, 
                 linetype = "dotted", color = "blue") +
      labs(x = "Mean", 
           y = sprintf("%s - %s", software1, software2)) +
      theme_minimal() +
      theme(legend.position = "none")
  })
}