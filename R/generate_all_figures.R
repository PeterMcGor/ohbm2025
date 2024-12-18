# Load libraries and helper functions
library(ggplot2)
source("R/data_preparation.R")
source("R/analysis.R")
source("R/subject_analysis.R")
source("app.R")



# Load and prepare data
data <- read.csv("data/synthseg_FS_concat_all_versions.csv")
processed_data <- prepare_data(data, normalize_by_icv = FALSE, scale_volumes = FALSE)
subject_analysis <- analyze_subject_differences(processed_data)

# Output folder for figures
output_dir <- "figures/"
dir.create(output_dir, showWarnings = FALSE)

# Input values
analysis_types <- c("software", "cv", "distances", "ratios")
regions <- unique(processed_data$region) # For 'software'
region_types <- c("lh", "rh", "other")   # For 'cv', 'distances', 'ratios'
distance_metrics <- c("mean_js_div", "mean_kendall_tau", "mean_wilcox_paired")
ratio_metrics <- c("js_ratio", "kendall_ratio", "wilcox_paired_ratio")

# Generate and save plots for each analysis type
for (analysis in analysis_types) {
  
  if (analysis == "software") {
    # Software comparison for each region
    for (region in regions) {
      cat("Generating software comparison plot for region:", region, "\n")
      results <- analyze_region(processed_data, region, p_adjust_method = "bonferroni")
      plot <- results[[4]]  # Assuming this returns the ggplot
      ggsave(filename = paste0(output_dir, "software_", region, ".pdf"), 
             plot = plot, width = 12, height = 16)
    }
    
  } else if (analysis == "cv") {
    # Coefficient of Variation plots
    for (region_type in region_types) {
      cat("Generating CV plot for region type:", region_type, "\n")
      plot <- create_single_metric_plot(subject_analysis$results, 
                                        "mean_within_cv",
                                        "Within-Subject Coefficient of Variation", 
                                        region_type)
      ggsave(filename = paste0(output_dir, "cv_", region_type, ".pdf"), 
             plot = plot, width = 12, height = 16)
    }
    
  } else if (analysis == "distances") {
    # Distance metrics plots
    for (metric in distance_metrics) {
      for (region_type in region_types) {
        cat("Generating distance metric plot:", metric, "for", region_type, "\n")
        plot <- create_single_metric_plot(subject_analysis$results, 
                                          metric,
                                          get_metric_title(metric), 
                                          region_type)
        ggsave(filename = paste0(output_dir, "distances_", metric, "_", region_type, ".pdf"), 
               plot = plot, width = 12, height = 16)
      }
    }
    
  } else if (analysis == "ratios") {
    # Discrimination ratio plots
    for (metric in ratio_metrics) {
      for (region_type in region_types) {
        cat("Generating ratio plot:", metric, "for", region_type, "\n")
        plot <- create_single_metric_plot(subject_analysis$results, 
                                          metric,
                                          get_metric_title(metric), 
                                          region_type)
        ggsave(filename = paste0(output_dir, "ratios_", metric, "_", region_type, ".pdf"), 
               plot = plot, width = 12, height = 16)
      }
    }
  }
}

cat("All figures have been successfully generated!\n")
