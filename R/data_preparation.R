#' Prepare Data
#'
#' @param data A data frame containing neuroimaging measurements
#' @param normalize_by_icv Logical, whether to normalize by intracranial volume
#' @param scale_volumes Logical, whether to scale the volumes
#' @param metadata_cols Character vector of metadata column names
#' @return A data frame in long format with processed measurements
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
  
  # Find common regions
  fs_regions <- names(fs_data)[colSums(!is.na(fs_data)) > 0]
  synth_regions <- names(synth_data)[colSums(!is.na(synth_data)) > 0]
  common_regions <- intersect(fs_regions, synth_regions)
  common_regions <- setdiff(common_regions, metadata_cols)
  
  # Keep only common regions and metadata
  keep_cols <- intersect(names(data), c(metadata_cols, common_regions))
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
    processed_data[common_regions] <- lapply(processed_data[common_regions], scale)
  }
  
  # Convert to long format
  processed_data <- processed_data %>%
    tidyr::pivot_longer(
      cols = all_of(common_regions),
      names_to = "region",
      values_to = "volume"
    )
  
  return(processed_data)
}