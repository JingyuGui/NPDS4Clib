#' Extract Lung Masks for Baseline and Follow-Up Sub-Images
#'
#' @description
#' The `get_segmented_lungs` function processes 3D sub-images of baseline and follow-up CT scans to extract lung masks. 
#' It iteratively applies the `get_segmented_lungs_in_CT_slice` function to each slice of the sub-images, 
#' segmenting the lung regions. The function updates the input `nodule_progress_detector` list with the processed 
#' sub-images and their corresponding binary lung masks.
#'
#' @param nodule_progress_detector A list returned by the `initialization` function. This list must include:
#' \describe{
#'   \item{\code{bf_sub_image}}{The subregion of the baseline CT scan, extracted based on the nodule's Z-axis range.}
#'   \item{\code{af_sub_image}}{The subregion of the follow-up CT scan, extracted based on the nodule's Z-axis range.}
#' }
#'
#' @return A modified version of the input list, with the following updated or added fields:
#' \describe{
#'   \item{\code{bf_sub_image}}{The processed subregion of the baseline CT scan, with non-lung regions set to zero.}
#'   \item{\code{af_sub_image}}{The processed subregion of the follow-up CT scan, with non-lung regions set to zero.}
#'   \item{\code{bf_sub_binary}}{A binary mask of the baseline CT scan subregion, indicating lung regions (\code{TRUE}) and non-lung regions (\code{FALSE}).}
#'   \item{\code{af_sub_binary}}{A binary mask of the follow-up CT scan subregion, indicating lung regions (\code{TRUE}) and non-lung regions (\code{FALSE}).}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the 3D sub-images (\code{bf_sub_image} and \code{af_sub_image}) from the input list.
#'   \item Iterates over each slice of the baseline sub-image (\code{bf_sub_image}) and:
#'         \enumerate{
#'           \item Segments the lungs in the slice using \code{get_segmented_lungs_in_CT_slice}.
#'           \item Updates the processed sub-image and binary mask for the slice.
#'         }
#'   \item Iterates over each slice of the follow-up sub-image (\code{af_sub_image}) and performs the same operations.
#'   \item Stores the processed sub-images and their binary masks in the input list.
#' }
#'
#' @examples
#' # Click “Run Example” and wait patiently, as it takes some time to execute.
#' # It is highly recommended to type "example("get_segmented_lungs", local = TRUE)"
#' # in the console for a more interactive and enhanced experience.
#' 
#' # Create a simulated nodule_progress_detector object with random 3D arrays
#' nodule_progress_detector <- list(
#'   bf_sub_image = array(rnorm(10 * 512 * 512, mean = -500, sd = 200), dim = c(10, 512, 512)),
#'   af_sub_image = array(rnorm(10 * 512 * 512, mean = -500, sd = 200), dim = c(10, 512, 512))
#' )
#' cat("Simulated nodule_progress_detector object initialized successfully.\n")
#'
#' # Run the lung mask extraction function
#' processed_data <- get_segmented_lungs(nodule_progress_detector)
#'
#' # Access the processed sub-images and binary masks
#' bf_processed <- processed_data$bf_sub_image
#' af_processed <- processed_data$af_sub_image
#' bf_binary_mask <- processed_data$bf_sub_binary
#' af_binary_mask <- processed_data$af_sub_binary
#' cat("Processed sub-images and binary masks retrieved successfully.\n")
#'
#' # Verify the dimensions of the results
#' cat("Baseline processed image dimensions:", dim(bf_processed), "\n")
#' cat("Follow-up processed image dimensions:", dim(af_processed), "\n")
#' cat("Baseline binary mask dimensions:", dim(bf_binary_mask), "\n")
#' cat("Follow-up binary mask dimensions:", dim(af_binary_mask), "\n")
#'
#' # Verify that processed sub-images match the binary masks
#' bf_match <- all((bf_processed == 0) == (bf_binary_mask == 0))
#' af_match <- all((af_processed == 0) == (af_binary_mask == 0))
#' cat("Baseline processed image matches binary mask:", bf_match, "\n")
#' cat("Follow-up processed image matches binary mask:", af_match, "\n")
#'
#' # Calculate the proportion of lung regions in binary masks
#' bf_lung_ratio <- sum(bf_binary_mask) / length(bf_binary_mask)
#' af_lung_ratio <- sum(af_binary_mask) / length(af_binary_mask)
#' cat("Proportion of lung area in baseline binary mask:", bf_lung_ratio, "\n")
#' cat("Proportion of lung area in follow-up binary mask:", af_lung_ratio, "\n")
#'
#' @seealso \code{\link{get_segmented_lungs_in_CT_slice}}, \code{\link{initialization}}
#' @export
get_segmented_lungs <- function(nodule_progress_detector) {
  
  #source("R/get_segmented_lungs_in_CT_slice.R")
  
  bf_sub_image <- nodule_progress_detector$bf_sub_image
  af_sub_image <- nodule_progress_detector$af_sub_image
  bf_sub_binary <- bf_sub_image
  af_sub_binary <- af_sub_image
  message("Running lung mask extraction ...")
  
  # Process each slice of bf_sub_image
  for (i in 1:dim(bf_sub_image)[1]) {
    cat("Processing bf_sub_image slice:", i, "\n")
    bf_sub_image[i,,] <- get_segmented_lungs_in_CT_slice(bf_sub_image[i,,])$im
    bf_sub_binary[i,,] <- get_segmented_lungs_in_CT_slice(bf_sub_image[i,,])$binary
  }
  
  # Process each slice of af_sub_image
  for (i in 1:dim(af_sub_image)[1]) {
    cat("Processing af_sub_image slice:", i, "\n")
    af_sub_image[i,,] <- get_segmented_lungs_in_CT_slice(af_sub_image[i,,])$im
    af_sub_binary[i,,] <- get_segmented_lungs_in_CT_slice(af_sub_image[i,,])$binary
  }
  
  # Save the processed binary images into new fields
  nodule_progress_detector$bf_sub_image <- bf_sub_image
  nodule_progress_detector$af_sub_image <- af_sub_image
  
  nodule_progress_detector$bf_sub_binary <- bf_sub_binary
  nodule_progress_detector$af_sub_binary <- af_sub_binary
  
  message("Lung mask extraction complete.")
  
  return(nodule_progress_detector)
}