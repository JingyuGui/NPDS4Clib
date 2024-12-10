#' Segment Lungs in a CT Slice
#'
#' @description
#' The `get_segmented_lungs_in_CT_slice` function processes a single CT slice to segment the lung regions.
#'
#' @param im A 2D numeric matrix representing a CT slice. Pixel intensity values are expected to be in Hounsfield Units (HU).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{im}}{The processed CT slice with non-lung regions set to zero.}
#'   \item{\code{binary}}{A binary mask where \code{TRUE} indicates lung regions and \code{FALSE} indicates non-lung regions.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Thresholds the CT slice to create a binary image, where pixels with values less than -400 HU are considered potential lung regions.
#'   \item Applies the \code{clear_border} function to remove objects touching the border of the image.
#'   \item Performs connected component labeling using \code{bwlabel} to identify distinct regions in the binary image.
#'   \item Extracts bounding box properties of labeled regions using \code{regionprops_bbox}.
#'   \item Processes the identified lung regions using \code{process_lung_regions} to remove small or irrelevant regions.
#'   \item Updates the binary mask to reflect the cleaned lung regions.
#'   \item Sets all non-lung regions in the original CT slice to zero, based on the updated binary mask.
#' }
#'
#' @examples
#' # Click “Run Example” and wait patiently, as it takes some time to execute.
#' # Alternatively, you can type the following command in the console:
#' # example("get_segmented_lungs_in_CT_slice", local = TRUE)
#' 
#' # Create a simulated CT slice with intensity values in Hounsfield Units (HU)
#' simulated_ct_slice <- matrix(rnorm(512 * 512, mean = -500, sd = 200), nrow = 512, ncol = 512)
#' cat("Simulated CT slice created successfully. Dimensions:", dim(simulated_ct_slice), "\n")
#'
#' # Segment the lungs
#' cat("Starting lung segmentation process...\n")
#' segmented_result <- get_segmented_lungs_in_CT_slice(simulated_ct_slice)
#' cat("Lung segmentation completed successfully.\n")
#'
#' # Access the processed CT slice and binary mask
#' processed_ct_slice <- segmented_result$im
#' binary_mask <- segmented_result$binary
#' cat("Processed CT slice and binary mask extracted successfully.\n")
#'
#' # Verify the results
#' # Check dimensions
#' cat("Dimensions of simulated_ct_slice:", dim(simulated_ct_slice), "\n")
#' cat("Dimensions of processed_ct_slice:", dim(processed_ct_slice), "\n")
#' cat("Dimensions of binary_mask:", dim(binary_mask), "\n")
#'
#' # Ensure binary_mask contains only 0 and 1
#' unique_values <- unique(binary_mask)
#' if (all(unique_values %in% c(0, 1))) {
#'   cat("Binary mask verification passed: contains only 0 and 1.\n")
#' } else {
#'   warning("Binary mask verification failed: contains unexpected values!\n")
#' }
#'
#' @useDynLib NPDS4Clib
#' @import Rcpp
#' @export
get_segmented_lungs_in_CT_slice <- function(im) {
  # Threshold the CT image to create a binary mask
  binary <- im < -400
  
  # Call the compiled clear_border function from the shared library
  cleared <- .Call("_NPDS4Clib_clear_border", binary, 0L, 0.0, PACKAGE = "NPDS4Clib")
  
  # Label connected regions in the cleared binary mask
  label_image <- bwlabel(cleared)
  
  # Call the compiled regionprops_bbox function to extract region properties
  regions <- .Call("_NPDS4Clib_regionprops_bbox", label_image, PACKAGE = "NPDS4Clib")
  
  # Call the compiled process_lung_regions function to clean and process lung regions
  cleaned_image <- .Call("_NPDS4Clib_process_lung_regions", label_image$labeled_image, regions, PACKAGE = "NPDS4Clib")
  
  # Update the binary mask with the processed lung regions
  binary <- cleaned_image$processed_image != 0
  
  # Set non-lung regions in the original image to 0
  im[binary == 0] <- 0
  
  # Return the processed image and the binary lung mask
  return(list(im = im, binary = binary))
}