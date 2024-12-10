#' Perform Rigid Registration of Baseline and Follow-Up CT Scans
#'
#' @description
#' The `registration_by_elastix` function performs rigid registration to align baseline and follow-up CT scans.
#' Using the `RNiftyReg` package, it applies a rigid transformation to the baseline CT image so that it aligns 
#' with the follow-up CT image, allowing for comparative analysis of nodules over time. The function expects 
#' as input a list returned by the `initialization` function and updates this list with the registration 
#' results, including transformed images and extracted subregions.
#'
#' @param input A list returned by the `initialization` function. This list must include the following elements:
#' \describe{
#'   \item{\code{bf_CT_nii}}{The baseline CT image in NIfTI format.}
#'   \item{\code{af_CT_nii}}{The follow-up CT image in NIfTI format.}
#'   \item{\code{z_start}}{The starting slice index for the Z-axis range of the nodule.}
#'   \item{\code{z_end}}{The ending slice index for the Z-axis range of the nodule.}
#' }
#'
#' @return A modified list with the following updated elements:
#' \describe{
#'   \item{\code{bf_CT_npy}}{The registered baseline CT image as a 3D data array, transposed for analysis.}
#'   \item{\code{bf_sub_image}}{A subregion of the registered baseline CT image extracted based on the Z-axis range.}
#'   \item{\code{registration_result}}{The result object returned by `RNiftyReg` containing details of the registration.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that the input is a list and contains all required elements.
#'   \item Uses \code{RNiftyReg::niftyreg} to perform rigid registration of the baseline CT image 
#'         (\code{bf_CT_nii}) to align it with the follow-up CT image (\code{af_CT_nii}).
#'   \item Transposes the registered baseline image to match the expected array structure.
#'   \item Extracts a subregion of the registered baseline image based on the Z-axis range 
#'         (\code{z_start} and \code{z_end}).
#'   \item Updates the input list with the registered image, extracted subregion, and registration results.
#' }
#'
#' @examples
#' # Click “Run Example” and wait oh-so-patiently! As it takes some time (10 min maybe) to execute.
#' # It is highly recommended to type "example("registration_by_elastix", local = TRUE)"
#' # in the console for a more interactive and enhanced experience.
#' 
#' # Initialize the data
#' cat("Initializing nodule imaging data...\n")
#' nodule_data <- initialization(
#'   X = 209,
#'   Y = 356,
#'   range_Z = '325-347',
#'   diameter = 12, # unit: mm
#'   baseline_CT_nii_path = system.file("extdata", "0002358111-20180516.nii.gz", package = "NPDS4Clib"),
#'   followup_CT_nii_path = system.file("extdata", "0002358111-20220707.nii.gz", package = "NPDS4Clib")
#' )
#'
#' # Perform rigid registration
#' cat("Performing rigid registration...\n")
#' registered_data <- registration_by_elastix(nodule_data)
#' cat("Registration complete.\n")
#'
#' # Validate dimensions of registered image and subregion
#' cat("Dimensions of registered baseline CT (bf_CT_npy):", dim(registered_data$bf_CT_npy), "\n")
#' cat("Dimensions of registered subregion (bf_sub_image):", dim(registered_data$bf_sub_image), "\n")
#'
#' # Check registration results
#' if (!is.null(registered_data$registration_result)) {
#'   cat("Registration result is available.\n")
#' } else {
#'   stop("Registration result is missing!")
#' }
#'
#' # Validate that the subregion is correctly extracted
#' z_start <- registered_data$z_start
#' z_end <- registered_data$z_end
#' expected_subregion_size <- z_end - z_start + 1
#' actual_subregion_size <- dim(registered_data$bf_sub_image)[1]
#' if (expected_subregion_size == actual_subregion_size) {
#'   cat("Subregion extraction size is correct.\n")
#' } else {
#'   stop("Subregion extraction size mismatch!")
#' }
#'
#' # Visualize the first slice of the registered subregion
#' slice <- registered_data$bf_sub_image[1, , ]
#' image(slice, main = "First Slice of Registered Subregion", col = gray.colors(256))
#' cat("Visualization complete.\n")
#'
#' @import oro.nifti
#' @import RNiftyReg
#' @export
registration_by_elastix <- function(input) {
  
  # Check if the input is a list
  if (!is.list(input)) {
    stop("Input must be a list returned by the initialization function.")
  }
  
  # Extract necessary elements from the input list
  bf_CT_nii <- input$bf_CT_nii
  af_CT_nii <- input$af_CT_nii
  z_start <- input$z_start
  z_end <- input$z_end
  
  # Perform rigid registration using RNiftyReg with baseline CT as the source image
  registration_result <- RNiftyReg::niftyreg(
    source = bf_CT_nii,   # Baseline CT image to be transformed
    target = af_CT_nii,   # Follow-up CT image as the target
    scope = "rigid",
    interpolation = 0
  )
  
  # Obtain the registered baseline image
  bf_CT_npy <- aperm(registration_result$image, c(3, 2, 1))
  bf_sub_image <- bf_CT_npy[(z_start + 1):(z_end + 1), , ]
  
  # Update elements in the input list
  input$bf_CT_npy <- bf_CT_npy
  input$bf_sub_image <- bf_sub_image
  
  # Print size information
  message("baseline CT size: ", paste(dim(input$bf_CT_npy), collapse = " x "))
  message("follow-up CT size: ", paste(dim(as.array(input$af_CT_npy)), collapse = " x "))
  message("Registration complete.")
  
  input$registration_result <- registration_result
  # Return the updated list
  return(input)
}