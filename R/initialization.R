#' Initialization of Nodule Imaging Data
#'
#' @description
#' This function initializes the imaging data for analyzing a nodule from baseline and follow-up CT scans. 
#' It processes the input parameters to extract and prepare the required imaging data for further analysis, 
#' including voxel coordinates, spatial dimensions, and sub-images around the nodule region.
#'
#' @param X The X coordinate of the center of the nodule in the slice from the follow-up CT scan where the maximum diameter of the nodule reaches its largest value.
#' @param Y The Y coordinate of the center of the nodule in the slice from the follow-up CT scan where the maximum diameter of the nodule reaches its largest value.
#' @param range_Z A string specifying the range of Z slices (e.g., '325-347') that include the nodule in follow-up CT scan.
#' @param diameter The largest value of the maximum diameter of the nodule across slices in the follow-up CT scan.
#' @param baseline_CT_nii_path File path to the baseline CT scan in `.nii` format.
#' @param followup_CT_nii_path File path to the follow-up CT scan in `.nii` format.
#' 
#' @return A list containing:
#' \describe{
#'
#'   \item{\strong{**Coordinate Information**}}{}
#'   \item{\code{coord_x}}{The processed nodule’s x-coordinate value for subsequent analysis.}
#'   \item{\code{coord_y}}{The processed nodule’s y-coordinate value for subsequent analysis.}
#'   \item{\code{range_z}}{The input z-coordinate range.}
#'   \item{\code{coord_z}}{The processed nodule’s z-coordinate value for subsequent analysis.}
#'   \item{\code{voxel_coord}}{The processed nodule’s voxel coordinates for subsequent analysis.}
#'
#'   \item{\strong{**Image Data**}}{}
#'   \item{\code{bf_CT_nii}}{Baseline CT image as a NIfTI object.}
#'   \item{\code{af_CT_nii}}{Follow-up CT image as a NIfTI object.}
#'   \item{\code{bf_CT_npy}}{Baseline CT image data array.}
#'   \item{\code{af_CT_npy}}{Follow-up CT image data array.}
#'   \item{\code{af_spacing}}{Spacing information of the follow-up CT image.}
#'
#'   \item{\strong{**Nodule Dimensions**}}{}
#'   \item{\code{diameter_pixel}}{The maximum value of the maximum diameter of the nodule in pixels.}
#'   \item{\code{diameter_mm}}{The maximum value of the maximum diameter of the nodule in millimeters.}
#'   \item{\code{diameter_z}}{Z-axis extent of the nodule in slices.}
#'
#'   \item{\strong{**Analysis Parameters**}}{}
#'   \item{\code{isflip}}{Logical indicating whether the image coordinates are flipped.}
#'   \item{\code{ClinvNod_NPDS_95th_percentiles}}{Predefined percentile values for subsequent analysis.}
#'   \item{\code{z_start}, \code{z_end}}{The starting and ending slice indices for the processed Z-axis range of the nodules.}
#'   \item{\code{split_size}}{Size of the region of interest for analysis.}
#'
#'   \item{\strong{**Processed Sub-Image**}}{}
#'   \item{\code{af_sub_image}}{Extracted sub-image from the follow-up CT data.}
#'   \item{\code{image_size}}{Width and height of the extracted sub-image.}
#' }
#' 
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Processes the input parameters, including the nodule's spatial coordinates and size.
#'   \item Reads the baseline and follow-up CT scans from the provided file paths using the \code{oro.nifti} package.
#'   \item Computes the Z-axis slice range for the nodule based on the input \code{range_Z}.
#'   \item Extracts spatial information, including spacing and voxel dimensions, from the CT images.
#'   \item Identifies whether the image coordinates need to be flipped, based on spacing values.
#'   \item Extracts a subregion of the follow-up CT image that contains the nodule, defined by the Z-axis range.
#'   \item Returns a list containing the processed information, including voxel coordinates, transposed image data, 
#'         and extracted subregions, to facilitate further analysis.
#' }
#' 
#' @examples
#' # Click “Run Example” and wait patiently, as it takes some time to execute.
#' # It is highly recommended to type "example("initialization", local = TRUE)"
#' # in the console for a more interactive and enhanced experience.
#' 
#' # Run the initialization function with example parameters
#' cat("Initializing nodule imaging data with example parameters...\n")
#' nodule_progress_detector <- initialization(
#'   X = 209,
#'   Y = 356,
#'   range_Z = '325-347',
#'   diameter = 12, # unit: mm
#'   baseline_CT_nii_path = system.file("extdata", "0002358111-20180516.nii.gz", package = "NPDS4Clib"),
#'   followup_CT_nii_path = system.file("extdata", "0002358111-20220707.nii.gz", package = "NPDS4Clib")
#' )
#'
#' # Validate dimensions of key components
#' cat("Baseline CT dimensions (bf_CT_npy):", dim(nodule_progress_detector$bf_CT_npy), "\n")
#' cat("Follow-up CT dimensions (af_CT_npy):", dim(nodule_progress_detector$af_CT_npy), "\n")
#' cat("Sub-image dimensions (af_sub_image):", dim(nodule_progress_detector$af_sub_image), "\n")
#'
#' # Check voxel coordinates and range
#' cat("Voxel coordinates:", nodule_progress_detector$voxel_coord, "\n")
#' cat("Z-axis range:", nodule_progress_detector$range_z, "\n")
#' cat("Z start and end indices:", nodule_progress_detector$z_start, nodule_progress_detector$z_end, "\n")
#'
#' # Check split size and image size
#' cat("Split size:", nodule_progress_detector$split_size, "\n")
#' cat("Image size:", nodule_progress_detector$image_size, "\n")
#'
#' # Verify spacing and diameter calculations
#' cat("Spacing of follow-up CT:", nodule_progress_detector$af_spacing, "\n")
#' cat("Diameter in pixels:", nodule_progress_detector$diameter_pixel, "\n")
#' cat("Diameter in millimeters:", nodule_progress_detector$diameter_mm, "\n")
#'
#' # Check for coordinate flipping
#' cat("Is flipped:", nodule_progress_detector$isflip, "\n")
#'
#' # Visualization
#' slice <- nodule_progress_detector$af_sub_image[1,,]
#' par(mfrow = c(1, 1))
#' image(slice, main = "Extracted Sub-Image Slice", col = gray.colors(256))
#' cat("Visualization complete.\n")
#' 
#' @import oro.nifti
#' @export
initialization <- function(X, Y, range_Z, diameter, baseline_CT_nii_path, followup_CT_nii_path) {
  # Load the oro.nifti package
  #library(oro.nifti)
  # Process the range_Z
  coord_x = X
  coord_y = Y
  range_z <- strsplit(range_Z, "-")[[1]]
  # Read baseline and follow-up CT images
  bf_CT_nii <- oro.nifti::readNIfTI(baseline_CT_nii_path, reorient = FALSE)
  af_CT_nii <- oro.nifti::readNIfTI(followup_CT_nii_path, reorient = FALSE)
  bf_CT_npy = aperm(bf_CT_nii@.Data, c(3, 2, 1))
  af_CT_npy = aperm(af_CT_nii@.Data, c(3, 2, 1))
  af_spacing = af_CT_nii@pixdim[2:4]
  diameter_pixel = diameter/af_CT_nii@pixdim[2]
  isflip = all(af_spacing[2:4] < 0)
  ClinvNod_NPDS_95th_percentiles = c(0.0011799599609374932, 0.005169005859374974, 0.0505342207031249, 0.10536974414062492)
  z_end = dim(af_CT_nii@.Data)[1] - as.integer(range_z[1])
  z_start = dim(af_CT_nii@.Data)[1] - as.integer(range_z[2])
  coord_z = as.integer((z_start + z_end) / 2.0)
  voxel_coord = c(coord_x,coord_y,coord_z)
  diameter_z = as.integer(as.integer(range_z[2]) - as.integer(range_z[1]))
  
  if (diameter_z > 32 || diameter_pixel > 32) {
    split_size <- 64
  } else {
    split_size <- 32
  }
  
  if (isflip) {
    voxel_coord <- c(512 - voxel_coord[1], 512 - voxel_coord[2], voxel_coord[3])
  }
  
  af_sub_image <- af_CT_npy[(z_start + 1):(z_end + 1), , ]
  image_size <- dim(af_sub_image)[2]
  
  # Print size information and initialization completion message
  cat("baseline CT size:", dim(bf_CT_npy), "\n")
  cat("follow-up CT size:", dim(af_CT_npy), "\n")
  cat("Initialization complete.\n")
  
  # Return a list of the coordinates and images in the specified order
  return(list(
    coord_x = X,
    coord_y = Y,
    range_z = range_z,
    bf_CT_nii = bf_CT_nii,
    af_CT_nii = af_CT_nii,
    bf_CT_npy = bf_CT_npy,
    af_CT_npy = af_CT_npy,
    af_spacing = af_spacing,
    diameter_pixel = diameter_pixel,
    diameter_mm = diameter,
    isflip = isflip,
    ClinvNod_NPDS_95th_percentiles = ClinvNod_NPDS_95th_percentiles,
    z_end = z_end,
    z_start = z_start,
    coord_z = coord_z,
    voxel_coord = voxel_coord,
    diameter_z = diameter_z,
    split_size = split_size,
    af_sub_image = af_sub_image,
    image_size = image_size
  )
  )
}