#' Calculate Nodule Progression Detection Score (NPDS)
#'
#' @description
#' The `NPDS_calculate` function calculates the Nodule Progression Detection Score (NPDS) based on the baseline 
#' and follow-up CT scans provided in the `nodule_progress_detector` list. This score quantifies the likelihood 
#' of nodule progression by analyzing lung tissue blocks and comparing Hounsfield Unit (HU) ratios between the 
#' scans. The function utilizes lung tissue blocks, a nodule block list, and a detection threshold sequence to 
#' compute the NPDS value.
#'
#' @param nodule_progress_detector A list returned by the `initialization` function, containing the following fields:
#' \describe{
#'   \item{\code{bf_sub_image}}{The baseline CT subregion containing the nodule.}
#'   \item{\code{af_sub_image}}{The follow-up CT subregion containing the nodule.}
#'   \item{\code{voxel_coord}}{The voxel coordinates of the nodule (\code{c(x, y, z)}).}
#'   \item{\code{split_size}}{The size of the region of interest for lung tissue blocks.}
#'   \item{\code{image_size}}{The width and height of the extracted subregion.}
#' }
#'
#' @return A modified version of the input list, with the following added field:
#' \describe{
#'   \item{\code{NPDS}}{The calculated Nodule Progression Detection Score (NPDS), which reflects the likelihood 
#'   of nodule progression. A positive score indicates progression, while a negative score suggests regression.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Generates lung tissue blocks for both baseline and follow-up subregions using 
#'         \code{generate_lung_tissue_blocks}.
#'   \item Creates a nodule block list by calling \code{generate_nodule_block_list}, focusing on the region of interest 
#'         around the nodule.
#'   \item Detects nodule progression using \code{HU_ratio_nodule_progression_detection}, which compares Hounsfield 
#'         Unit (HU) ratios between corresponding blocks from the baseline and follow-up scans.
#'   \item Computes the detection matrix and detection list based on a range of detection thresholds (\code{detection_lambda}).
#'   \item Integrates the detection results over the threshold range using the trapezoidal rule to calculate 
#'         the NPDS values.
#'   \item Determines the final NPDS value by comparing the mean of positive and negative NPDS values:
#'         \enumerate{
#'           \item If the mean of positive NPDS values is larger in magnitude, the maximum NPDS value is selected.
#'           \item Otherwise, the minimum NPDS value is selected.
#'         }
#' }
#'
#' @examples
#' # Click “Run Example” and wait patiently, as it takes some time to execute.
#' # Alternatively, you can type the following command in the console:
#' # example("NPDS_calculate", local = TRUE)
#' 
#' # Create a simulated nodule_progress_detector object with random 3D arrays
#' nodule_progress_detector <- list(
#'   bf_sub_image = array(rnorm(10 * 512 * 512, mean = -500, sd = 200), dim = c(10, 512, 512)),
#'   af_sub_image = array(rnorm(10 * 512 * 512, mean = -500, sd = 200), dim = c(10, 512, 512)),
#'   voxel_coord = c(256, 256, 5),
#'   split_size = 32,
#'   image_size = 512
#' )
#' cat("Simulated nodule_progress_detector object initialized successfully.\n")
#'
#' # Define placeholder functions for dependency simulation
#' generate_lung_tissue_blocks <- function(image, split_size, image_size) {
#'   # Simulate lung tissue blocks as random 3D arrays
#'   array(rnorm(split_size * split_size * dim(image)[1]), dim = c(dim(image)[1], split_size, split_size))
#' }
#'
#' generate_nodule_block_list <- function(bf_image, af_image, x, y, split_size) {
#'   # Simulate nodule block list as random 3D arrays
#'   list(baseline_blocks = array(rnorm(split_size * split_size), dim = c(split_size, split_size)),
#'        followup_blocks = array(rnorm(split_size * split_size), dim = c(split_size, split_size)))
#' }
#'
#' HU_ratio_nodule_progression_detection <- function(A1, A2, nodule_block_list, split_size, detection_threshold) {
#'   # Simulate detection results
#'   detection_matrix <- matrix(rnorm(length(detection_threshold) * 10), nrow = 10)
#'   detection_list <- matrix(runif(length(detection_threshold) * 10), nrow = 10)
#'   list(detection_matrix = detection_matrix, detection_list = detection_list)
#' }
#' cat("Placeholder functions defined successfully.\n")
#'
#' # Calculate the NPDS score
#' result <- NPDS_calculate(nodule_progress_detector)
#' cat("NPDS calculation completed successfully.\n")
#'
#' # Access and print the calculated NPDS value
#' npds_score <- result$NPDS
#' cat(sprintf("NPDS Score: %.4f\n", npds_score))
#'
#' @importFrom pracma trapz
#' @export
NPDS_calculate <- function(nodule_progress_detector){
  detection_lambda <- seq(1, 100) / 100.0
  #source("R/generate_lung_tissue_blocks.R")
  A1 <- generate_lung_tissue_blocks(nodule_progress_detector$bf_sub_image,
                                    split_size=nodule_progress_detector$split_size,
                                    image_size=nodule_progress_detector$image_size)
  A2 <- generate_lung_tissue_blocks(nodule_progress_detector$af_sub_image,
                                    split_size=nodule_progress_detector$split_size,
                                    image_size=nodule_progress_detector$image_size)
  
  #source("R/generate_nodule_block_list.R")
  nodule_block_list = generate_nodule_block_list(nodule_progress_detector$bf_sub_image,
                                                 nodule_progress_detector$af_sub_image,
                                                 nodule_progress_detector$voxel_coord[1],
                                                 nodule_progress_detector$voxel_coord[2],
                                                 split_size=nodule_progress_detector$split_size)
  
  #source("R/HU_ratio_nodule_progression_detection.R")
  detection <- HU_ratio_nodule_progression_detection(A1,
                                                     A2,
                                                     nodule_block_list=nodule_block_list,
                                                     split_size=nodule_progress_detector$split_size,
                                                     detection_threshold=detection_lambda)
  discrete_stat_matrix <- detection$detection_matrix
  NPDSt_lambda_list <- detection$detection_list
  NPDSt <- array(0,dim = dim(NPDSt_lambda_list)[1])
  for (i in 1:dim(NPDSt_lambda_list)[1]) {
    detection_value <- NPDSt_lambda_list[i, ]
    NPDSt[i] <- pracma::trapz(detection_lambda, detection_value)
  }
  
  max_NPDSt = max(NPDSt)
  min_NPDSt = min(NPDSt)
  
  pos_NPDSt_values <- NPDSt[NPDSt > 0]
  neg_NPDSt_values <- NPDSt[NPDSt < 0]
  
  if(length(pos_NPDSt_values) > 0){
    mean_pos_NPDSt <- mean(pos_NPDSt_values)
  } else {
    mean_pos_NPDSt <- 0
  }
  
  if(length(neg_NPDSt_values) > 0){
    mean_neg_NPDSt <- mean(neg_NPDSt_values)
  } else {
    mean_neg_NPDSt <- 0
  }
  
  if(abs(mean_pos_NPDSt) > abs(mean_neg_NPDSt)){
    nodule_progress_detector$NPDS <- max_NPDSt
  } else {
    nodule_progress_detector$NPDS <- min_NPDSt
  }
  
  return(nodule_progress_detector)
}