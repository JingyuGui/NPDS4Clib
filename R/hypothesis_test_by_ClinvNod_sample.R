#' Perform Hypothesis Testing for Nodule Progression Prediction
#'
#' @description
#' The `hypothesis_test_by_ClinvNod_sample` function performs a hypothesis test to assess the progression likelihood 
#' of a nodule based on the Nodule Progression Detection Score (NPDS). It compares the calculated NPDS against 
#' pre-defined percentile thresholds for different nodule size groups, and evaluates its statistical significance 
#' using sample distributions from reference clinical data.
#'
#' @param nodule_progress_detector A list containing the necessary parameters for hypothesis testing, including:
#' \describe{
#'   \item{\code{diameter_mm}}{The maximum diameter of the nodule in millimeters.}
#'   \item{\code{NPDS}}{The calculated Nodule Progression Detection Score (NPDS).}
#'   \item{\code{ClinvNod_NPDS_95th_percentiles}}{A vector of pre-defined 95th percentile thresholds for different 
#'   nodule size groups (e.g., <=5 mm, <=10 mm, <=15 mm, >15 mm).}
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{NPDS}}{The calculated NPDS value.}
#'   \item{\code{Progression}}{A logical value indicating whether the nodule progression is predicted (\code{TRUE}) 
#'   or not (\code{FALSE}).}
#'   \item{\code{p_value}}{The p-value indicating the statistical significance of the NPDS compared to the reference 
#'   clinical sample distribution.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Loads the diameter, NPDS, and 95th percentile thresholds from the input list.
#'   \item Determines the nodule size group based on \code{diameter_mm}:
#'         \enumerate{
#'           \item Group 1: Diameter <= 5 mm.
#'           \item Group 2: Diameter <= 10 mm.
#'           \item Group 3: Diameter <= 15 mm.
#'           \item Group 4: Diameter > 15 mm.
#'         }
#'   \item Compares the NPDS value against the corresponding 95th percentile threshold for the determined group.
#'   \item Reads the reference clinical sample data for the determined group from a CSV file.
#'   \item Calculates the p-value as the proportion of sample NPDS values greater than the calculated NPDS.
#' }
#'
#' @examples
#' # Click “Run Example” and wait patiently, as it takes some time to execute.
#' # Alternatively, you can type the following command in the console:
#' # example("hypothesis_test_by_ClinvNod_sample", local = TRUE)
#' 
#' # Simulate a nodule_progress_detector object
#' nodule_progress_detector <- list(
#'   diameter_mm = 12,
#'   NPDS = 0.45,
#'   ClinvNod_NPDS_95th_percentiles = c(0.30, 0.35, 0.40, 0.50)
#' )
#' cat("Simulation completed. Nodule parameters initialized.\n")
#'
#' # Mock the CSV file reading function for demonstration purposes
#' read.csv <- function(file_path, stringsAsFactors = FALSE) {
#'   data.frame(S = runif(100, min = 0.2, max = 0.6)) # Simulated reference sample data
#' }
#' cat("Mock function defined successfully.\n")
#'
#' # Perform the hypothesis test
#' result <- hypothesis_test_by_ClinvNod_sample(nodule_progress_detector)
#' cat("Hypothesis test completed successfully.\n")
#'
#' # Access and print the result with explanations
#' cat("Results of the hypothesis test:\n")
#' cat("Calculated NPDS:", result$NPDS, "\n")
#' cat("Progression Prediction Result:", ifelse(result$Progression, "Progression Detected", "No Progression"), "\n")
#' cat("Calculated p-value:", result$p_value, "\n")
#'
#' # Provide interpretation of the results
#' if (result$p_value < 0.05) {
#'   cat("The result is statistically significant (p-value < 0.05).\n")
#' } else {
#'   cat("The result is not statistically significant (p-value >= 0.05).\n")
#' }
#'
#' @seealso \code{\link{NPDS_calculate}}
#' @importFrom utils read.csv
#' @export
hypothesis_test_by_ClinvNod_sample <- function(nodule_progress_detector) {
  # Load parameters
  diameter_mm <- nodule_progress_detector$diameter_mm
  NPDS <- nodule_progress_detector$NPDS
  ClinvNod_NPDS_95th_percentiles <- nodule_progress_detector$ClinvNod_NPDS_95th_percentiles
  
  # Determine the group
  if (diameter_mm <= 5) {
    group_reference_num <- 1
    progress <- NPDS > ClinvNod_NPDS_95th_percentiles[1]
  } else if (diameter_mm <= 10) {
    group_reference_num <- 2
    progress <- NPDS > ClinvNod_NPDS_95th_percentiles[2]
  } else if (diameter_mm <= 15) {
    group_reference_num <- 3
    progress <- NPDS > ClinvNod_NPDS_95th_percentiles[3]
  } else {
    group_reference_num <- 4
    progress <- NPDS > ClinvNod_NPDS_95th_percentiles[4]
  }
  
  # Construct the file path
  file_name <- paste0("ClinvSample_NPDS_G", group_reference_num, ".csv")
  file_path <- system.file("extdata", file_name, package = "NPDS4Clib")
  
  # Read the CSV file
  ClinvNod_NPDSs <- read.csv(file_path, stringsAsFactors = FALSE)$S
  
  # Calculate significance
  ClinvNod_sig <- sapply(ClinvNod_NPDSs, function(x) ifelse(x > NPDS, 1, 0))
  p_value <- sum(ClinvNod_sig) / length(ClinvNod_sig)
  
  # Print the result
  cat(sprintf("NPDS: %.10f\nProgression Prediction Result: %s\np_value: %.10f\n",
              NPDS, progress, p_value))
  
  # Return the result as a list
  return(list(NPDS = NPDS, Progression = progress, p_value = p_value))
}