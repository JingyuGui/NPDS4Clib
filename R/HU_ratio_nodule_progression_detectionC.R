#' @keywords internal
HU_ratio_nodule_progression_detectionC <- function(
    A1,#generate_lung_tissue_blocks调用基线图像的return
    A2,#generate_lung_tissue_blocks调用随访图像的return
    anno_i=NULL,
    anno_j=NULL,
    nodule_block_list=NULL,#nodule_block_list的return
    split_size=32,
    image_size=512,
    detection_threshold){
  M <- dim(A1)[1]
  R <- length(detection_threshold)
  split_num <- floor(image_size / split_size)
  block_num <- split_num^2

  #change_ratio_matrix <- array(0, dim=c(M,block_num))
  detection_matrix <- array(0, dim=c(R,M,block_num))
  detection_list <- array(0, dim=c(M,R))

  for (m in 1:M){
    detection_loop <- HU_ratio_nodule_progression_detection_slice_cpp(
      A1[m, , ],#Numeric Matrix
      A2[m, , ],#Numeric Matrix
      anno_i,
      anno_j,
      nodule_block_list[m, , ],
      split_num,
      split_num^2,
      detection_threshold)

    detection_matrix[ ,m, ] <- detection_loop$detection_matrix_slice
    detection_list[m, ] <- detection_loop$detection_list_slice
  }
  detection <- list(detection_matrix = detection_matrix,
                    detection_list = detection_list)
  return(detection)
}
