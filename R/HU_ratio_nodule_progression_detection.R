#' @keywords internal
HU_ratio_nodule_progression_detection <- function(
    A1,#generate_lung_tissue_blocks调用基线图像的return
    A2,#generate_lung_tissue_blocks调用随访图像的return
    anno_i=NULL,
    anno_j=NULL,
    nodule_block_list=NULL,#nodule_block_list的return
    split_size=32,
    image_size=512,
    detection_threshold=c(0.1)){
  M <- dim(A1)[1]
  R <- length(detection_threshold)
  split_num <- floor(image_size / split_size)
  block_num <- split_num^2

  change_ratio_matrix <- array(0, dim=c(M,block_num))
  detection_matrix <- array(0, dim=c(R,M,block_num))
  detection_list <- array(0, dim=c(M,R))
  for (m in 1:M){
    if (!is.null(anno_i) && !is.null(anno_j)) {
      nodule_block_1 <- A1[m, anno_i * split_num + anno_j, ]
      nodule_block_2 <- A2[m, anno_i * split_num + anno_j, ]
    } else {
      nodule_block_1 <- nodule_block_list[m, 1, ]
      nodule_block_2 <- nodule_block_list[m, 2, ]
    }
    for (i in 1:split_num){
      for (j in 1:split_num){
        block_1d_index = (i-1) * split_num + j
        block_1_now = A1[m, block_1d_index, ]
        block_2_now = A2[m, block_1d_index, ]
        ratio_1 = nodule_block_1 / abs(block_1_now + 0.1)
        ratio_2 = nodule_block_2 / abs(block_2_now + 0.1)
        change = (mean(ratio_2) - mean(ratio_1)) / abs(mean(ratio_1))
        change_ratio_matrix[m, block_1d_index] <- change
        for (r in 1:R){
          if(change > detection_threshold[r]){
            detection_matrix[r, m, block_1d_index] <- 1.0
          }
          else if(change < -detection_threshold[r]){
            detection_matrix[r, m, block_1d_index] <- -1.0
          }
          else{
            detection_matrix[r, m, block_1d_index] <- 0.0
          }

        }
      }
    }
    for (r in 1:R){
      detection_list[m, r] = sum(detection_matrix[r, m, ]) / block_num
    }
  }
  detection <- list(detection_matrix = detection_matrix,
                    detection_list = detection_list)
  return(detection)
}
