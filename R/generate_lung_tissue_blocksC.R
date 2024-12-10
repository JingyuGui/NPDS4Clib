#' @keywords internal
generate_lung_tissue_blocksC <- function(image, split_size = 32, image_size = 512){

  # 创建空的3D数组 A
  A <- array(0, dim = c(dim(image)[1], floor(image_size / split_size)^2, split_size^2))

  # 遍历图像的每一张切片
  for (m in 1:dim(image)[1]) {
    A[m, , ] <- generate_lung_tissue_blocks_slice_cpp(image[m, ,], image_size, split_size)
  }

  return(A)
}
