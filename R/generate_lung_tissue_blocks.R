#' @keywords internal
generate_lung_tissue_blocks <- function(image, split_size = 32, image_size = 512){

  # 创建空的3D数组 A
  A <- array(0, dim = c(dim(image)[1], floor(image_size / split_size)^2, split_size^2))

  # 遍历图像的每一张切片
  for (m in 1:dim(image)[1]) {
    for (i in 1:floor(image_size / split_size)) {
      for (j in 1:floor(image_size / split_size)) {

        # 计算当前的 split_index
        split_index <- (i - 1) * floor(image_size / split_size) + (j - 1) + 1

        # 提取出当前的 32x32 的子矩阵
        sub_matrix <- image[m, ((i - 1) * split_size + 1):(i * split_size), ((j - 1) * split_size + 1):(j * split_size)]

        # 将子矩阵赋值到 A 中，展平子矩阵为一维
        for (k in 1:split_size) {
          for (l in 1:split_size) {
            pixel_index <- (k - 1) * split_size + l
            A[m, split_index, pixel_index] <- sub_matrix[k, l]
          }
        }
      }
    }
  }

  return(A)
}


