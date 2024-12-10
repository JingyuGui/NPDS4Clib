generate_nodule_block_listC2 <- function(image, image_reg, X, Y, split_size = 32) {
  n_slices <- dim(image)[1]
  n_pixels <- split_size^2
  nodule_block_list <- numeric(n_slices * 2 * n_pixels)

  x_start <- floor(X - split_size / 2) - 1
  x_end <- floor(X + split_size / 2) - 1
  y_start <- floor(Y - split_size / 2) - 1
  y_end <- floor(Y + split_size / 2) - 1

  for (m in 1:n_slices) {
    current_slice <- generate_nodule_block_list_slice_cpp_v2(
      image[m, , ],
      image_reg[m, , ],
      x_start, x_end, y_start, y_end, split_size
    )
    nodule_block_list[((m-1)*2*n_pixels+1):(m*2*n_pixels)] <- current_slice
  }

  # 原始一维向量
  dim(nodule_block_list) <- c(n_pixels, 2, n_slices) # 按列优先整形
  nodule_block_list <- aperm(nodule_block_list, c(3, 2, 1)) # 调整维度顺序

  return(nodule_block_list)
}
