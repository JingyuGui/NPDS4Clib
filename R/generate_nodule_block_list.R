#' @keywords internal
generate_nodule_block_list <- function(
    image,#基线CT bf_sub_image
    image_reg,#随访CT
    X,#结节X坐标
    Y,#结节Y坐标
    split_size =32){
  nodule_block_list <- array(0,
                             dim = c(dim(image)[1],
                                     2,
                                     split_size^2))
  x_start = floor(X - split_size / 2)
  x_end = floor(X + split_size / 2)
  y_start = floor(Y - split_size / 2)
  y_end = floor(Y + split_size / 2)

  # 遍历每个切片
  for (m in 1:dim(image)[1]) {  # 第一维是切片
    # 提取当前切片对应的两个图像（基线CT和随访CT）的nodule块
    nodule_block_now1 <- image[m, y_start:y_end, x_start:x_end]
    nodule_block_now2 <- image_reg[m, y_start:y_end, x_start:x_end]

    # 遍历每个小块的像素
    for (k in 1:split_size) {  # 第一维是行
      for (l in 1:split_size) {  # 第二维是列
        pixel_index <- (k - 1) * split_size + l  # 计算当前像素的索引
        # 将像素值赋给nodule_block_list
        nodule_block_list[m, 1, pixel_index] <- nodule_block_now1[k, l]  # 基线CT
        nodule_block_list[m, 2, pixel_index] <- nodule_block_now2[k, l]  # 随访CT
      }
    }
  }
  return(nodule_block_list)
}


