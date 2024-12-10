#' @keywords internal
generate_nodule_block_listC <- function(
    image,#基线CT bf_sub_image
    image_reg,#随访CT
    X,#结节X坐标
    Y,#结节Y坐标
    split_size =32){
  nodule_block_list <- array(0,
                             dim = c(dim(image)[1],
                                     2,
                                     split_size^2))
  x_start = floor(X - split_size / 2)-1
  x_end = floor(X + split_size / 2)-1
  y_start = floor(Y - split_size / 2)-1
  y_end = floor(Y + split_size / 2)-1 #因为是要传给Rcpp的参数所以下标需要调整

  # 遍历每个切片
  for (m in 1:dim(image)[1]) {# 第一维是切片
    nodule_block_list[m, , ] <- generate_nodule_block_list_slice_cpp(image[m, , ],
                                                                     image_reg[m, , ],
                                                                     x_start,
                                                                     x_end,
                                                                     y_start,
                                                                     y_end,
                                                                     split_size)
  }
  return(nodule_block_list)
}
