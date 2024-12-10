#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix generate_nodule_block_list_slice_cpp(
    NumericMatrix image_slice,
    NumericMatrix image_reg_slice,
    int x_start, int x_end,
    int y_start, int y_end,
    int split_size) {

  // 创建一个二维矩阵，用来存储结节块
  NumericMatrix nodule_block_list_slice(2, split_size * split_size);

  // 遍历每个小块的像素并填充nodule_block_list
  int pixel_index = 0;
  for (int k = y_start; k < y_end; k++) {
    for (int l = x_start; l < x_end; l++) {
      // 基线CT和随访CT的像素值分别填充到nodule_block_list
      nodule_block_list_slice(0, pixel_index) = image_slice(k, l);
      nodule_block_list_slice(1, pixel_index) = image_reg_slice(k, l);
      pixel_index++;
    }
  }

  return nodule_block_list_slice;
}



