#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix generate_lung_tissue_blocks_slice_cpp(NumericMatrix image_slice, int image_size, int split_size) {

  int num_splits = floor(image_size / split_size);  // 切分数
  NumericMatrix result(num_splits * num_splits, split_size * split_size);  // 用来存放拆分后的子块

  // 遍历切片的每个小块
  for (int i = 0; i < num_splits; i++) {
    for (int j = 0; j < num_splits; j++) {

      int split_index = i * num_splits + j;
      NumericMatrix sub_matrix(split_size, split_size);  // 用来存储当前的32x32子矩阵

      // 提取当前小块的子矩阵
      for (int k = 0; k < split_size; k++) {
        for (int l = 0; l < split_size; l++) {
          sub_matrix(k, l) = image_slice(i * split_size + k, j * split_size + l);
        }
      }

      // 将子矩阵展平并赋值到结果矩阵的相应位置
      for (int k = 0; k < split_size; k++) {
        for (int l = 0; l < split_size; l++) {
          int pixel_index = k * split_size + l;
          result(split_index, pixel_index) = sub_matrix(k, l);
        }
      }
    }
  }

  return result;
}
