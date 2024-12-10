#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector generate_nodule_block_list_slice_cpp_v2(
    NumericMatrix image_slice,
    NumericMatrix image_reg_slice,
    int x_start, int x_end,
    int y_start, int y_end,
    int split_size) {
  NumericVector current_slice(2 * split_size * split_size);
  int pixel_index = 0;
  for (int k = y_start; k < y_end; k++) {
    for (int l = x_start; l < x_end; l++) {
      current_slice[pixel_index] = image_slice(k, l); // 基线CT
      current_slice[split_size * split_size + pixel_index] = image_reg_slice(k, l); // 随访CT
      pixel_index++;
    }
  }
  return current_slice;
}
