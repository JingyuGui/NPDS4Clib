#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List HU_ratio_nodule_progression_detection_slice_cpp(
    NumericMatrix A1_slice,  // 基线图像的切片
    NumericMatrix A2_slice,  // 随访图像的切片
    Nullable<int> anno_i,    // 可选的i坐标
    Nullable<int> anno_j,    // 可选的j坐标
    NumericMatrix nodule_block_list_slice, // 结节块列表
    int split_num,  // 每行每列的分块数
    int block_num,  // 总块数
    NumericVector detection_threshold // 阈值
) {
  // 动态计算 R
  int R = detection_threshold.size();

  // 初始化返回值矩阵
  NumericMatrix detection_matrix_slice(R, split_num * split_num);  // 检测矩阵
  NumericVector detection_list_slice(R);      // 检测列表
  NumericVector change_ratio_matrix_slice(split_num * split_num); // 变化率矩阵

  NumericVector nodule_block_1, nodule_block_2;

  // 判断是否有指定的 i 和 j 坐标
  if (anno_i.isNotNull() && anno_j.isNotNull()) {
    int i = as<int>(anno_i);
    int j = as<int>(anno_j);

    // 计算切片的索引
    int index = (i - 1) * split_num;

    nodule_block_1 = A1_slice.row(index);  // 赋值而非重新声明
    nodule_block_2 = A2_slice.row(index);  // 赋值而非重新声明
  } else {
    nodule_block_1 = nodule_block_list_slice.row(0);  // 赋值
    nodule_block_2 = nodule_block_list_slice.row(1);  // 赋值
  }

  // 用于调试的计数器
  int count = 0;

  // 遍历所有的分块
  for (int i = 0; i < split_num; i++) {
    for (int j = 0; j < split_num; j++) {
      // 计算一维索引
      int block_1d_index = i * split_num + j;

      // 检查 block_1d_index 是否超出 block_num 范围
      if (block_1d_index >= block_num) {
        Rcpp::Rcout << "Error: block_1d_index exceeds block_num: "
                    << block_1d_index << std::endl;
        stop("block_1d_index exceeds block_num.");
      }

      // 打印调试信息，查看索引是否正确
      //Rcpp::Rcout << "Processing block: " << count++ << ", block_1d_index: "
                  //<< block_1d_index << std::endl;

      // 从 A1 和 A2 中提取当前块
      NumericVector block_1_now = A1_slice.row(block_1d_index);
      NumericVector block_2_now = A2_slice.row(block_1d_index);

      // 计算比例
      NumericVector ratio_1 = nodule_block_1 / abs(block_1_now + 0.1);
      NumericVector ratio_2 = nodule_block_2 / abs(block_2_now + 0.1);

      // 计算变化率
      double mean_ratio_1 = mean(ratio_1);
      double mean_ratio_2 = mean(ratio_2);
      double change = (mean_ratio_2 - mean_ratio_1) / abs(mean_ratio_1);

      // 更新 change_ratio_matrix
      change_ratio_matrix_slice[block_1d_index] = change;

      // 根据变化率更新 detection_matrix
      for (int r = 0; r < R; r++) {
        if (change > detection_threshold[r]) {
          detection_matrix_slice(r, block_1d_index) = 1.0;
        } else if (change < -detection_threshold[r]) {
          detection_matrix_slice(r, block_1d_index) = -1.0;
        } else {
          detection_matrix_slice(r, block_1d_index) = 0.0;
        }
      }
    }
  }

  // 计算每个块的检测值
  for (int r = 0; r < R; r++) {
    detection_list_slice[r] = sum(detection_matrix_slice(r,_)) / block_num;  // 每个块的平均检测值
  }





  // 返回检测矩阵和检测列表
  return List::create(Named("detection_matrix_slice") = detection_matrix_slice,
                      Named("detection_list_slice") = detection_list_slice);
}
