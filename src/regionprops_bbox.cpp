#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix regionprops_bbox(List input) {
  IntegerMatrix labeled_image = input["labeled_image"];
  int num_labels = input["num_labels"];

  // 初始化一个二维矩阵来存储每个标签的边界框
  IntegerMatrix bboxes(num_labels, 4);

  // 遍历每个标签（排除背景标签 0）
  for (int label = 1; label <= num_labels; label++) {
    int x_min = labeled_image.nrow();
    int x_max = -1;
    int y_min = labeled_image.ncol();
    int y_max = -1;

    // 遍历图像的所有像素
    for (int i = 0; i < labeled_image.nrow(); i++) {
      for (int j = 0; j < labeled_image.ncol(); j++) {
        if (labeled_image(i, j) == label) {
          if (i < x_min) x_min = i;
          if (i > x_max) x_max = i;
          if (j < y_min) y_min = j;
          if (j > y_max) y_max = j;
        }
      }
    }

    // 将边界框信息存储到矩阵中
    if (x_max != -1) {
      bboxes(label - 1, 0) = x_min;
      bboxes(label - 1, 1) = x_max;
      bboxes(label - 1, 2) = y_min;
      bboxes(label - 1, 3) = y_max;
    } else {
      // 如果标签没有有效像素点，用 NA 表示
      bboxes(label - 1, 0) = NA_INTEGER;
      bboxes(label - 1, 1) = NA_INTEGER;
      bboxes(label - 1, 2) = NA_INTEGER;
      bboxes(label - 1, 3) = NA_INTEGER;
    }
  }

  return bboxes;
}

