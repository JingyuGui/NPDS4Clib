#include <Rcpp.h>
using namespace Rcpp;

// 防御性编程方式的有效区域过滤函数
std::vector<int> filter_lung_regions(const std::vector<int>& valid_regions, const std::vector<int>& region_areas) {
  // 检查输入参数是否为空
  if (valid_regions.empty() || region_areas.empty()) {
    Rcpp::stop("filter_lung_regions: Empty input vectors detected.");
  }

  // 检查区域标签的合法性
  for (int label : valid_regions) {
    if (label <= 0 || label > static_cast<int>(region_areas.size())) {
      Rcpp::stop("filter_lung_regions: Invalid label found in valid_regions.");
    }
  }

  // 过滤有效区域，找出前两个最大的
  std::vector<int> top_two_labels;
  int max_area_1 = -1;
  int max_area_2 = -1;
  int max_label_1 = -1;
  int max_label_2 = -1;

  for (int label : valid_regions) {
    int area = region_areas[label - 1];  // 获取标签对应的区域面积

    // 判断是否为最大或者第二大的区域
    if (area > max_area_1) {
      // 更新最大和第二大
      max_area_2 = max_area_1;
      max_label_2 = max_label_1;

      max_area_1 = area;
      max_label_1 = label;
    } else if (area > max_area_2) {
      max_area_2 = area;
      max_label_2 = label;
    }
  }

  // 确保找到了两个有效区域
  if (max_label_1 != -1) {
    top_two_labels.push_back(max_label_1);
  }
  if (max_label_2 != -1) {
    top_two_labels.push_back(max_label_2);
  }

  // 如果最终的有效区域不足两个，给出警告
  if (top_two_labels.size() < 2) {
    Rcpp::warning("filter_lung_regions: Less than two valid lung regions were found.");
  }

  return top_two_labels;
}

// [[Rcpp::export]]
List process_lung_regions(IntegerMatrix label_image, IntegerMatrix regions) {
  int n_labels = regions.nrow();
  if (n_labels == 0) {
    Rcpp::stop("process_lung_regions: No regions provided.");
  }

  std::vector<int> region_areas(n_labels, 0);

  // 统计每个区域的实际面积
  for (int i = 0; i < label_image.nrow(); i++) {
    for (int j = 0; j < label_image.ncol(); j++) {
      int label = label_image(i, j);
      if (label > 0) {
        if (label > n_labels) {
          Rcpp::stop("process_lung_regions: Detected a label greater than available regions.");
        }
        region_areas[label - 1] += 1; // label 是从1开始的，但索引从0开始
      }
    }
  }

  // 筛选出有效的肺区域
  std::vector<int> valid_regions;
  for (int i = 0; i < n_labels; i++) {
    int x_min = regions(i, 0);
    int x_max = regions(i, 1);
    int y_min = regions(i, 2);
    int y_max = regions(i, 3);
    int x_range = x_max - x_min;
    int y_range = y_max - y_min;

    if (x_range < 350 && y_range < 350) {
      valid_regions.push_back(i + 1); // 有效区域标签，+1表示区域标记
    }
  }

  // 判断有效区域数量并过滤最大两个有效区域
  std::vector<int> top_two_labels;
  if (valid_regions.size() > 2) {
    top_two_labels = filter_lung_regions(valid_regions, region_areas);
  } else {
    top_two_labels = valid_regions;
  }

  // 遍历图像，过滤掉非最大两个区域，将其设置为背景
  for (int i = 0; i < label_image.nrow(); i++) {
    for (int j = 0; j < label_image.ncol(); j++) {
      int label = label_image(i, j);
      if (label > 0 && std::find(top_two_labels.begin(), top_two_labels.end(), label) == top_two_labels.end()) {
        label_image(i, j) = 0;  // 设置为背景
      }
    }
  }

  // 返回处理后的图像和有效区域
  return List::create(
    _["processed_image"] = label_image,
    _["valid_regions"] = top_two_labels
  );
}
