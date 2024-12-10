#include <Rcpp.h>
#include <queue>
#include <stack>
#include <cmath>
#include <unordered_set>

using namespace Rcpp;

struct XYPoint {
  int x, y;
  XYPoint() {}
  XYPoint(int xx, int yy) : x(xx), y(yy) {}
};

// flood fill 模板函数
template <class T>
void _floodFill(T *m, XYPoint size, XYPoint xy, T rc, double tol = 1e-3) {
  std::stack<XYPoint> s;
  XYPoint pt = xy;
  bool spanLeft, spanRight;
  T tc = m[pt.x + pt.y * size.x];

  if (std::fabs(static_cast<double>(tc - rc)) <= tol) {
    rc = static_cast<T>(rc + tol + 1);
  }

  s.push(pt);

  while (!s.empty()) {
    pt = s.top();
    s.pop();

    while (pt.y >= 0 && std::fabs(static_cast<double>(m[pt.x + pt.y * size.x] - tc)) <= tol) {
      pt.y--;
    }
    pt.y++;
    spanLeft = spanRight = false;

    while (pt.y < size.y && std::fabs(static_cast<double>(m[pt.x + pt.y * size.x] - tc)) <= tol) {
      m[pt.x + pt.y * size.x] = rc;

      if (!spanLeft && pt.x > 0 && std::fabs(static_cast<double>(m[(pt.x - 1) + pt.y * size.x] - tc)) <= tol) {
        s.push(XYPoint(pt.x - 1, pt.y));
        spanLeft = true;
      } else if (spanLeft && pt.x > 0 && std::fabs(static_cast<double>(m[(pt.x - 1) + pt.y * size.x] - tc)) > tol) {
        spanLeft = false;
      }

      if (!spanRight && pt.x < size.x - 1 && std::fabs(static_cast<double>(m[(pt.x + 1) + pt.y * size.x] - tc)) <= tol) {
        s.push(XYPoint(pt.x + 1, pt.y));
        spanRight = true;
      } else if (spanRight && pt.x < size.x - 1 && std::fabs(static_cast<double>(m[(pt.x + 1) + pt.y * size.x] - tc)) > tol) {
        spanRight = false;
      }
      pt.y++;
    }
  }
}

// _bwlabel 模板函数
template <class T>
int _bwlabel(T *src, int *res, XYPoint size) {
  XYPoint pt;
  int pos = 0;
  int idx = 1;

  for (int i = 0; i < size.x * size.y; i++) {
    res[i] = (src[i] == 0.0) ? 0 : -1;
  }

  for (int ky = 0; ky < size.y; ky++) {
    for (int kx = 0; kx < size.x; kx++, pos++) {
      if (res[pos] == -1) {
        pt.x = kx;
        pt.y = ky;
        _floodFill<int>(res, size, pt, idx, 0);
        idx++;
      }
    }
  }
  return idx - 1;  // 返回连通区域数量
}

// bwlabel 函数入口

/*
 * bwlabel - Label connected components in a binary image.
 *
 * This function is based on the `bwlabel` implementation from the EBImage package,
 * originally developed by Andrzej Oleś, Gregoire Pau, Mike Smith, Oleg Sklyar, Wolfgang Huber,
 * Joseph Barry, and Philip A. Marais.
 *
 * License: LGPL (Lesser General Public License)
 *
 * URL: https://github.com/aoles/EBImage
 */

// [[Rcpp::export]]
List bwlabel(NumericMatrix x) {
  int nrow = x.nrow();
  int ncol = x.ncol();
  XYPoint size = {nrow, ncol};

  IntegerMatrix res(nrow, ncol);

  // 调用 _bwlabel，获取连通区域数量
  int num_labels = _bwlabel(REAL(x), INTEGER(res), size);

  // 返回包含标记图像和区域数量的 List
  return List::create(Named("labeled_image") = res,
                      Named("num_labels") = num_labels);
}

// [[Rcpp::export]]
IntegerVector get_border_indices(IntegerMatrix labels, LogicalMatrix borders) {
  std::unordered_set<int> unique_labels;

  // 遍历矩阵，将边界上标签收集到集合中
  for (int i = 0; i < labels.nrow(); i++) {
    for (int j = 0; j < labels.ncol(); j++) {
      if (borders(i, j)) {
        unique_labels.insert(labels(i, j));
      }
    }
  }

  // 将集合转化为 Rcpp 的 IntegerVector
  IntegerVector border_indices(unique_labels.begin(), unique_labels.end());
  return border_indices;
}


// [[Rcpp::export]]
LogicalVector create_label_mask(IntegerVector indices, IntegerVector borders_indices) {
  std::unordered_set<int> borders_set(borders_indices.begin(), borders_indices.end());
  LogicalVector label_mask(indices.size());

  // 遍历 indices，检查每个元素是否在 borders_set 中
  for (int i = 0; i < indices.size(); i++) {
    label_mask[i] = borders_set.count(indices[i]) > 0;
  }

  return label_mask;
}



// [[Rcpp::export]]
LogicalMatrix create_clear_mask(IntegerMatrix labels, LogicalVector label_mask) {
  int nrow = labels.nrow();
  int ncol = labels.ncol();

  // 创建 mask 矩阵
  LogicalMatrix mask(nrow, ncol);

  // 逐元素遍历 labels 矩阵，根据 label_mask 生成 mask
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int label = labels(i, j);
      mask(i, j) = label_mask[label];  // 如果 label_mask[label] 为 true，则该位置需要清除
    }
  }

  return mask;
}


// [[Rcpp::export]]
NumericMatrix clear_border_pixels(NumericMatrix out, LogicalMatrix mask, double bgval) {
  int nrow = out.nrow();
  int ncol = out.ncol();

  // 遍历 out 和 mask，设置需要清除的像素为背景值
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (mask(i, j)) {
        out(i, j) = bgval;  // 设置为背景值
      }
    }
  }

  return out;
}



// [[Rcpp::export]]
NumericMatrix clear_border(NumericMatrix labels, int buffer_size = 0, double bgval = 0) {
  // 克隆输入矩阵，创建一个副本
  NumericMatrix out = clone(labels);

  // 设置扩展范围
  int ext = buffer_size + 1;

  // 获取矩阵的行数和列数
  int nrow = out.nrow();
  int ncol = out.ncol();

  // 初始化边界矩阵为 false
  LogicalMatrix borders(nrow, ncol);

  // 设置行边界（顶部和底部）
  for (int i = 0; i < ext; ++i) {
    for (int j = 0; j < ncol; ++j) {
      borders(i, j) = true;                // 顶部边界
      borders(nrow - i - 1, j) = true;     // 底部边界
    }
  }

  // 设置列边界（左侧和右侧）
  for (int i = 0; i < ext; ++i) {
    for (int j = 0; j < nrow; ++j) {
      borders(j, i) = true;                // 左侧边界
      borders(j, ncol - i - 1) = true;     // 右侧边界
    }
  }

  // 调用 bwlabel 以标记清理后的图像
  List labeled_result = bwlabel(out);
  IntegerMatrix labeled_image = labeled_result["labeled_image"];
  int num_labels = labeled_result["num_labels"];

  // 获取与边界相连的标签
  IntegerVector border_indices = get_border_indices(labeled_image, borders);
  IntegerVector indices = seq(0, num_labels); // 创建 indices 数组
  LogicalVector label_mask = create_label_mask(indices, border_indices);
  LogicalMatrix mask = create_clear_mask(labeled_image, label_mask);

  out = clear_border_pixels(out, mask, bgval);

  return out;
}

