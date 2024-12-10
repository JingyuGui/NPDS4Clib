#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double trapz_rcpp(NumericVector x, NumericVector y) {
  int m = x.size();

  // 检查 x 和 y 的长度是否一致
  if (y.size() != m) {
    stop("Arguments 'x' and 'y' must have the same length.");
  }

  // 如果 x 或 y 的长度小于等于 0，返回 0.0
  if (m <= 0) return 0.0;

  // 构造 xp 和 yp 向量
  NumericVector xp(m * 2); // xp 长度为 2*m
  NumericVector yp(m * 2); // yp 长度为 2*m


  // 填充 xp 和 yp
  for (int i = 0; i < m; ++i) {
    xp[i] = x[i];
    xp[m + i] = x[m - i - 1]; // 逆序填充
    yp[i] = 0.0;
    yp[m + i] = y[m - i - 1]; // 逆序填充 y
  }



  // 计算 p1 和 p2
  double p1 = 0.0, p2 = 0.0;

  int n = 2*m;

  for (int i = 0;i < n-1; ++i) {
    p1 += xp[i] * yp[i + 1];
    p2 += xp[i + 1] * yp[i];
  }

  p1 += xp[n-1] * yp[0];
  p2 += xp[0] * yp[n-1];

  // 返回最终结果
  return 0.5 * (p1 - p2);
}
