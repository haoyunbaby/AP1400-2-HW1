#include "hw1.h"

#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

using namespace std;

namespace algebra {
Matrix zeros(size_t n, size_t m) {
  Matrix mat(n, vector<double>(m, 0));
  return mat;
}
Matrix ones(size_t n, size_t m) {
  Matrix mat(n, vector<double>(m, 1));
  return mat;
}
Matrix random(size_t n, size_t m, double min, double max) {
  if (min > max) {
    throw logic_error("random");
  }
  Matrix mat(n, vector<double>(m));
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> dis(min, max);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      mat[i][j] = dis(gen);
    }
  }
  return mat;
}
void show(const Matrix& matrix) {
  cout << fixed << setprecision(3);
  int n = matrix.size();
  int m = matrix[0].size();
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      cout << matrix[i][j] << ' ';
    }
    cout << '\n';
  }
}
Matrix multiply(const Matrix& matrix, double c) {
  Matrix mat = matrix;
  int n = mat.size();
  int m = mat[0].size();
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      mat[i][j] *= c;
    }
  }
  return mat;
}
Matrix multiply(const Matrix& matrix1, const Matrix& matrix2) {
  if (matrix1.empty() || matrix2.empty()) {
    return Matrix{};
  }
  int h = matrix1.size(), j = matrix1[0].size(), k = matrix2.size(),
      l = matrix2[0].size();
  if (j != k) {
    throw logic_error("multiply");
  }
  Matrix mat = zeros(h, l);
  for (size_t row = 0; row < h; row++) {
    for (size_t col = 0; col < l; col++) {
      double sum = 0;
      for (size_t t = 0; t < j; t++) {
        sum += matrix1[row][t] * matrix2[t][col];
      }
      mat[row][col] = sum;
    }
  }
  return mat;
}
Matrix sum(const Matrix& matrix, double c) {
  if (matrix.empty()) {
    return Matrix{};
  }
  Matrix mat = matrix;
  int n = mat.size();
  int m = mat[0].size();
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      mat[i][j] += c;
    }
  }
  return mat;
}
Matrix sum(const Matrix& matrix1, const Matrix& matrix2) {
  if (matrix1.empty() && matrix2.empty()) {
    return Matrix{};
  }
  if (matrix1.empty() || matrix2.empty()) {
    throw logic_error("sum");
  }
  int h = matrix1.size(), j = matrix1[0].size(), k = matrix2.size(),
      l = matrix2[0].size();

  if (h != k || j != l) {
    throw logic_error("sum");
  }
  Matrix mat = zeros(h, j);
  for (size_t row = 0; row < h; row++) {
    for (size_t col = 0; col < j; col++) {
      mat[row][col] = matrix1[row][col] + matrix2[row][col];
    }
  }
  return mat;
}
Matrix transpose(const Matrix& matrix) {
  if (matrix.empty()) {
    return Matrix{};
  }
  int n = matrix.size();
  int m = matrix[0].size();
  Matrix mat = zeros(m, n);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      mat[j][i] = matrix[i][j];
    }
  }
  return mat;
}
Matrix minor(const Matrix& matrix, size_t n, size_t m) {
  int nr = matrix.size();
  int nc = matrix[0].size();
  Matrix mat = zeros(nr - 1, nc - 1);
  int cr = 0;
  for (size_t i = 0; i < nr; i++) {
    if (i == n) {
      continue;
    }
    int cc = 0;
    for (size_t j = 0; j < nc; j++) {
      if (j == m) {
        continue;
      }
      mat[cr][cc] = matrix[i][j];
      cc++;
    }
    cr++;
  }
  return mat;
}
double determinant(const Matrix& matrix) {
  if (matrix.empty()) {
    return 1;
  }
  int n = matrix.size();
  int m = matrix[0].size();
  if (n != m) {
    throw logic_error("determinant");
  }
  if (n == 1) {
    return matrix[0][0];
  }
  double ans = 0;
  for (size_t j = 0; j < m; j++) {
    if (j % 2) {
      ans -= matrix[0][j] * determinant(minor(matrix, 0, j));
    } else {
      ans += matrix[0][j] * determinant(minor(matrix, 0, j));
    }
  }
  return ans;
}
Matrix inverse(const Matrix& matrix) {
  if (matrix.empty()) {
    return Matrix{};
  }
  double det = determinant(matrix);
  if (!det) {
    throw logic_error("inverse");
  }
  int n = matrix.size();
  Matrix mat = zeros(n, n);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      mat[i][j] = determinant(minor(matrix, i, j));
    }
  }
  mat = transpose(mat);
  mat = multiply(mat, 1 / det);
  return mat;
}
Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis = 0) {
  if (matrix1.empty() && matrix2.empty()) {
    return Matrix{};
  }
  if (matrix1.empty()) {
    Matrix mat = matrix2;
    return mat;
  }
  if (matrix2.empty()) {
    Matrix mat = matrix1;
    return mat;
  }
  int h = matrix1.size(), j = matrix1[0].size(), k = matrix2.size(),
      l = matrix2[0].size();
  if (axis) {
    if (h != k) {
      throw logic_error("concatenate");
    }
    int target = j + l;
    Matrix mat = zeros(h, target);
    for (size_t row = 0; row < h; row++) {
      for (size_t col = 0; col < target; col++) {
        if (col < j) {
          mat[row][col] = matrix1[row][col];
        } else {
          mat[row][col] = matrix2[row][col - j];
        }
      }
    }
    return mat;
  } else {
    if (j != l) {
      throw logic_error("concatenate");
    }
    int target = h + k;
    Matrix mat = zeros(target, j);
    for (size_t row = 0; row < target; row++) {
      if (row < h) {
        mat[row] = matrix1[row];
      } else {
        mat[row] = matrix2[row - h];
      }
    }
    return mat;
  }
}
Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2) {
  int n = matrix.size();
  if (r1 < 0 || r1 >= n) {
    throw logic_error("ero_swap");
  }
  if (r2 < 0 || r2 >= n) {
    throw logic_error("ero_swap");
  }
  Matrix mat = matrix;
  swap(mat[r1], mat[r2]);
  return mat;
}
Matrix ero_multiply(const Matrix& matrix, size_t r, double c) {
  Matrix mat = matrix;
  int m = matrix[r].size();
  for (size_t j = 0; j < m; j++) {
    mat[r][j] *= c;
  }
  return mat;
}
Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2) {
  Matrix mat = matrix;
  int m = matrix[0].size();
  for (size_t j = 0; j < m; j++) {
    mat[r2][j] += mat[r1][j] * c;
  }
  return mat;
}
Matrix upper_triangular(const Matrix& matrix) {
  Matrix mat = matrix;
  if (mat.empty()) {
    return mat;
  }
  int n = matrix.size();
  int m = matrix[0].size();
  if (n != m) {
    throw logic_error("upper_triangular");
  }
  for (size_t i = 0; i < n - 1; i++) {
    if (!mat[i][i]) {
      for (size_t j = i + 1; j < n; j++) {
        if (mat[j][i]) {
          mat = ero_swap(mat, i, j);
          break;
        }
      }
    }
    for (size_t j = i + 1; j < n; j++) {
      double factor = -mat[j][i] / mat[i][i];
      mat = ero_sum(mat, i, factor, j);
    }
  }
  return mat;
}
}  // namespace algebra
