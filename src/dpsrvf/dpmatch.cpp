#include "dpmatch.h"
#include <math.h>
#include <stdlib.h>
#include <vector>

dpmatch::dpmatch() = default;

dpmatch::~dpmatch() = default;

float *dpmatch::match(int n, int T, float *q1, float *q2) {
  // n is the dimension, i.e. R^n.
  // T is the size (num of points along the shape).

  // Initialization with experimentally determined weights.
  constexpr int NBR_SIZ = 63;
  constexpr int Nbrs[NBR_SIZ][2] = {
      {1, 1}, {1, 2},  {1, 3},  {1, 4},  {1, 5},  {1, 6},  {1, 7},  {1, 8},
      {1, 9}, {1, 10}, {2, 1},  {2, 3},  {2, 5},  {2, 7},  {2, 9},  {3, 1},
      {3, 2}, {3, 4},  {3, 5},  {3, 7},  {3, 8},  {3, 10}, {4, 1},  {4, 3},
      {4, 5}, {4, 7},  {4, 9},  {5, 1},  {5, 2},  {5, 3},  {5, 4},  {5, 6},
      {5, 7}, {5, 8},  {5, 9},  {6, 1},  {6, 5},  {6, 7},  {7, 1},  {7, 2},
      {7, 3}, {7, 4},  {7, 5},  {7, 6},  {7, 8},  {7, 9},  {7, 10}, {8, 1},
      {8, 3}, {8, 5},  {8, 7},  {8, 9},  {9, 1},  {9, 2},  {9, 4},  {9, 5},
      {9, 7}, {9, 8},  {9, 10}, {10, 1}, {10, 3}, {10, 7}, {10, 9}};


  // Initialize different paths.
  std::vector<std::vector<float>> Path_x(T, std::vector<float>(T, 0));
  std::vector<std::vector<float>> Path_y(T, std::vector<float>(T, 0));

  // Initialize Energies in cache friendly way.
  std::vector<std::vector<float>> Energy(T, std::vector<float>(T));
  for (int i = 0; i < T; i++) {
    Energy[i][0] = 50000000000;
  }
  for (int i = 0; i< T; i++){
    Energy[0][i] = 50000000000;
  } Energy[0][0] = 0;

  float CandE[NBR_SIZ]{0};
  for (int i = 1; i < T; i++) {
    for (int j = 1; j < T; j++) {
      float minCandE{100000};
      int minCandE_idx = 0;
      for (int Num = 0; Num < NBR_SIZ; Num++) {
        int k = i - Nbrs[Num][0];
        int l = j - Nbrs[Num][1];
        if (k >= 0 && l >= 0) {
          CandE[Num] = Energy[k][l] + DPcost(q1, q2, n, T, k, l, i, j);
          if (Num == 0 || CandE[Num] < minCandE) {
            minCandE = CandE[Num];
            minCandE_idx = Num;
          }
        }
        Energy[i][j] = minCandE;
        Path_x[i][j] = i - Nbrs[minCandE_idx][0];
        Path_y[i][j] = j - Nbrs[minCandE_idx][1];
      }
    }
  }

  std::vector<float> x(T);
  std::vector<float> y(T);
  x[0] = T - 1;
  y[0] = T - 1;
  int cnt = 0;
  while (x[cnt] > 0) {
    int i = (int)y[cnt];
    int j = (int)x[cnt];
    y[cnt + 1] = Path_x[i][j];
    x[cnt + 1] = Path_y[i][j];
    cnt++;
  }

  std::vector<float> xnew(T);
  std::vector<float> ynew(T);
  for (int i = 0; i < cnt; i++) {
    xnew[i] = (x[cnt - i - 1] - x[cnt - 1]) / (x[0] - x[cnt - 1]);
    ynew[i] = (y[cnt - i - 1] - y[cnt - 1]) / (y[0] - y[cnt - 1]);
  }

  xnew[cnt - 1] = 1;
  ynew[cnt - 1] = 1;
  xnew[0] = 0;
  ynew[0] = 0;

  float *gamma = new float[T]; // Store values for gamma.
  float *xx1 = new float[T];
  for (int i = 0; i < T; i++) {
    xx1[i] = (float)i / (T - 1);
  }
  // Perform linear interpolation.
  linint(xnew, ynew, cnt, xx1, gamma, T);
  delete[] xx1;

  return (gamma);
}

float dpmatch::DPcost(float *q1, float *q2, int n, int T, int k, int l, int i,
                      int j) {
  float E2 = 0;
  std::vector<float> vecarray(n);

  for (int x = l; x <= j; ++x) {
    float slope = (float)(i - k) / (j - l);
    float y = k + (x - l) * slope;
    int y1 = (int)floorf(y);
    int y2 = (int)ceilf(y);
    float f = y - y1;

    for (int kk = 0; kk < n; ++kk) {
      vecarray[kk] = (f * q2[kk * T + y2] +
                      (1 - f) * q2[kk * T + y1] * sqrt(slope));
      E2 = E2 +
           (q1[kk * T + x] - vecarray[kk]) * (q1[kk * T + x] - vecarray[kk]);
    }
  }
  float E = E2 / T;
  return E;
}
/*
Precondition:
x_new: sorted, xx: sorted
*/
void dpmatch::linint(const std::vector<float> &xnew,
                     const std::vector<float> &ynew, int cnt, float *xx,
                     float *yy, int n) {
  // Find the interval where xx[0] is located.
  float slope = 0;
  for (int i = 0; i < n; i++) {
    // TODO(D-moe) Determine if initialization shoud be inside for loop.
    int idx = 0;
    while (idx < cnt - 1)
    {
      if (xx[i] >= xnew[idx] && xx[i] <= xnew[idx + 1]) {
        yy[i] = ynew[idx] + (xx[i] - xnew[idx]) * (ynew[idx + 1] - ynew[idx]) /
                                (xnew[idx + 1] - xnew[idx]);
        slope = (ynew[idx + 1] - ynew[idx])/
                (xnew[idx + 1] - xnew[idx]);
        break;
      } else if (xx[i] <= xnew[idx + 1]) {
        // Need to extrapolate.
        // Use the slope 	of points corresponding to idx and idx + 1.
        slope = (ynew[idx + 1] - ynew[idx]) / (xnew[idx + 1] - xnew[idx]);
        yy[i] = ynew[idx] + slope * (xx[i] - xnew[idx]);
        break;
      } else if (xx[i] >= xnew[cnt - 1]) {
        // This is beyond the boundary of xnew.
        // Need to extrapolate.
        // Use the slope of points corresponding to cnt -1 and cnt -2.
        slope =
            (ynew[cnt - 1] - ynew[cnt - 2]) / (xnew[cnt - 1] - xnew[cnt - 2]);
        yy[i] = ynew[cnt - 1] + slope * (xx[i] - xnew[cnt - 1]);
        break;
      }
      // Else.
      idx++;
    }
  }
}
