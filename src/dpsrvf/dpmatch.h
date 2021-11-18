#ifndef DYNAMIC_PROGRAMMING_LIB
#define DYNAMIC_PROGRAMMING_LIB

#include<vector>

class dpmatch {

public:

    dpmatch();

    ~dpmatch();

    float* match(int n, int T, float* q1, float* q2);

    float DPcost(float *q1, float *q2, int n, int T, int k, int l, int i, int j);

    void linint(const std::vector<float>& xnew, const std::vector<float>& ynew, int cnt, float* xx, float* yy, int n);
};
#endif
