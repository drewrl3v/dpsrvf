#include <cstdio>
#include <vector>

class dpmatch {

public:

    dpmatch();

    ~dpmatch();

    float* match(int n, int T, float* q1, float* q2);

    float DPcost(float *q1, float *q2, int n, int T, int k, int l, int i, int j);

    void linint(const std::vector<float>& xnew, const std::vector<float>& ynew, int cnt, float *xx, float *yy, int n); // (float *xnew, float *ynew, int cnt, float *xx, float *yy, int n)



};


class shape
{
public:
    //Member variables
    float *m_pfPhi;
    float *m_pfTheta;
    float *m_v11;
    float *m_v12;
    float *m_v13;
    float **v;
    float *h;
    int m_iT;
    int m_n;

    //Member functions
    shape();//default constructor
    shape(int);
    ~shape();
    float **arr;
};
