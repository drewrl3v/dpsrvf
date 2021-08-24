#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "dpmatch.h"
using namespace std;

dpmatch::dpmatch()
{

}

dpmatch::~dpmatch()
{

}


float* dpmatch::match(int n, int T, float *q1, float *q2)
{
    // n is the dimension, i.e. R^n
    // T is the size (num of points along the shape)

    constexpr int NBR_SIZ = 63;
    constexpr int Nbrs[NBR_SIZ][2] = {{10,9}, {10,7}, {10,3}, {10,1}, {9,10},
                                      {9,8}, {9,7}, {9,5}, {9,4}, {9,2},
                                      {9,1}, {8,9}, {8,7}, {8,5}, {8,3},
                                      {8,1}, {7,10}, {7,9}, {7,8}, {7,6}, 
                                      {7,5}, {7,4}, {7,3}, {7,2}, {7,1}, 
                                      {6,7}, {6,5}, {6,1}, {5,9}, {5,8}, 
                                      {5,7}, {5,6}, {5,4}, {5,3}, {5,2}, 
                                      {5,1}, {4,9}, {4,7}, {4,5}, {4,3}, 
                                      {4,1}, {3,10}, {3,8}, {3,7}, {3,5}, 
                                      {3,4}, {3,2}, {3,1}, {2,9}, {2,7}, 
                                      {2,5}, {2,3}, {2,1}, {1,10}, {1,9}, 
                                      {1,8}, {1,7}, {1,6}, {1,5}, {1,4}, {1,3}, {1,2}, {1,1}};

    //constexpr int NBR_SIZ = 7;
    //constexpr int Nbrs[NBR_SIZ][2] = {{1, 1}, {1, 2}, {2, 1}, {2, 3}, {3, 2},{ 1, 3},{ 3, 1}};

    //constexpr int NBR_SIZ = 7;
    //constexpr int Nbrs[NBR_SIZ][2] = {{3,2},{3,1},{2,3},{1,3},{2,1},{1,2},{1,1}};

    //const int NBR_SIZ = 23;
    //int Nbrs[NBR_SIZ][2] = {{1, 1}, {1, 2}, {2, 1}, {2, 3}, {3, 2},{ 1, 3},{ 3, 1}, {1,4}, {3,4}, {4,3}, {4,1}, {1,5}, {2,5}, {3,5}, {4,5}, {5,4}, {5,3},{5,2}, {5,1}, {1,6}, {5,6}, {6,5}, {6,1}};

    // Looks like this works
    //constexpr int NBR_SIZ = 23;
    //constexpr int Nbrs[NBR_SIZ][2] = {{6,1}, {6,5}, {5,6}, {1,6}, {5,1}, {5,2}, {5,3}, {5,4}, {4,5}, {3,5}, {2,5}, {1,5}, {4,1}, {4,3}, {3,4}, {1,4}, {3,1}, {1,3}, {3,2}, {2,3}, {2,1}, {1,2}, {1,1}};


    //Forming energies associated with different paths
    float** Path_x = {(float **)malloc(T*sizeof(float *))};
    float** Path_y = {(float **)malloc(T*sizeof(float *))};
    float** Energy = {(float **)malloc(T*sizeof(float *))};
    for(int i = 0; i < T; i++)
    {
        Energy[i] = (float *)calloc(T,sizeof(float));
        Path_x[i] = (float *)calloc(T,sizeof(float));
        Path_y[i] = (float *)calloc(T,sizeof(float));
        Energy[0][i] = 5000000;
        Energy[i][0] = 5000000;
    }
    Energy[0][0] = 0; // This is the starting point of the dynamic time warping path

    float CandE[NBR_SIZ]{0};
    int minCandE_idx = 0;
    float* xx1 = (float *) malloc(T*sizeof(float));
    xx1[0] = 0;
    for(int i = 1 ; i < T; i ++)
    {
        for(int j = 1; j < T ; j ++)
        {
            float minCandE{10000};
            for(int Num = 0; Num < NBR_SIZ; Num++)
            {
                int k = {i - Nbrs[Num][0]};
                int l = {j - Nbrs[Num][1]};
                if(k >= 0 && l >= 0)
                {
                    CandE[Num] = Energy[k][l] + DPcost(q1, q2,n, T, k,l,i,j);
                }
                else
                {
                    CandE[Num] = 10000;//5000000;//10000;
                }
                if(CandE[Num] < minCandE )
                {
                    minCandE = CandE[Num];
                    minCandE_idx = Num;
                }
                Energy[i][j] = minCandE;

                Path_x[i][j] = i - Nbrs[minCandE_idx][0];
                Path_y[i][j] = j - Nbrs[minCandE_idx][1];
            }
        }
        xx1[i] = (float ) i/(T - 1);
    }


    float* x = {(float *)malloc(T*sizeof(float))};
    float* y = {(float *)malloc(T*sizeof(float))};
    x[0] = T-1;
    y[0] = T-1;
    int cnt{0};
    while ( x[cnt] > 0 )
    {
        int i = {(int) y[cnt]};
        int j = {(int) x[cnt]};
        y[cnt + 1] = Path_x[i][j];
        x[cnt + 1] = Path_y[i][j];
        cnt++;
    }

    float* xnew = {(float *)malloc(T*sizeof(float))};
    float* ynew = {(float *)malloc(T*sizeof(float))};
    for (int i = 0; i < cnt; i ++)
    {
        xnew[i] = (x[cnt-i-1] -x[cnt-1] )/(x[0] - x[cnt - 1]);
        ynew[i] = (y[cnt-i-1] -y[cnt-1] )/(y[0] - y[cnt - 1]);
    }

    xnew[cnt-1] = 1;
    ynew[cnt-1] = 1;
    xnew[0] = 0;
    ynew[0] = 0;
    float* gamma = new float [T];
    linint(xnew, ynew, cnt, xx1, gamma, T);

    free(Path_x);
    free(Path_y);
    free(Energy);

    free(x);
    free(y);
    free(xx1);
    free(xnew);
    free(ynew);

    return(gamma);
}


float dpmatch::DPcost(float *q1, float *q2, int n, int T, int k, int l, int i, int j)
{
    //float slope{(float) (i - k)/(j - l)};
    float E2{0};
    float* vecarray = {(float* )malloc(n * sizeof(float))};

    for(int x = l; x <= j; x++)
    {
        float slope{(float) (i - k)/(j - l)};
        float y = {k + (x - l) * slope};
        int y1 = {(int)floorf(y)};
        int y2 = {(int)ceilf(y)};
        float f = {y - y1};

        for (int kk = 0; kk < n; kk++)
        {
            vecarray[kk] = (f*q2[kk*T + y2] + (1 - f)*q2[kk*T + y1])*sqrt(slope);
            E2 = E2 + (q1[kk*T + x] - vecarray[kk])*(q1[kk*T + x] - vecarray[kk]);
        }
    }
    float E{E2/T};
    free(vecarray);
    vecarray = NULL;
    return E;
}

void dpmatch::linint(float *xnew, float *ynew, int cnt, float *xx, float *yy, int n)
{
	//Assume xnew and xx are sorted. 
	//Find the interval where xx[0] is located
	float slope{0};
  int idx{0};
	for (int i = 0; i < n; i++)
	{
    //int idx{0}; // Not sure if this shoud be inside the for-loop
		while(idx < cnt - 1)
		{
			if(xx[i] >= xnew[idx] && xx[i] <= xnew[idx+1] )
			{
				yy[i] = ynew[idx] + (xx[i] - xnew[idx]) *( ynew[idx+1] - ynew[idx] ) / ( xnew[idx+1] - xnew[idx] );
				break;
			}
			else if ( xx[i] <= xnew[idx+1] )
			{
				// Need to extrapolate
				// Use the slope 	of points corresponding to idx and idx + 1
				slope = (ynew[idx+1]  - ynew[idx])/(xnew[idx+1]  - xnew[idx]);
				yy[i] = ynew[idx] + slope * ( xx[i] - xnew[idx]);
				break;
			}
			else if(xx[i] >= xnew[cnt-1] )
			{
				//This is beyond the boundary of xnew
				// Need to extrapolate
				// Use the slope 	of points corresponding to cnt -1 and cnt -2
				slope = (ynew[cnt-1]  - ynew[cnt-2])/(xnew[cnt-1]  - xnew[cnt-2]);
				yy[i] = ynew[cnt - 1] + slope * ( xx[i] - xnew[cnt-1]);			
				break;
			}
			else 
				idx++;
			}
		}
}

shape::shape( int v_iT)
{
    m_pfPhi = (float *)malloc(v_iT*sizeof(float));
    m_pfTheta = (float *)malloc(v_iT*sizeof(float));
    m_v11 = (float *)malloc(v_iT*sizeof(float));
    m_v12 = (float *)malloc(v_iT*sizeof(float));
    m_v13 = (float *)malloc(v_iT*sizeof(float));
    m_iT = v_iT;

}

shape::~shape()
{
    free(m_pfPhi);
    free(m_pfTheta);
    for (int i=0; i<m_n; i++)
        free(arr[i]);
    free(arr);

}
