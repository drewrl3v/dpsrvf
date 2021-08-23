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
    float** Energy = nullptr;

/*
    constexpr int NBR_SIZ = 63;
    constexpr int Nbrs[NBR_SIZ][2] = { // the dynamic time warping neighborhood 
            {  1,  1 }, {  1,  2 }, {  1,  3 }, {  1,  4 }, {  1,  5 }, {  1,  6 }, {  1,  7 }, {  1,  8 }, {  1,  9 }, {  1, 10 },
            {  2,  1 }, {  2,  3 }, {  2,  5 }, {  2,  7 }, {  2,  9 }, {  3,  1 }, {  3,  2 }, {  3,  4 }, {  3,  5 }, {  3,  7 },
            {  3,  8 }, {  3, 10 }, {  4,  1 }, {  4,  3 }, {  4,  5 }, {  4,  7 }, {  4,  9 }, {  5,  1 }, {  5,  2 }, {  5,  3 },
            {  5,  4 }, {  5,  6 }, {  5,  7 }, {  5,  8 }, {  5,  9 }, {  6,  1 }, {  6,  5 }, {  6,  7 }, {  7,  1 }, {  7,  2 },
            {  7,  3 }, {  7,  4 }, {  7,  5 }, {  7,  6 }, {  7,  8 }, {  7,  9 }, {  7, 10 }, {  8,  1 }, {  8,  3 }, {  8,  5 },
            {  8,  7 }, {  8,  9 }, {  9,  1 }, {  9,  2 }, {  9,  4 }, {  9,  5 }, {  9,  7 }, {  9,  8 }, {  9, 10 }, { 10,  1 },
            { 10,  3 }, { 10,  7 }, { 10,  9 }
    };
 
*/


    //constexpr int NBR_SIZ = 7;
    //constexpr int Nbrs[NBR_SIZ][2] = {{1, 1}, {1, 2}, {2, 1}, {2, 3}, {3, 2},{ 1, 3},{ 3, 1}};

    //constexpr int NBR_SIZ = 7;
    //constexpr int Nbrs[NBR_SIZ][2] = {{3,2},{3,1},{2,3},{1,3},{2,1},{1,2},{1,1}};

    //constexpr int NBR_SIZ = 13;
    //constexpr int Nbrs[NBR_SIZ][2] = {{3,4},{2,4},{1,4},{4,3},{4,2},{4,1},{3,2},{3,1},{2,3},{1,3},{2,1},{1,2},{1,1}};

    //const int NBR_SIZ = 23;
    //int Nbrs[NBR_SIZ][2] = {{1, 1}, {1, 2}, {2, 1}, {2, 3}, {3, 2},{ 1, 3},{ 3, 1}, {1,4}, {3,4}, {4,3}, {4,1}, {1,5}, {2,5}, {3,5}, {4,5}, {5,4}, {5,3},{5,2}, {5,1}, {1,6}, {5,6}, {6,5}, {6,1}};

    constexpr int NBR_SIZ = 23;
    constexpr int Nbrs[NBR_SIZ][2] = {{6,1}, {6,5}, {5,6}, {1,6}, {5,1}, {5,2}, {5,3}, {5,4}, {4,5}, {3,5}, {2,5}, {1,5}, {4,1}, {4,3}, {3,4}, {1,4}, {3,1}, {1,3}, {3,2}, {2,3}, {2,1}, {1,2}, {1,1}};

    // drastic  
//         int Nbrs[][2] = 
//         {{1, 1},{1, 2},{2, 1},{2, 3},{3, 2},{1, 3},{3, 1},{1, 4},{3, 4},{4, 3},{4, 1},{1, 5},{2, 5},{3, 5},{4, 5}, 
//          {5, 4},{5, 3},{5, 2},{5, 1},{1, 6},{5, 6},{6, 5},{6, 1},{1, 7},{2, 7},{3, 7},{4, 7},{5, 7},{6, 7},{1, 6},
//          {1, 7},{1, 8},{1, 9},{1, 10},{7, 1},{8, 1},{9, 1},{10, 1},{1, 20}, {2, 20},{3, 30},
//         };
//       const int NBR_SIZ = 41;

    int i = 0;
    int j = 0;
    int Num = 0;
    float CandE[NBR_SIZ];
    int k = 0,l = 0;
    float minCandE = 10000;
    int minCandE_idx = 0;
    float** Path_x = nullptr;
    float** Path_y = nullptr;
    Path_x = (float **)malloc(T*sizeof(float *));
    //Path_x = new float** (T * sizeof(float*)) // float pointer set to value T * siazeof(float*)
    Path_y = (float **)malloc(T*sizeof(float *));

    float *x = nullptr;
    float *y = nullptr;
    float *xnew = nullptr;
    float *ynew = nullptr;
    x = (float *)malloc(T*sizeof(float));
    y = (float *)malloc(T*sizeof(float));
    float *xx1 = (float *) malloc(T*sizeof(float) );

    xnew = (float *)malloc(T*sizeof(float));
    ynew = (float *)malloc(T*sizeof(float));

    int cnt = 0;

    Energy = (float **)malloc(T*sizeof(float *));
    for(i = 0; i < T; i++)
    {
        Energy[i] = (float *)calloc(T,sizeof(float));
        Path_x[i] = (float *)calloc(T,sizeof(float));
        Path_y[i] = (float *)calloc(T,sizeof(float));

        //Forming energies associated with different paths
        Energy[0][i] = 5000000;
        Energy[i][0] = 5000000;
    }

    Energy[0][0] = 0;
    xx1[0] = 0;
    for(i = 1 ; i < T; i ++)
    {
        fflush(stdout);
        for(j = 1; j < T ; j ++)
        {
            minCandE = 10000;
            for(Num = 0; Num < NBR_SIZ; Num++)
            {
                k = i - Nbrs[Num][0];
                l = j - Nbrs[Num][1];
                if(k >= 0 && l >= 0)
                {
                    CandE[Num] = Energy[k][l] + DPcost(q1, q2,n, T, k,l,i,j);
                }
                else
                {
                    CandE[Num] = 5000000;//10000;
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

    x[0] = T-1; y[0] = T-1;
    while ( x[cnt] > 0 )
    {
        i = (int) y[cnt];
        j = (int) x[cnt];
        y[cnt + 1] = Path_x[i] [j]  ;
        x[cnt + 1] = Path_y[i] [j] ;
        cnt ++;
    }

    for (i = 0; i < cnt; i ++)
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

    //cost = 0;
    for(i = 0;i < T;i ++)
    {
        free(Energy[i]);
        free(Path_x[i]);
        free(Path_y[i]);
    }

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
    int x;
    float y;
    int y1,y2;
    float slope = 0;
    float E,E1,E2;
    float f;
    float vec11, vec12, vec21, vec22,vec23;

    E1 = 0; E2 = 0;
    E = 0;
    slope =(float ) ( i - k)/(j - l);
    if (slope == 0)
        printf("\nslope zero\n");
    float *vecarray = nullptr;
    vecarray = (float *)malloc(n * sizeof(float));

    for(x = l; x <= j ; x ++)
    {
        y = k + (x - l) * slope;
        y1 = (int )floorf(y);
        y2 = (int )ceilf(y);
        f = y - y1;
        for (int kk = 0; kk < n; kk++)
        {
            vecarray[kk] = (f*q2[kk*T + y2] + (1 - f)*q2[kk*T + y1])*sqrt(slope);
            E2 = E2 + (q1[kk*T + x] - vecarray[kk])*(q1[kk*T + x] - vecarray[kk]);
        }
    }
    E = E2/T;
    free(vecarray);
    vecarray = nullptr;
    return E;
}

void dpmatch::linint(float *xnew, float *ynew, int cnt, float *xx, float *yy, int n)
{
	int i = 0;
	int idx = 0;
	//Assume xnew and xx are sorted. 
	//Find the interval where xx[0] is located
	float m = 0;
	for (i = 0; i < n; i ++)
	{
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
				m = (ynew[idx+1]  - ynew[idx])/(xnew[idx+1]  - xnew[idx]);
				yy[i] = ynew[idx] + m * ( xx[i] - xnew[idx]);
				break;
			}
			else if(xx[i] >= xnew[cnt-1] )
			{
				//This is beyond the boundary of xnew
				// Need to extrapolate
				// Use the slope 	of points corresponding to cnt -1 and cnt -2
				m = (ynew[cnt-1]  - ynew[cnt-2])/(xnew[cnt-1]  - xnew[cnt-2]);
				yy[i] = ynew[cnt - 1] + m * ( xx[i] - xnew[cnt-1]);			
				break;
			}
			else 
				idx ++;
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
/*
float dpmatch::CostFn2(float *q1L, float *q2L,int n, int scl, int k, int l, int i, int j) {
    float m = (i-k)/(float)(j-l), sqrtm = sqrt(m), E = 0, y, tmp, ip, fp;
    int x, idx, d, iL=i*scl, kL=k*scl, lL=l*scl;

    for (x = kL; x <= iL; ++x) {
        y = (x-kL)*m + lL;
        fp = modf(y, &ip);
        idx = (int)(ip + (fp >= 0.5));

        for (d = 0; d < n; ++d) {
            tmp = q1L[n*x + d] - sqrtm*q2L[n*idx + d];
            E += tmp*tmp;
        }
    }

    return E;
}
*/
