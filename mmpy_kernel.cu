// Matrix multiply device code
#include <assert.h>
#include <math.h>
#include "utils.h"
#include "types.h"
using namespace std;

// define size of shared memory
#define TW 32

__global__ void matMul(int N, _DOUBLE_ *C, _DOUBLE_ *A, _DOUBLE_ *B) {

    __shared__ _DOUBLE_  a[TW][TW], b[TW][TW];

    int ty = threadIdx.y, tx = threadIdx.x;
    int by = blockIdx.y, bx = blockIdx.x;
    int I = by*TW + ty; int J= bx*TW + tx;

    int total = N/TW;
    if (N%TW!=0) total +=1;

    _DOUBLE_ Cij[4] = {0,0,0,0};

if(by < N/TW){ 
    for (int kk=0; kk<total; kk++){
        a[ty][tx] = A[I * N + kk*TW+tx];
        b[ty][tx] = B[N*(kk*TW+ty)+J];
        a[ty+TW/4][tx] = A[I * N + kk*TW+tx + N*(TW/4)];
        b[ty+TW/4][tx] = B[N*(kk*TW+ty)+J + N*(TW/4)];
        a[ty+TW/2][tx] = A[I * N + kk*TW+tx + N*(TW/2)];
        b[ty+TW/2][tx] = B[N*(kk*TW+ty)+J + N*(TW/2)];
        a[ty+TW*3/4][tx] = A[I * N + kk*TW+tx + N*(TW*3/4)];
        b[ty+TW*3/4][tx] = B[N*(kk*TW+ty)+J + N*(TW*3/4)];
        __syncthreads();

    # pragma unroll
        for (int k=0; k<TW; k++){
            Cij[0]+= a[ty][k] * b[k][tx];
            Cij[1]+= a[ty+TW/4][k] * b[k][tx];
            Cij[2]+= a[ty+TW/2][k] * b[k][tx];
            Cij[3]+= a[ty+TW*3/4][k] * b[k][tx];
        }
        __syncthreads();
    }
    C[I*N + J] = Cij[0];
    C[I*N + J + N*(TW/4)] = Cij[1];
    C[I*N + J + N*(TW/2)] = Cij[2];
    C[I*N + J + N*(TW*3/4)] = Cij[3];
}
}
