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

    _DOUBLE_ Cij = 0;
    for (int kk=0; kk<N/TW; kk++){
        a[ty][tx] = A[I * N + kk*TW+tx];
        b[ty][tx] = A[N*(kk*TW+ty)+J];
        __syncthreads();
        for (int k=0; k<TW; k++)
            Cij+= a[ty][k] * b[k][tx];
        __syncthreads();
    }
    C[I*N + J] = Cij;
}

