// Matrix multiply device code
#include <assert.h>
#include <math.h>
#include "utils.h"
#include "types.h"
using namespace std;

// define size of blocked MM, C[TA,TB]
#define TA 64
#define TB 64
#define TW 16
#define NIOUT ((TA)/(BLOCKDIM_Y))
#define NJOUT ((TB)/(BLOCKDIM_X))
#define NOUT (NIOUT*NJOUT)
// nOUT = TA*TB / (nT in a tblock)
#define min(a,b) (((a)<(b))?(a):(b))

__global__ void matMul(int N, _DOUBLE_ *C, _DOUBLE_ *A, _DOUBLE_ *B) {

    int ty = threadIdx.y, tx = threadIdx.x;
    int by = blockIdx.y, bx = blockIdx.x;
    int I = by*TA + ty; int J= bx*TB + tx;

    int total = N/TW;
    if (N%TW!=0) total +=1;

if(by*TA>=N||bx*TB>=N)return; // no need

// currently do not handle corner case

    _DOUBLE_ Cij[NOUT] = {0};

    // each matrix block of C located in I, J (left upper corner)
    for (int kk=0; kk<total; kk++){
        
        __shared__ _DOUBLE_  a[TA][TW], b[TW][TB]; 
        
        // one thread block handle a matrix block，multiple ops per thread

        // load A 

        //# pragma unroll
        for(int istep=0;istep<NIOUT;istep++){

            //# pragma unroll
            for(int jstep=0;jstep<(TW/BLOCKDIM_X);jstep++){
                //assert(ty+istep<TA);
                //assert(tx+jstep<TW);    
                //assert(I+istep<N);
                //assert(kk*TW+tx+jstep<N); 
                //if(kk*TW+tx+jstep>=N) break;
                a[ty+istep*BLOCKDIM_Y][tx+jstep*BLOCKDIM_X] = A[(I+istep*BLOCKDIM_Y)*N+(kk*TW+tx+jstep*BLOCKDIM_X)];
            }
        }

        // load B 
        //# pragma unroll
        for(int istep=0;istep<TW/BLOCKDIM_Y;istep++){


            //# pragma unroll
            for(int jstep=0;jstep<NJOUT;jstep++){
                //assert(ty+istep<TW);
                //assert(tx+jstep<TB);  
                //assert(kk*TW+ty+istep<N);
                //assert(J+jstep<N);
                //if(J+jstep>=N) break;
                b[ty+istep*BLOCKDIM_Y][tx+jstep*BLOCKDIM_X] = B[(kk*TW+ty+istep*BLOCKDIM_Y)*N+(J+jstep*BLOCKDIM_X)];
            }
        }

        __syncthreads();

    // calculate multiple values of CIJ block
        //# pragma unroll
        for (int k=0; k<TW; k++){
            //# pragma unroll
            for(int i=0;i<NIOUT;i++){
                //# pragma unroll
                for(int j =0;j<NJOUT;j++){
                        //assert(ty+i*BLOCKDIM_Y<TA);
                        //assert(tx+j*BLOCKDIM_X<TB);  
                        Cij[i*NJOUT+j]+=a[ty+i*BLOCKDIM_Y][k] * b[k][tx+j*BLOCKDIM_X];
                }
            }
        }
        __syncthreads();
    }
        //# pragma unroll
        for(int i=0;i<NIOUT;i++){
           // # pragma unroll
            for(int j =0;j<NJOUT;j++){
            if (I+i*BLOCKDIM_Y<N && J+j*BLOCKDIM_X<N)
                C[I*N + J + N*(i*BLOCKDIM_Y) + j*BLOCKDIM_X ]=Cij[i*NJOUT+j];
            }
        }

}

