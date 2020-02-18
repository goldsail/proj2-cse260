#define TA 64
#define TB 64
#define TW 16

void setGrid(int n, dim3 &blockDim, dim3 &gridDim)
{
   // set your block dimensions and grid dimensions here
   gridDim.x = n / TB;
   gridDim.y = n / TA;
   if(n % TB != 0)
   	gridDim.x++;
   if(n % TA != 0)
    	gridDim.y++;
}
