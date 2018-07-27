#include<stdlib.h>
#include<R.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// stage2 Molecular Complex Prediction
//     procedure MCODE-FIND-COMPLEX
void complex(int* neighbor,int* neighbor_indx,float* vertex_weight,float* D,int* seed_vertex,int* seen,int* COMPLEX)
{
    int seed,i;
    seed=*seed_vertex;
    if(!seen[seed]){
        seen[seed]=1;
        COMPLEX[seed]=1;
        if((neighbor_indx[seed+1]-neighbor_indx[seed])>0){
            //// neighborhood of vertex 'seed_vertex'
            for(i=neighbor_indx[seed];i<neighbor_indx[seed+1];i++){
                if(vertex_weight[neighbor[i]] >= (vertex_weight[seed] * (1-*D))){
                    complex(neighbor,neighbor_indx,vertex_weight,D,&(neighbor[i]),seen,COMPLEX);
                }
            }
        }
    }
}
