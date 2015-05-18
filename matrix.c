#include<stdlib.h>
#include<matrix.h>
#include<stdio.h>
float **matrix(long nrl, long nrh, long ncl, long nch){
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */

  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m){
    printf("allocation failure 1 in matrix()");
    exit(1);
  }
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl]){
    printf("allocation failure 2 in matrix()");
    exit(1);
  }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void matrix_free(float** m, long nrl, long nrh, long ncl, long nch) {
    free(m[nrl]-NR_END);
    free(m-NR_END);
}
