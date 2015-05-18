/*
 Copyright (c) 2015 UChicago Argonne, LLC

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/

#include <stdlib.h>
#include <stdio.h>
#include "comm.h"

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
