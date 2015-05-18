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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "comm.h"

void read_scattering_matrix(char *sm_file, float **scattering_matrix, int nb){
    int i,j, n, m;
    float val;
    FILE * file;

    file=fopen(sm_file,"r");
    fscanf(file,"%d",&n);
    fscanf(file,"%d",&m);
    assert (n == nb);
    assert (n == m);

    for(i=0;i<n;i++) {
      for(j=0;j<m;j++) {
        fscanf(file,"%f",&val);
        scattering_matrix[i][j] = val;
        printf("%f ", scattering_matrix[i][j]);
      }
      printf("\n");
    }

    fclose(file);
}

/*---------------------------------------------------
    - void init_xsdata(inout xsdata, in nl)
    intialiaze cross section data. for purposes
    of this kernel app the values are arbitrary, only
    size of data structure matters.
    ------------------------------------------------------*/
void init_xsdata(float *xsdata, int nl){
    int i;
    for (i = 0; i < nl; ++i)
        xsdata[i] = 1.0;
}

void init_particles(Particle *p, long npl, int mype){
    long i;
    for (i=0;i<npl;++i){
      p[i].energy = 2.0;
      p[i].absorbed = FALSE;
      p[i].band = 0;
      p[i].proc = mype;
    }
}

void set_energy_ranges(Range *energy_ranges, int nm){
    int i;
    float delta = 2.0/nm;
    energy_ranges[0].lower=0.0;
    energy_ranges[0].upper=delta;
    for (i=1;i<nm;++i){
      energy_ranges[i].lower = energy_ranges[i-1].upper;
      energy_ranges[i].upper = energy_ranges[i].lower + delta;
    }
}

long tot(long f[], long n){
    long sum = 0;
    int i;
    for (i=0;i<n;++i)
      sum += f[i];
    return sum;
}

double sum(double* f, int n){
    double sum = 0;
    int i;
    for (i=0;i<n;++i)
      sum += f[i];
    return sum;
}

void create_scattering_matrix(float **scattering_matrix, int nb){
    //create pdf
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j < nb; ++j) {
            if (i > j)
                scattering_matrix[i][j] = 0.0;
            else
                scattering_matrix[i][j] = 1.0/(nb-i);
        }
    }

    //convert to cdf for mini-app
    for (int i = 0; i < nb; ++i) {
        for (int j = i + 1; j < nb; ++j) {
            scattering_matrix[i][j] += scattering_matrix[i][j-1];
        }
        scattering_matrix[i][nb-1] = 1.0; /* So as to avoid roundoff error */
    }
}

