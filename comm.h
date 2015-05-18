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

#ifndef _COMM_H_
#define _COMM_H_

#define MAX_FILENAME_CHARS = 128
typedef enum boolean{
    FALSE,
    TRUE
} Boolean;

typedef struct range{
    float lower;
    float upper;
} Range;

typedef struct proc{
    int lrank;
    int grank;
    int type;
} Proc;

typedef struct Particle_type{
    int band;
    double energy;
    Boolean absorbed;
    int proc;
} Particle;

typedef struct {
    double val;
    int rank;
} DoubleIntPair;

extern float **matrix(long nrl, long nrh, long ncl, long nch);
extern void create_scattering_matrix(float **scattering_matrix, int nb);
extern void read_scattering_matrix(char* sm_file, float **scattering_matrix, int nb);
extern void init_xsdata(float *xsdata, int nl);
extern void init_particles(Particle *p, long npl, int mype);
extern void set_energy_ranges(Range*, int);
extern long tot(long f[], long n);
extern double sum(double* f, int n);

extern void process_input(int argc, char **argv,
    int* nm, int* nb,
    long* gsizeb,
    long* npg,
    float* tracking_rate, int *use_file,
    char *sm_file);

#define NR_END 1
extern float **matrix(long nrl, long nrh, long ncl, long nch);
extern void matrix_free(float** m, long nrl, long nrh, long ncl, long nch);

extern int read_params(char *paramfile);
extern void list_params();
extern int get_param( char *name, void *val);

#endif
