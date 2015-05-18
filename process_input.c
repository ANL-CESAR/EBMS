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

#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include<mpi.h>
#include "comm.h"

void process_input(int argc, char **argv,
		   int* r, int* nb,
		   long* gsizekb,
		   long* npg,
		   float* tracking_rate, int *use_file,
		   char *sm_file){

  int nprocs;
  int ret;
  char *param_file;
  long npgk;
  long gsizemb;

  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  if (argc != 2){
    printf("must enter input parameter file\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  param_file = argv[1];

  read_params(param_file);
  //  list_params();

  /* check that all expected params are defined */
  assert(get_param("r", r) != -1);
  assert(get_param("nb", nb) != -1);

  assert(get_param("gsizemb", &gsizemb) != -1);
  *gsizekb = gsizemb * 1024;

  assert(get_param("npgk", &npgk) != -1);
  *npg = npgk * 1000;

  assert(get_param("tracking_rate", tracking_rate) != -1);
  ret = get_param("sm_file", sm_file);
  if (ret == -1)              /* do some consistency checks */
    *use_file = 0;
  else
    *use_file = 1;

  return;

}
