#include "mpi.h"
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include "runtime_parameters.h"

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
