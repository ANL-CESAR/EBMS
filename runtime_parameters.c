#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include "runtime_parameters.h"

#define MAX_PARMS  500
#define MAX_LINE   100
#define MAX_TOKENS 100


typedef enum vartype_{
  TYPE_INT, TYPE_LONG, TYPE_FLOAT, TYPE_DOUBLE, TYPE_STRING} vartype;

struct{
  vartype  type;
  char     name[32];
  void     *val;
} param[MAX_PARMS];

static int get_tokens(char *s, char *tokens[]);
static int nparms = 0;

/*

int main(){
  int nx, ny;
  char filehead[128];
  void *dat;
  read_params("test.par");
  list_params();
  get_param("nx", &nx);
  get_param("ny", &ny);
  get_param("filehead", filehead);
  printf("nx: %d\n", nx);
  printf("ny: %d\n", ny);
  printf("filehead: %s\n", filehead);
}
*/

int read_params(char *paramfile){
  char   line  [MAX_LINE];
  char*  tokens[MAX_TOKENS];
  int    ntokens;
  FILE*  f;
  char*  last_token;
  int    len_last_token;

  f = fopen(paramfile, "r");

  assert(f);

  while (fgets(line, MAX_LINE, f) != NULL){
    ntokens = get_tokens(line, tokens);

    /* hack to "remove" trailing whitepace */

    last_token = tokens[ntokens - 1];
    len_last_token  = strlen(last_token);
    last_token[len_last_token - 1] = '\0';

    /*    assert(ntokens <= MAX_TOKENS);*/

    if  ( strcmp(tokens[0], "int") == 0 ){
      param[nparms].type   = TYPE_INT;
      strcpy(param[nparms].name,tokens[1]);
      param[nparms].val    = malloc(sizeof(int));
      * ((int *) param[nparms].val) = atoi(tokens[2]);
      ++nparms;
    }
    else if ( strcmp(tokens[0], "long") == 0 ){
      param[nparms].type   = TYPE_LONG;
      strcpy(param[nparms].name,tokens[1]);
      param[nparms].val    = malloc(sizeof(long));
      * ((long *) param[nparms].val) = atol(tokens[2]);
      ++nparms;
    }
    else if ( strcmp(tokens[0], "float") == 0 ){
      param[nparms].type   = TYPE_FLOAT;
      strcpy(param[nparms].name,tokens[1]);
      param[nparms].val    = malloc(sizeof(float));
      * ((float *) param[nparms].val) = atof(tokens[2]);
      ++nparms;
    }

    else if ( strcmp(tokens[0], "double") == 0 ){
      param[nparms].type   = TYPE_DOUBLE;
      strcpy(param[nparms].name,tokens[1]);
      param[nparms].val    = malloc(sizeof(double));
      * ((double *) param[nparms].val) = strtod(tokens[2],NULL);
      //      * ((double *) param[nparms].val) = atof(tokens[2]);
      ++nparms;
    }

    else if ( strcmp(tokens[0], "string") == 0 ){
      param[nparms].type = TYPE_STRING;
      strcpy(param[nparms].name,tokens[1]);
      param[nparms].val    = malloc(strlen(tokens[2]));
      strcpy( param[nparms].val, tokens[2] );
      ++nparms;
    }
    else if (tokens[0][0] == '#') {
      //printf("ignoring comment line: %s\n",line);
    } else {
      //printf("ignoring invalid line: %s\n", line);
    }
  }

  fclose(f);
  return(nparms);
}


int get_param( char *name, void* val){
  int i;
  for (i = 0; i < nparms; ++i){

    if (strcmp(param[i].name, name) == 0){
      if (param[i].type == TYPE_INT){
	* ( (int *) val ) = * ( (int *) param[i].val);
	return(0);
      }
      if (param[i].type == TYPE_LONG){
	* ( (long *) val ) = * ( (long *) param[i].val);
	return(0);
      }
      else if (param[i].type == TYPE_FLOAT){
	* ( (float *) val ) = * ( (float *) param[i].val);
	return(0);
      }

      else if (param[i].type == TYPE_DOUBLE){
	* ( (double *) val ) = * ( (double *) param[i].val);
	return(0);
      }

      else if (param[i].type == TYPE_STRING){
	strcpy( (char *) val, (char *) param[i].val );
	return(0);
      }
    }

  }
  return(-1);
}


void list_params(){
  int i;
  char *fmt;

  for (i = 0; i < nparms; ++i){
    if      (param[i].type == TYPE_INT){
      fmt = "%s %d\n";
      printf(fmt, param[i].name, * ( (int *) param[i].val));
    }
    if      (param[i].type == TYPE_LONG){
      fmt = "%s %ld\n";
      printf(fmt, param[i].name, * ( (long *) param[i].val));
    }
    else if (param[i].type == TYPE_FLOAT){
      fmt = "%s %f\n";
      printf(fmt, param[i].name, * ( (float *) param[i].val));
    }
    else if (param[i].type == TYPE_DOUBLE){
      fmt = "%s %f\n";
      printf(fmt, param[i].name, * ( (double *) param[i].val));
    }
    else if (param[i].type == TYPE_STRING){
      fmt = "%s %s\n";
      printf(fmt, param[i].name, (char *) param[i].val);
    }
  }
  return;
}


static int get_tokens(char *s, char *tokens[]){
  char *tok;
  int i = 0;

  tok = strtok(s," ");
  while ( tok != NULL ){
    tokens[i++] = tok;
    tok = strtok(NULL," ");
  }

  return(i);

}




