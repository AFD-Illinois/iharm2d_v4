/*---------------------------------------------------------------------------------

  PARAMETERS.C
 
  -Declare "param" struct
  -Read and set parameters  

-----------------------------------------------------------------------------------*/

#include "decs.h"
#include <ctype.h>

// O(n) dictionary for converting strings to pointers to global variables
struct paramtable_entry 
{
  char *key;
  void *data;
};

#define MAXPARAMS (1000)
static struct paramtable_entry paramtable[MAXPARAMS];
static int nparam = 0;
static int nparamset = 0;

// Function to set paramtable keys and variable pointers
void set_param(char *key, void *data)
{
  paramtable[nparam].key = key;
  paramtable[nparam].data = data;
  nparam++;
}

// Check if param in param.dat matches any paramtable key 
int get_param(char *key, void **data)
{
  int n = 0;
  while (strncmp(key, paramtable[n].key, strlen(key)) != 0)
  {
    n++;
    if (n >= nparam)
    {
      *data = NULL;
      return 0;
    }
  }
  *data = paramtable[n].data;
  return 1;
}

// Assign param names and pointers
void set_core_params()
{
  set_param("tf", &tf);
  set_param("dt", &dt);
#if METRIC == MINKOWSKI
  set_param("x1Min", &x1Min);
  set_param("x1Max", &x1Max);
  set_param("x2Min", &x2Min);
  set_param("x2Max", &x2Max);
#elif METRIC == MKS
  set_param("a", &a);
  set_param("hslope", &hslope);
#if DEREFINE_POLES
  set_param("poly_xt", &poly_xt);
  set_param("poly_alpha", &poly_alpha);
  set_param("mks_smooth", &mks_smooth);
#endif
  if (N2 < NG) hslope = 1.;
  set_param("Rout", &Rout);
#endif
  set_param("cour", &cour);
  set_param("gam", &gam);
#if ELECTRONS
  set_param("game", &game);
  set_param("gamp", &gamp);
  set_param("fel0", &fel0);
  set_param("tptemin", &tptemin);
  set_param("tptemax", &tptemax);
#endif
  set_param("DTd", &DTd);
  set_param("DTf", &DTf);
  set_param("DTl", &DTl);
  set_param("DTr", &DTr);
  set_param("DTp", &DTp);
}

// Set runtime parameters from param.dat
void read_params(char* pfname)
{
  void *ptr;
  
  // Try reading the parameter file
  FILE *fp = fopen(pfname, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Cannot open parameter file: %s\n", pfname);
    exit(-1);
  }
  
  // Loop through lines in the parameter file
  char line[STRLEN];
  while (fgets(line, STRLEN, fp))
  {
    // Ignore comments, newlines and leading whitespace
    if (line[0] == '#' || line[0] == '\n' || isspace(line[0]))
      continue;
    
    // Check if param is not whitespace or null character
    char test[STRLEN], key[STRLEN];
    test[0] = '\0';
    sscanf(line, "%*s %s %*s %s", key, test);
    char *word = test;
    while (isspace(*word))
    {
      word++;
    }
    if (word[0] == '\0')
      continue;

    // Read in parameter depending on datatype
    char type[6];
    strncpy(type, line, 5);
    type[5] = 0;
    if (get_param(key, &ptr))
    {
      if (!strncmp(type, "[int]", 5))
      {
        int buf;
        sscanf(line, "%*s %s %*s %d", key, &buf);
        *((int*)ptr) = buf;
        nparamset++;
      } else if (!strncmp(type, "[dbl]", 5))
      {
        double buf;
        sscanf(line, "%*s %s %*s %lf", key, &buf);
        *((double*)ptr) = buf;
        nparamset++;
      } else if (!strncmp(type, "[str]", 5)) 
      {
        char buf[STRLEN];
        sscanf(line, "%*s %s %*s %s", key, buf);
        strcpy((char*)ptr, buf);
        nparamset++;
      }
    }
  }

  if (nparamset != nparam)
  {
    fprintf(stderr, "Set %i parameters, needed %i!\n", nparamset, nparam);
    exit(-1);
  }

  fclose(fp);
  fprintf(stdout, "Parameter file read\n\n");
}
