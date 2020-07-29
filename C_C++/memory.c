/* memory.c */

#include "header.h"

int *AllocVecI (int size) 
{
    int *v;

    v = NULL;
    if (size == 0) {
      v = NULL;
      return v;
    } else {
      v = (int *) malloc (size * sizeof (int));
      if (!v) { printf("Allocation Failed!\n"); exit(1); }
    }
    return (v - 1);
}

int *ReallocVecI (int *v, int size) 
{
    if (size == 0) {
      free (v + 1);
      return NULL;
    } else {
      v = (int *) realloc (&v[1], size * sizeof (int));
      if (!v) { printf("Reallocation Failed!\n"); exit(1); }
    }
    return (v - 1);
}
  
void FreeVecI (int *v) 
{
    if (v != NULL)
      free (v + 1);
}

int **AllocMatI (int size1, int size2) 
{
    int **v, k;

    v = NULL;
    if ((size1 == 0) || (size2 == 0)) {
      v = NULL;
      return v;
    } else {
      v = (int **) malloc (size1 * sizeof (int *));
      v[0] = (int *) malloc (size1 * size2 * sizeof (int)) - 1;
      for (k = 1; k < size1; k ++) v[k] = v[k - 1] + size2;
      if (!v) { printf("Allocation Failed!\n"); exit(1); } 
    }
    return (v - 1);
}

int **ReallocMatI (int **v, int size1, int size2)
{
    int k;

    if ((size1 == 0) || (size2 == 0)) {
      v = NULL;
      return v;
    } else {
      v = (int **) realloc (v, size1 * sizeof (int *));
      v[0] = (int *) realloc (v[0], size1 * size2 * sizeof (int)) - 1;
      for (k = 1; k < size1; k ++) v[k] = v[k - 1] + size2;
      if (!v) { printf("Allocation Failed!\n"); exit(1); } 
    }
    return (v - 1);
}

void FreeMatI (int **v) 
{
    if (v != NULL) {
      free ((v + 1)[0] + 1);
      free (v + 1);
    }
}
  

real *AllocVecR (int size) 
{
    real *v;

    v = NULL;
    if (size == 0) {
      v = NULL;
      return v;
    } else {
      v = (real *) malloc (size * sizeof (real));
      if (!v) { printf("Allocation Failed!\n"); exit(1); }
    }
    return (v - 1);
}

void FreeVecR (real *v) 
{
    if (v != NULL)
      free (v + 1);
}
  
real **AllocMatR (int size1, int size2) 
{
    real **v;
    int k;

    v = NULL;
    if ((size1 == 0) || (size2 == 0)) {
      v = NULL;
      return v;
    } else {
      v = (real **) malloc (size1 * sizeof (real *));
      v[0] = (real *) malloc (size1 * size2 * sizeof (real)) - 1;
      for (k = 1; k < size1; k ++) v[k] = v[k - 1] + size2;
      if (!v) { printf("Allocation Failed!\n"); exit(1); }
    }
    return (v - 1);
}

void FreeMatR (real **v) 
{
    if (v != NULL) {
      free ((v + 1)[0] + 1);
      free (v + 1);
    }
}

void ErrExit (char *s) 
{
    printf ("Error: %s\n", s);
    exit (0);
}

void *AllocMemory (unsigned nelem, unsigned elsize)
{
    void *p;

    p=NULL;
    if ((nelem==0) || (elsize==0))
      p=NULL;
    else
      {
        if ((p=calloc((size_t)nelem,(size_t)elsize))==NULL) {
          fprintf(stdout,"Allocation Failed!\n"); fflush(stdout);
          exit(1);
        }
      }
    return p;
}

char **AllocMatC (int size1, int size2) {
    char **v;
    int k;
    v = (char **) malloc (size1 * sizeof (char *));
    v = v - 1;
    v[1] = (char *) malloc (size1 * size2 * sizeof (char));
    if (!v[1]) { printf("Allocation Failed!\n"); exit(1); }
    for (k = 2; k <= size1; k ++) v[k] = v[k - 1] + size2;
    if (!v) { printf("Allocation Failed!\n"); exit(1); }
    return (v);
}

void FreeMatC (char **v)
{
    free (v[1]);
    free (v + 1);
}

