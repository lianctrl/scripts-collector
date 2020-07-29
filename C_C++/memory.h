/* memory.h */

int   *AllocVecI (int size);
int *ReallocVecI (int *v, int size);
void FreeVecI (int *v);

int  **AllocMatI (int size1, int size2);
int **ReallocMatI (int **v, int size1, int size2);
void FreeMatI (int **v);

real  *AllocVecR (int size);
void FreeVecR (real *v);

real **AllocMatR (int size1, int size2);
void FreeMatR (real **v);

char *AllocVecC (int size);
void FreeVecC (char *v);

char **AllocMatC (int size1, int size2);
void FreeMatC (char **v);

real RandR (int *seed);
void RandVec3 (real *p, int *pSeed);
void ErrExit (char *s);
