#ifndef Vector_h
#define Vector_h
#include <vector>
#include "Complex.h"
using std::vector;

complex conj(complex A);
double VectorNorm(int const num,complex*const B);
double VectorNorm(complex*const B,int const num);
double VectorNorm(double*const b,int const num);

void ConjMatrix(complex**const a,complex** b,int const m,int const n);
void ConjMatrix(complex**const a,complex** b,int const m);

//void ConjVector(complex*const a,complex** b,int const m);
complex* ConjVector(complex*const a,int const m);

void Inner_product(complex**const a,complex*const b,complex* c,int const m,int const n);
double Inner_product(complex* B,int n);
double Inner_product(double *a,double* b,int n);
//complex* Inner_product(complex**const a,complex*const b,int const m,int const n);

void MulMatrix(complex**const a,complex**const b,complex **C,int const m,int const n);
void MulMatrix(complex**const a,complex**const b,complex **C,int const m);

complex Vectorneiji(int const num,complex*const A,complex*const B);
void Inner_product(complex** a, vector<complex>& b, complex* c, int m, int n);        //¸´Êýc=A*b;

#endif // Vector_h
