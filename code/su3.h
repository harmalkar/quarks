#ifndef SU3_H
#define SU3_H

#include <Eigen/Dense>
#include <stdlib.h>

using namespace Eigen;
using namespace std;

typedef complex<double> cd;

extern int D, *L, V, M;
extern double beta;
extern Matrix3cd *U, *x;

int step(int i, int d, int s);
void plaquette(Matrix3cd *g, int i, int d1, int d2);
void polyakov(Matrix3cd *g, int i);
void genX(ranlux48& rnd, double e);
void initialize(ranlux48& rnd, double e);
void su3(Matrix3cd *g, double e, ranlux48& rnd);
void su2(Matrix2cd *ret, double e, double r0, double r1, double r2
, double r3);

#endif /* ndef SU3_H */
