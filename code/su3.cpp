#include "su3.h"

#include <stdlib.h>
#include <math.h>
#include <iostream>

int D, *L, V, M;
double beta;
Matrix3cd *U, *x;

Matrix2cd s_x = (Matrix2cd() << cd(0,0), cd(1,0), cd(1,0), cd(0,0)).finished();
Matrix2cd s_y = (Matrix2cd() << cd(0,0), cd(0,-1), cd(0,1), cd(0,0)).finished();
Matrix2cd s_z = (Matrix2cd() << cd(1,0), cd(0,0), cd(0,0), cd(-1,0)).finished();

int step(int i, int d, int s) {
        int under = 1;
	for (int i = 0; i < d; i++) {
		under *= L[i];
	}
        return (under*L[d])*(i/(under*L[d])) + (i+under*s+abs(s)*under*L[d])%(under*L[d]); 
}

void plaquette(Matrix3cd *g, int i, int d1, int d2) {
        *g = Matrix3cd::Identity();
        *g = *g*U[i*D+d1];
        *g = *g*U[step(i,d1,1)*D+d2];
        *g = *g*(U[step(i,d2,1)*D+d1].inverse());
        *g = *g*(U[i*D+d2].inverse());
} 

void polyakov(Matrix3cd *g, int i) {
        *g = Matrix3cd::Identity();
        int i0 = i;
        do {
                *g = *g*U[i*D+0];
                i = step(i,0,1);
        } while (i0 != i);
}

void genX(ranlux48& rnd, double e){
	x = new Matrix3cd[2*M];
	for(int i = 0; i < M; i++){
		su3(&x[i], e, rnd);
		x[i+M] = x[i].conjugate().eval().transpose();	
	}
}

void initialize(ranlux48& rnd, double e){
	V = 1;
	for(int i = 0; i < D; i++) V *= L[i];
	M = 3*V;
	U = new Matrix3cd[V*D];
	for(int i = 0; i < V*D; i++){
		su3(&U[i], e, rnd);
	}
}

void su3(Matrix3cd *g, double e, ranlux48& rnd){
	uniform_real_distribution<> rdist(-0.5,0.5);
	Matrix2cd r, s, t;
	Matrix3cd R, S, T;
	
	su2(&r, e, rdist(rnd), rdist(rnd), rdist(rnd), rdist(rnd));
	su2(&s, e, rdist(rnd), rdist(rnd), rdist(rnd), rdist(rnd));
	su2(&t, e, rdist(rnd), rdist(rnd), rdist(rnd), rdist(rnd));

	R.block<2,2>(0,0) = r; R(2,2) = 1;
	S(0,0) = s(0,0); S(0,2) = s(0,1); S(1,1) = cd(1,0); S(2,0) = s(1,0); S(2,2) = s(1,1);
	T.block<2,2>(1,1) = t; T(0,0) = 1;

	*g = R*S*T; 
}

void su2(Matrix2cd *ret, double e, double r0, double r1, double r2, double r3){
	double mag = sqrt(r1*r1 + r2*r2 + r3*r3);
	*ret = Matrix2cd::Identity();
	*ret *= abs(r0)/r0*sqrt(1-e*e); 
	*ret += e/mag*cd(0,1)*(r1*s_x + r2*s_y + r3*s_z);
}
