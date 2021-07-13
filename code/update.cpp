#include "metropolis.h"

int update(ranlux48& rnd) {
	Matrix3cd updater, neighbors;
	int accepted = 0;
	for(int i = 0; i < V; i++) {
		for(int mu = 0; mu < D; mu++) {
			neighbors = Matrix3cd::Zero();
			for(int nu = 0; nu < D; nu++){
				if(nu == mu) continue;
				neighbors += U[step(i,mu,1)*D+nu]*(U[step(i,nu,1)*D+mu].inverse())*(U[i*D+nu].inverse());
				neighbors += U[step(step(i,mu,1),nu,-1)*D+nu].inverse()*(U[step(i,nu,-1)*D+mu].inverse())*U[step(i,nu,-1)*D+nu];	
			}
			updater = x[idist(rnd)];
			double prob = exp(beta/3*((updater-Matrix3cd::Identity())*U[i*D+mu]*neighbors).trace().real());
			if(rdist(rnd) < prob){
				accepted++;
				U[i*D+mu] = updater*U[i*D+mu];
			}
		}
	}
	return accepted;
}


