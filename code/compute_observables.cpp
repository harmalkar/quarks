#include "compute_observables.h"
#include <iostream>

enum {
	PLAQUETTE,
	POLYAKOV,
	NOBS
};

double meas[NOBS][NMAX];
double mean[NOBS], stdev[NOBS];
double re_val, im_val;

int main(int argc, char *argv[]) {
	/* Seed the PRNG. */
	srand(time(NULL) + getpid());

	scanf("%*s\nD=%d",&D);
	L = new int[D];
	for(int i = 0; i < D; i++)
		scanf(" L%*d=%d",&L[i]);
	
	V = 1;
	for (int i = 0; i < D; i++)
		V *= L[i];

	U = new Matrix3cd[V*D];
	Matrix3cd res;

	int N = 0;
	while(read_configuration()) {
		meas[PLAQUETTE][N] = 0.;
		for (int i = 0; i < V; i++)
			for (int d = 0; d < D; d++)
				for (int dp = 0; dp < d; dp++){
					plaquette(&res, i, d, dp);
					meas[PLAQUETTE][N] += res.trace().real() / D / (D-1) * 2. / V;
				}
		for (int j = 0; j < V/L[0]; j++) {
			int k = j;
			int i = 0;
			for (int d = 1; d < D; d++) {
				i = step(i, d, k%L[0]);
				k /= L[0];
			}
			polyakov(&res, i);
			meas[POLYAKOV][N] += res.trace().real() / (V/L[0]);
		}
		N++;
	}

	for (int n = 0; n < NOBS; n++)
		bootstrap(&mean[n], &stdev[n], meas[n], N);

	printf("%d", N);
	for (int n = 0; n < NOBS; n++)
		printf(" %g %g", mean[n], stdev[n]);
	printf("\n");

	return 0;
}

bool read_configuration() {
	char c, c1, c2, c3;
	
	// skip all lines until gauge configuration is found 

	while(true){
		if(feof(stdin)) return 0;
		if((fscanf(stdin, "%c%c%c", &c1, &c2, &c3) == 3) && c1 == 'U' && c2 == ':' && c3 == ' ') {
			break;
		} else {
			do{
				c = fgetc(stdin);
			}
			while (c != '\n' && !feof(stdin));
		}
	}

	// read gauge configuration
	
	for (int i = 0; i < V*D; i++)
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++){
				if (fscanf(stdin, "%lf %lf ", &re_val, &im_val) != 2) return 0;
				U[i](j,k) = cd(re_val, im_val);
			}
	return 1;
}

// TODO we need blocking -- maybe even done automagically

void bootstrap(double *mean, double *stdev, double *dat, int len) {
	*mean = 0;
	for (int i = 0; i < len; i++)
		*mean += dat[i]/len;
	double var = 0.;
	for (int b = 0; b < B; b++) {
		double meanp = 0.;
		for (int i = 0; i < len; i++)
			meanp += dat[rand()%len]/len;
		var += (meanp-*mean)*(meanp-*mean)/B;
	}
	*stdev = sqrt(var);
}
