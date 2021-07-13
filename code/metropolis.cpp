#include "metropolis.h"

uniform_real_distribution<> rdist;
uniform_int_distribution<int> idist;
int accepted = 0, cur_accepted = 0;

/* Implement update method here:
int update(ranlux48& rnd){

}
*/

int main(int argc, char** argv){
	int iseed, nthermo, nskip, nmeas;
	double e;

	scanf("%d %d %d %d %lf %lf\n",&nthermo,&nskip,&nmeas,&iseed,&beta,&e);
	scanf("%d ",&D);
	L = new int[D];
	for(int i = 0; i < D; i++)
		scanf("%d",&L[i]); 

	printf("nthermo=%d nskip=%d nmeas=%d D=%d L=%d iseed=%d beta=%e e=%e\n",nthermo,nskip,nmeas,D,L,iseed,beta,e);

	ranlux48 rnd(iseed);
	initialize(rnd, e);	
	genX(rnd, e);
	rdist = uniform_real_distribution<>(0,1);
	idist = uniform_int_distribution<int>(0,2*M-1);
	
	printf("STARTING CONFIG: \n"); print();

	cur_accepted = 0;
	for(int i = 0; i < nthermo; i++){
		cur_accepted += update(rnd);
		genX(rnd,e);
		printf("THERMALIZING CONFIG: \n"); print();
		printf("THERMALIZING ACCEPTANCE RATE: %e\n", (double) cur_accepted/nskip/V/D);
		cur_accepted = 0;
	}

	cur_accepted = 0;
	for(int i = 0; i < nmeas; i++){
		for(int j = 0; j < nskip; j++) cur_accepted += update(rnd);
		genX(rnd,e);
		printf("U: "); print();
		
		printf("METROPOLIS ACCEPTANCE RATE: %e\n", (double) cur_accepted/nskip/V/D);
		cur_accepted = 0;
	}

	return 0;	 
}

void print() {
	for(int i = 0; i < V; i++){
		for(int j = 0; j < D; j++){
			for(int k = 0; k < 3; k++){
				for(int l = 0; l < 3; l++){
					printf("%e %e ",U[i*D+j](k,l).real(),U[i*D+j](k,l).imag());
				}
			}
		}
	}
	printf("\n");
}
