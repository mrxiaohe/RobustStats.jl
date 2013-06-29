
#include <math.h>

void akerd_C(double *pts, int npts, double *x, int n, double *alam, double hval, double *output) {
 	double temp[5];
 	double sqrt5=pow(5, 0.5);
 	int i, j;
 	for(i=0; i<npts; i++){
 		double epan_alam_hval=0.0;
 		for(j=0; j<n; j++){
 			temp[j]=(pts[i]-x[j])/(hval*alam[j]);
 			if(fabs(temp[j])<sqrt5){
 				epan_alam_hval+=(0.75*(1.0-0.2*temp[j]*temp[j])/sqrt5)/(alam[j]*hval);
 			}
 		}
 		output[i]=epan_alam_hval/npts;
 	}
}

void alam_C(double *fhat, int n, double gm, double aval, double *output){
	int i;
	double numerator, power;
	for(i=0; i<n; i++){
		numerator=fhat[i]/gm;
		power=0.0-aval;
		output[i]=pow(numerator, power);
	}
}
