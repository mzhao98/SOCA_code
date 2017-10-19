#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "AA_class.h"
using namespace std;



void AA_class::AA_bid(double *c, double *p, int p_size_input, double *p_old, int *j, int *m, int *count){
	double eps = 0.001;
	int i;
	//int p_size = sizeof(p)/sizeof(p[0]);
	int p_size = p_size_input;
	if(fabs(p[*j]) > p_old[*j]){
		int max = *m;
		int length = 0;
		for(i = 0; i < p_size; i++){
			if (fabs(p[i]) > 0 && i > *m){
				max = i;
			}
			if (p[i] > 0){
				length += 1;
			}
		}
		*m = max;
		if (length >= *m)
		{
			*m = length + 1;
			for (i = 0; i < *m; i++){
				p[i] = -1 * (fabs(p[i]) + eps);
			}
		}
		double v = 100000;

		for (i = 0; i < *m; i++){
			if((fabs(p[i]) + c[i]) < v){
				v = fabs(p[i]) + c[i];
				*j = i;
			}
		}
		double w = 100000;
		for (i = 0; i < *m; i++){
			if(i != *j){
				if(fabs(p[i]+c[i]) < w){
					w = fabs(p[i]+c[i]);
				}
			}
		}
		double gamma = w - v + eps;
		p[*j] = fabs(p[*j]) + gamma;
		*count = 0;
	}
	
	else {
		int foundAny = 0;
		for(i = 0; i < p_size; i++){
			if(p[i] != p_old[i]){
				foundAny = 1;
				break;
			}
		}
		//printf("foundAny = %d\n", foundAny);
		if(foundAny == 1){
			*count = 0;
			int max = *m;
			for(i = 0; i < p_size; i++){
				if (fabs(p[i]) > 0 && i > *m){
					max = i;
				}
			}
			*m = max;
		}
		else {
			*count += 1;
		}
	}
}