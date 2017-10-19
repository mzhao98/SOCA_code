//global double p[3];
//double j;
//int m_bool;
//int count;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "AA_bid.h"
using namespace std;

void AA_bid(double *c, double *p, int p_size_input, double *p_old, int *j, int *m, int *count){
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
		printf("foundAny = %d\n", foundAny);
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
/*
int main(){
	double c[] = {-0.6867, 1.4183, -1.0614, 1.4172, 0.9309, -0.2425};
	double p[] = {0, 0, 0, 0, 0, 0};
	double p_old[] = {-1, -1, -1, -1, -1, -1};
	int j = 1;
	int m = 2;
	
	int count = 0;
	AA_bid(c, p, 3, p_old, &j, &m, &count);

	int i;
	for(i = 0; i < 3; i++){
		printf("c %d: %f\n", i, c[i]);
	}
	for(i = 0; i < 3; i++){
		printf("p %d: %f\n", i, p[i]);
	}
	for(i = 0; i < 3; i++){
		printf("p_old %d: %f\n", i, p_old[i]);
	}
	printf("int j = %d\n", j);
	printf("int m = %d\n", m);
	
	printf("int count = %d\n", count);
}

*/
