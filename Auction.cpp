#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cmath>
using namespace std;

struct AA_param_struct{
	double N;
	double R_comm;
	double max_comm_iter;
};

int getEuclideanNorm(int input[], int size){
	int sum = 0;
	int i;
	for (i = 0; i < size; i++){
		sum += (input[i] * input[i]);
	}
	return sqrt(sum);
}


int* auction(double *current_pos, double *xf, int xf_cols, int xf_rows, int *m, double *p, int *j, AA_param_struct AA_param){
	double N = AA_param.N;
	double R_comm = AA_param.R_comm;
	double max_comm_iter = AA_param.max_comm_iter;
	int M = xf_cols;
	int i;
	int k;
	if(p == NULL){
		p = new int[N][M];
		p_old = new int[N][M];
		for (i = 0; i < N; i++){
			for (k = 0; i < M; i++){
				p[i][k] = 0;
				p_old[i][k] = -1;
			}
		}

	}
	else{
		p_old = p;
	}
	if(j == NULL){
		j = new int[N];
		for (i = 0; i < N; i++){
			j[i] = 1;
		}
	}

	c = new int[N][M];
	p_temp = new int[N][M];
	done = new int[N];
	count = new int[N];
	m_new = new int[N];

	for (i = 0; i < N; i++){
		done[i] = 0;
		count[i] = 0;
		m_new[i] = 0;
		for (k = 0; i < M; i++){
			c[i][k] = 0;
			p_temp[i][k] = -1;
		}
	}
	int norm_input[2];
	int l;
	for (i = 0; i < N; i++){
		for (k = 0; i < M; i++){
			for(l = 0; l < xf_rows; l++){
				norm_input[l] = current_pos[l][i] - xf[l][k];
			}
			c[i][k] = pow(getEuclideanNorm(norm_input),2);
		}
	}
	int done[3] = {0, 0, 0};
	while(1){
		int anyNotDone = 0;
		for (i = 0; i < 3; i++){
			if (done[i]==0){
				anyNotDone = 1;
				break;
			}
		}
		if(anyNotDone == 0){
			break;
		}
		for(i = 0; i < N; i++){
			if (done[i]==0){
				int p_input[N] = p[i];
				int c_input[N] = c[i];
				int p_old_input[N] = c[i];
				int j_input = j[i];
				int m_input = m[i];
				int m_new_input = m_new[i];
				int count_input = count[i];
				AA_bid(c_input, p_input, p_old_input, j_input, m_input, m_new_input, count_input);
				p[i] = p_input[N];
				j[i] = j_input;
				m[i] = m_input;
				m_new[i] = m_new_input;
				count[i] = count_input;

				p_old[i] = p[i];
				done[i] = (count[i]>max_comm_iter);
			}
		}

		for(i = 0; i < N; i++){
			int r;
			int Neighbors[N] = new int[N];
			for (r = 0; r < N; r++){
				Neighbors[r] = pow((current_pos[1][i] - current_pos[1][r]),2);
			}



		}

	}




	return xf;	

}

int main(){
	//double* current_pos = malloc(sizeof(double)*2*3);
	double current_pos[2][3] = {1, 2, 3}, {2, 3, 4};
	double xf[2][3] = {2, 3, 4}, {3, 4, 5};
	int m[3] = {2, 2, 2};
	double p[3][3] = {0.003, 0.005, -0.001}, {0.003, 0.005, -0.001}, {0.003, 0.005, -0.001};
	int j[3] = {1, 1, 1};
	AA_param_struct AA_param; 
	AA_param.N = 3;
	AA_param.R_comm = 10;
	AA_param.max_comm_iter = 5;
	double* final_pos = Auction(current_pos,xf, 3, 2, m,p,j,AA_param);
	int i;
	int j;
	for(i = 0; i < 2; i++){
		for(j = 0; j < 3; i++){
			printf("final_pos %d: %f\n", i, final_pos[i][j]);
		}
	}
	for(i = 0; i < 3; i++){
		printf("m %d: %f\n", i, m[i]);
	}
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; i++){
			printf("p %d: %f\n", i, p[i][j]);
		}
	}
}



