#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include "AA_class.h"
//#include "EuclideanNorm.h"
using namespace std;

/*struct rod_t{
	double docks[2];
	//string shape;
	//double size[2];
};
struct con_t{
	double docks[6];
	//string shape;
	//double size;
};*/

struct AA_param_struct{
	double N[2];
	double R_comm;
	double max_comm_iter;
	double dock_const;
	double rod_docks[2];
	double con_docks[6];
	
};

/*double getEuclideanNorm(int input[], int size){
	double sum = 0;
	int i;
	for (i = 0; i < size; i++){
		sum += (input[i] * input[i]);
	}
	return sqrt(sum);
}*/
//extern void AA_bid(double *c, double *p, int p_size_input, double *p_old, int *j, int *m, int *count);


double** auction(double **current_pos, double **xf, int xf_cols, int *m, int m_length, double **p, int p_empty, int *j, int j_empty, AA_param_struct AA_param, double *nDocks, double *angDocks){
	int i;
	int k;
	int w;
	double** final_pos = 0;
	final_pos = new double*[3];

      for (int i = 0; i < 3; i++)
      {
            final_pos[i] = new double[6];

            for (int w = 0; w < 6; w++)
            {
                  // fill in some initial values
                  // (filling in zeros would be more logic, but this is just for the example)
                  final_pos[i][w] = 0;
            }
      }
	//double final_pos[3][6];
	int N = AA_param.N[1]+AA_param.N[0];
	double R_comm = AA_param.R_comm;
	double max_comm_iter = AA_param.max_comm_iter;
	int M = 100;
	double p_old[6][6];
	
	for (i = 0; i < m_length; i++){
		if(m[i]<M){
			M = m[i];
		}
	}
	if(xf_cols < M){
		M = xf_cols;
	}
	if(p_empty == 1){
		for (k = 0; k<6; k++){
			for(w = 0; w<6; w++){
				p[k][w] = 0;
				p_old[k][w] = -1;
			}
		}
		//p[6][6] = {{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0}};
		//double p_old[6][6] = {{-1, -1, -1, -1, -1, -1},{-1, -1, -1, -1, -1, -1},{-1, -1, -1, -1, -1, -1},{-1, -1, -1, -1, -1, -1},{-1, -1, -1, -1, -1, -1},{-1, -1, -1, -1, -1, -1}};
	}
	else{
		
		double p_old[6][6];
		for (int i = 0; i < N; ++i){
			for (int k = 0; k < M; ++k){
				p_old[i][k] = p[i][k];
			}
		}
	}
	if(j_empty == 1){
		for(i = 0; i<6; i++){
			j[i] = 1;
		}
	}

	double c[6][6] = {{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0}};
	double p_temp[6][6] = {{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0}};
	double done[6]= {0, 0, 0, 0, 0, 0};
	int count[6] = {0, 0, 0, 0, 0, 0};
	int maxDocks;
	int h;
	double angPenalty;
	double dockPenalty;
	double docks[6] = {0, 0, 0, 0, 0, 0};
	for (i = 0; i < N; i++){
		if(i <= AA_param.N[0]){
			//double docks[6];
			for (h = 0; h<6; h++){
				docks[h] = AA_param.con_docks[h];
			}
			maxDocks = 6; //Intentionally made length(docks)=6, account for this later.
		}
		else{
			//double docks[2];
			for (h = 0; h<2; h++){
				docks[h] = AA_param.rod_docks[h];
			}
			maxDocks = 2;
		}
		double angs[3];
		double norm_input[3];
		for(k = 0; k < M; k++){
			for (h = 0; h < 6; h++){
				angs[h] = angDocks[k] - docks[h];
			}
			int isAny = 0;
			for (h = 0; h < 6; h++){
				if(fabs(angs[i])<0.0001){
					isAny = 1;
					break;
				}
			}
			if(isAny == 1){
				angPenalty = 0;
			}
			else{
				angPenalty = 0; //changed this because didn't understand
			}
			dockPenalty = maxDocks - nDocks[k] + 1;
			if (dockPenalty < 0){
				dockPenalty = 0;
			}
			
			for (h = 0; h < 3; h++){
				norm_input[h] = current_pos[h][i] - xf[h][k];
				double sum = 0;
				int o;
				for (o = 0; o < 2; o++){
					sum += (norm_input[o] * norm_input[o]);
				}
				double EuclideanNorm = sqrt(sum);
				c[i][k] = 1000 * pow(EuclideanNorm,2) - dockPenalty - angPenalty;
			}
		}
	}

	
	while(1){
		int anyNotDone = 0;
		for (i = 0; i < 6; i++){
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
				/*int p_input[N] = p[i];
				int c_input[N] = c[i];
				int p_old_input[N] = c[i];
				int j_input = j[i];
				int m_input = m[i];
				int m_new_input = m_new[i];
				int count_input = count[i];*/
				AA_class AA_val;
				AA_val.AA_bid(c[i], p[i], 6, p_old[i], &j[i], &m[i], &count[i]);
				/*p[i] = p_input[N];
				j[i] = j_input;
				m[i] = m_input;
				m_new[i] = m_new_input;
				count[i] = count_input;*/
				int l;
				for (l = 0; l < 6; l++){
					p_old[i][l] = p[i][l];
				}
				done[i] = (count[i]>max_comm_iter);
			}
		}

		
		double p_max;
		int r;
		double Neighbors[6];
		int y; 
		for(i = 0; i < N; i++){	
			for (r = 0; r < N; r++){
				Neighbors[r] = (pow((current_pos[1][i] - current_pos[1][r]),0.5) + pow((current_pos[2][i] - current_pos[2][r]),2)) <= R_comm;
			}
			
			for(r = 0; r < 6; r++){
				p_max = -100;
				if(Neighbors[r]!=0){
					for (y = 0; y < 6; y++){
						if(p[y][r] > p_max){
							p_max = p[y][r];
						}
					}
				}
				p_temp[y][r] = p_max;
			}
			for(r = 0; r < 6; r++){
				for (y = 0; y < 6; y++){
					p[r][y] = p_temp[r][y];
				}
			}
			//p = p_temp;
		}
		//double final_pos[3][6];
		
		
		
	
		for (i = 0; i < 6; i++){
			int col_use = j[i];
			for(r = 0; r < 3; r++){
				final_pos[r][i] = xf[r][col_use];
			}
		}
	}
	return final_pos;	
}


int main(){
	int h;
	int w;

	double** current_pos = 0;
      current_pos = new double*[3];

      for (h = 0; h < 3; h++)
      {
            current_pos[h] = new double[6];
      }

      current_pos[0][0] = -0.6867;
      current_pos[0][1] =1.4183;
      current_pos[0][2] =-1.0614;
      current_pos[0][3] =1.4172;
      current_pos[0][4] =0.9309;
      current_pos[0][5] =-0.2425;

      current_pos[1][0] =0.1453;
      current_pos[1][1] =1.4412;
      current_pos[1][2] =1.4588;
      current_pos[1][3] =-0.0453;
      current_pos[1][4] =-1.1102;
      current_pos[1][5] =1.2888;

      current_pos[2][0] =1.8360;
      current_pos[2][1] =2.8871;
      current_pos[2][2] =0.9785;
      current_pos[2][3] =-2.9172;
      current_pos[2][4] =2.1936;
      current_pos[2][5] =2.7269;

      double** xf = 0;
      xf = new double*[3];

      for (h = 0; h < 3; h++)
      {
            xf[h] = new double[6];
      }
      xf[0][0] = 0;
      xf[0][1] = 0.3464;
      xf[0][2] = 0.3464;
      xf[0][3] = 0.3464;
      xf[0][4] = 0.1732;
      xf[0][5] = 0.1732;

      xf[1][0] = 0;
      xf[1][1] = 0.2000;
      xf[1][2] = -0.2000;
      xf[1][3] = 0;
      xf[1][4] = 0.1000;
      xf[1][5] = -0.1000;

      xf[2][0] = 0;
      xf[2][1] = 0;
      xf[2][2] = 0;
      xf[2][3] = 0;
      xf[2][4] = 0;
      xf[2][5] = 0;

     double** p = 0;
      p = new double*[6];

      for (h = 0; h < 6; h++)
      {
            p[h] = new double[6];

            for (w = 0; w < 6; w++)
            {
                  // fill in some initial values
                  // (filling in zeros would be more logic, but this is just for the example)
                  p[h][w] = 0;
            }
      }

	/*double current_pos[3][6] = {{-0.6867, 1.4183, -1.0614, 1.4172, 0.9309, -0.2425}, 
		{0.1453, 1.4412, 1.4588, -0.0453, -1.1102, 1.2888},
		{1.8360, 2.8871, 0.9785, -2.9172, 2.1936, 2.7269}};
	double xf[3][6] = {{0, 0.3464, 0.3464, 0.3464, 0.1732, 0.1732}, 
		{0, 0.2000, -0.2000, 0, 0.1000, -0.1000}, {0, 0, 0, 0, 0, 0}};*/
	
	//double *m;
	//m = new double[6];
	int m[6] = {6, 6, 6, 6, 6, 6};
	//double p[6][6];
	int j[6];
	AA_param_struct AA_param; 
	AA_param.N[0]=3;
	AA_param.N[1]=3;
	AA_param.R_comm = 0.75;
	AA_param.max_comm_iter = 12;
	AA_param.dock_const = 2.5;
	AA_param.rod_docks[0] = 2;
	AA_param.rod_docks[1] = 2;
	AA_param.con_docks[0] = 2;
	AA_param.con_docks[1] = 2;

	

	//double *nDocks;
	//nDocks = new double[6];
	double nDocks[6] = {2,2,2,2,2,2};
	//double *angDocks;
	//angDocks = new double[6];
	double angDocks[6] = {1.0472, 1.0472, 1.0472, 3.1416, 3.1416, 3.1416};
	//double** final_pos = auction(current_pos,xf, 6, m, p, 1, j, 1, AA_param, nDocks, angDocks);
	printf("Hello\n");
	double** final_pos = auction(current_pos, xf, 6, m, 6, p, 1, j, 1, AA_param, nDocks, angDocks);
	
	
	for (h = 0; h<3; h++){
		for (w = 0; w<6; w++){
			printf("final_pos: %f\n", final_pos[h][w]);
		}
	}
	//printf("final_pos: %f", final_pos[5][5]);
	return 0;
}



