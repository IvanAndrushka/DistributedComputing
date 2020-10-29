#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005

typedef struct particle_t{
	double x;
	double y;
	double vx;
	double vy;
	double m;
} particle_t;

typedef struct cell{
	double x;
	double y;
	double m;
} cell;

void init_particles(long seed, long ncside, long long n_part, particle_t *par)
{
    long long i;

    srandom(seed);

    for(i = 0; i < n_part; i++)
    {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;

        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}

int main(int argc, char *argv[]){ // seed random generator 1; size of grid (nxn) 2; num particles 3; num time-stamps 4

	int seed = atoi(argv[1]);
	int ncside = atoi(argv[2]);
	int n_part = atoi(argv[3]);
	int time_stamps = atoi(argv[4]);

	int aux_x, aux_y, i, j, t, p, linha, coluna, flag_coluna, flag_linha;
	double Fax, Fay, d_2, x_final=0, y_final=0, m_final=0, d_x, d_y;

	particle_t *par = calloc(n_part, sizeof(particle_t));

	cell *grid = (cell*) malloc(ncside * ncside * sizeof(cell));

	//inicializacao de particulas
	init_particles(seed, ncside, n_part, par);

	for(t = 0; t < time_stamps; t++)
	{
		#pragma omp parallel
		{
			//inicializacao da grid
			#pragma omp for
			for(i = 0; i < ncside * ncside; i++)
			{
				grid[i].x = 0;
				grid[i].y = 0;
				grid[i].m = 0;
			}

			//calculo do somatorio para a media do centro de massa
			#pragma omp for
			for(i = 0; i < n_part; i++)
			{
				aux_x = floor(par[i].x*ncside);
				aux_y = floor(par[i].y*ncside);

				#pragma omp atomic
				grid[aux_x * ncside + aux_y].m += par[i].m;
				#pragma omp atomic
				grid[aux_x * ncside + aux_y].x += par[i].x*par[i].m;
				#pragma omp atomic
				grid[aux_x * ncside + aux_y].y += par[i].y*par[i].m;
			}

			//calculo do centro de massa
			#pragma omp for
			for(i = 0; i < ncside * ncside; i++)
			{
					if(grid[i].m == 0)
						continue;

					grid[i].x = grid[i].x / grid[i].m;
					grid[i].y = grid[i].y / grid[i].m;
			}

			//calculo das forcas e posicoes
			#pragma omp for private(Fax,Fay, i, j, flag_linha, flag_coluna, linha, coluna, d_y, d_x, d_2)
			for(p = 0; p < n_part; p++){
				Fax = 0;
				Fay = 0;

				//percorer linhas das celulas adjacentes
				for(i=0, linha = floor(par[p].x*ncside) -1; i < 3; i++, linha++){
					flag_linha = 0;

					//celula adjacente fora da grid
					if(linha < 0){
						flag_linha = 1;
						linha = ncside - 1;
					}
					else if(linha >= ncside){
						flag_linha = 2;
						linha = 0;
					}

					//percorrer colunas das celulas adjacentes
					for(j = 0, coluna = floor(par[p].y*ncside) - 1; j < 3; j++, coluna++){
						flag_coluna=0;

						//celula adjacente fora da grid
						if(coluna < 0){
							flag_coluna = 1;
							coluna= ncside - 1;
						}
						else if(coluna >= ncside){
							flag_coluna = 2;
							coluna = 0;
						}

						//calculo da distancia entre celula e centro de massa
						if(flag_linha == 1)
							d_y = (grid[linha * ncside + coluna].y - par[p].y) - 1;
						else if(flag_linha == 2)
							d_y = 1 - (grid[linha * ncside + coluna].y - par[p].y);
						else
							d_y = grid[linha * ncside + coluna].y - par[p].y;

						if(flag_coluna == 1)
							d_x = (grid[linha * ncside + coluna].x - par[p].x) - 1;
						else if(flag_coluna == 2)
							d_x = 1 - (grid[linha * ncside + coluna].x - par[p].x);
						else
							d_x = grid[linha * ncside + coluna].x - par[p].x;

						//calcular da distancia
						d_2 = sqrt(d_x * d_x + d_y * d_y);

						//calculo da força resultante
						if(d_2 >= EPSLON){
							Fax += G*(par[p].m * grid[linha * ncside + coluna].m)*(d_x)/(d_2 * d_2 * d_2);
							Fay += G*(par[p].m * grid[linha * ncside + coluna].m)*(d_y)/(d_2 * d_2 * d_2);
						}
					}
				}

				//calcular aceleraçao
				Fax = Fax / par[p].m;
				Fay = Fay / par[p].m;

				//calcular velocidade
				par[p].vx = par[p].vx + Fax;
				par[p].vy = par[p].vy + Fay;

				//calcular posiçao
				par[p].x = par[p].x + par[p].vx + Fax/2;
				while(par[p].x >= 1)
					par[p].x = par[p].x - 1;
				while(par[p].x < 0)
					par[p].x = par[p].x + 1;

				par[p].y = par[p].y + par[p].vy + Fay/2;
				while(par[p].y >= 1)
					par[p].y = par[p].y - 1;
				while(par[p].y < 0)
					par[p].y = par[p].y + 1;

			}
			//printf("time:%d\n", t);
		}
	}

	//print da particula 0
	printf("%.2f %.2f\n", par[0].x, par[0].y);

	//calculo do centro de massa total
	#pragma omp parallel for reduction(+:m_final, x_final, y_final)
	for(i=0; i<n_part; i++)
	{
		m_final += par[i].m;
		x_final += par[i].x*par[i].m;
		y_final += par[i].y*par[i].m;
	}

	x_final = x_final / m_final;
	y_final = y_final / m_final;

	//print do centro de massa total
	printf("%.2f %.2f\n", x_final, y_final);

	//free da memoria alocada
	free(par);
	free(grid);

	exit(0);
}
