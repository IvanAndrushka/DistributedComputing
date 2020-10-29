#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include <mpi.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005
#define N_SIZE n_part*2/p

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )

typedef struct particle_t {
  double x;
  double y;
  double vx;
  double vy;
  double m;
  bool tag;
  bool send;
} particle_t;

typedef struct cell {
  double x;
  double y;
  double m;
} cell;

typedef struct Matrix {
  int id;
  int col;
  int row;
  double i;
  double j;
} Matrix;

void init_particles(long seed, long ncside, long long n_part, long long n_part_total, particle_t *par)
{
  long long i;

  for (i = 0; i < n_part; i++)
  {
    par[i].x = RND0_1;
    par[i].y = RND0_1;
    par[i].vx = RND0_1 / ncside / 10.0;
    par[i].vy = RND0_1 / ncside / 10.0;

    par[i].m = RND0_1 * ncside / (G * 1e6 * n_part_total);
    par[i].tag = false;
    par[i].send = false;
  }
}

int get_index(int size, int x, int y) {
  return size * y + x;
}

int get_proc(int i, int j, int col, int row, Matrix *matrix, int p, int std_col, int std_row) {
  if (i < 0)
    i += row;
  else if (i >= row)
    i -= row;

  if (j < 0)
    j += col;
  else if (j >= col)
    j -= col;

  for (int k = 0; k < p; k++) {
    if (matrix[k].j / std_col == j && matrix[k].i / std_row == i)
      return k;
  }
  return -1;
}

void centro_calc( cell *in, cell *inout, int *len, MPI_Datatype *type ){
    inout->x += in->x;
    inout->y += in->y;
    inout->m += in->m;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) { // seed random generator 1; size of grid (nxn) 2; num particles 3; num time-stamps 4

  if (argc < 4) {
    printf("not enough parameters\n");
    exit(0);
  }

  int seed = atoi(argv[1]); srandom(seed);
  int ncside = atoi(argv[2]);
  int n_part = atoi(argv[3]), n_part_total = n_part;
  int time_stamps = atoi(argv[4]);
  int id, p, counter = 0, counter_aux=0;
  int col_p, row_p, std_col, std_row;
  int col_row[2];
  cell * grid;

  Matrix *matrix;

  int aux_x, aux_y;
  double Fax, Fay, d_2, d_x, d_y, flag_linha, flag_coluna;

  // init MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &p);

  if(p > ncside*ncside){
    p = ncside*ncside;
  }

  if(id < p){
  // enviar tamanhos das matrizes
    matrix = (Matrix *) malloc(p * sizeof(Matrix));
    if(matrix == NULL){
      printf("%d-- error: malloc matrix\n", id);fflush(stdout);
      MPI_Finalize();
      exit(1);
    }

    col_p = sqrt(p);
    while(p % col_p){
      col_p--;
    }
    int ncside_aux = ncside;
    while (ncside_aux % col_p) {
      ncside_aux--;
    }
    row_p = p / col_p;
    col_row[0] = ncside_aux / col_p;
    col_row[1] = ncside_aux / row_p;
    std_col = col_row[0];
    std_row = col_row[1];

    for (int i = 0, proc = 0; i < row_p; i++) {
      if (i == row_p - 1) {
        col_row[1] += ncside % row_p;
      }

      for (int j = 0; j < col_p; j++) {

        if (j == col_p - 1) {
          col_row[0] += ncside % col_p;
        }

        matrix[proc].id = proc;
        matrix[proc].col = col_row[0];
        matrix[proc].row = col_row[1];
        matrix[proc].i = i * std_row;
        matrix[proc].j = j * std_col;
        proc++;
      }
      col_row[0]-=ncside % col_p;
    }

    //inicializacao de particulas
    particle_t *par = (particle_t *)malloc(N_SIZE * sizeof(particle_t));
    if(par == NULL){
      printf("%d-- error: malloc par\n", id);fflush(stdout);
      MPI_Finalize();
      exit(1);
    }
    particle_t *send_par = (particle_t *)malloc((n_part / p + n_part % p) * sizeof(particle_t));
    if(send_par == NULL){
      printf("%d-- error: malloc send_par\n", id);fflush(stdout);
      MPI_Finalize();
      exit(1);
    }

    int send_npart = n_part/p;

    for (int i = 0; i < p; i++) {
      if(i == p-1){
        send_npart+= n_part % p;
      }
      if (id == 0)
      {
        init_particles(seed, ncside, send_npart,n_part, send_par);
        if(i == 0){
          send_par[0].tag = true;
        }
      }
      MPI_Bcast(send_par, send_npart * sizeof(particle_t), MPI_BYTE, 0, MPI_COMM_WORLD);
      //guardar proprias particulas
      for (int j = 0; j < send_npart; j++) {
        if ((send_par[j].x >= (matrix[id].j) / ((double)ncside)) && (send_par[j].y >= (matrix[id].i) / ((double)ncside)) && (send_par[j].x < (matrix[id].j) / ((double)ncside) + (matrix[id].col) / ((double)ncside)) && (send_par[j].y < (matrix[id].i) / ((double)ncside) + (matrix[id].row) / ((double)ncside))) {
          par[counter] = send_par[j];
          counter++;
        }
      }
    }

    //alocar matriz centros de massa
    grid = (cell*) malloc((matrix[id].col + 2) * (matrix[id].row + 2) * sizeof(cell));
    if(grid == NULL){
      printf("%d-- error: malloc grid\n", id);fflush(stdout);
      MPI_Finalize();
      exit(1);
    }
    MPI_Status status;
    int dest[8];
    cell* grid_send = (cell*) malloc(max(matrix[id].row, matrix[id].col) * sizeof(cell));
    if(grid_send == NULL){
      printf("%d-- error: malloc grid_send\n", id);fflush(stdout);
      MPI_Finalize();
      exit(1);
    }

    // determinar os processos de destino
    dest[0] = get_proc(matrix[id].i / std_row - 1, matrix[id].j / std_col - 1, col_p, row_p, matrix, p, std_col, std_row);
    dest[1] = get_proc(matrix[id].i / std_row - 1, matrix[id].j / std_col, col_p, row_p, matrix, p, std_col, std_row);
    dest[2] = get_proc(matrix[id].i / std_row - 1, matrix[id].j / std_col + 1, col_p, row_p, matrix, p, std_col, std_row);
    dest[3] = get_proc(matrix[id].i / std_row, matrix[id].j / std_col - 1, col_p, row_p, matrix, p, std_col, std_row);
    dest[4] = get_proc(matrix[id].i / std_row, matrix[id].j / std_col + 1, col_p, row_p, matrix, p, std_col, std_row);
    dest[5] = get_proc(matrix[id].i / std_row + 1, matrix[id].j / std_col - 1, col_p, row_p, matrix, p, std_col, std_row);
    dest[6] = get_proc(matrix[id].i / std_row + 1, matrix[id].j / std_col, col_p, row_p, matrix, p, std_col, std_row);
    dest[7] = get_proc(matrix[id].i / std_row + 1, matrix[id].j / std_col + 1, col_p, row_p, matrix, p, std_col, std_row);

    //enviar particulas que escaparam
    particle_t **esc_par = (particle_t **)malloc(8 * sizeof(particle_t*));
    if(esc_par == NULL){
      printf("%d-- error: malloc sec_par\n", id);fflush(stdout);
      MPI_Finalize();
      exit(1);
    }

    for(int i=0; i<8; i++){
      esc_par[i] = (particle_t *)malloc((N_SIZE)*sizeof(particle_t));
      if(esc_par[i] == NULL){
        printf("%d-- error: malloc esc_par[%d]\n", id, i);fflush(stdout);
        MPI_Finalize();
        exit(1);
      }
    }

    int esc_counter[8] = {0};

    for (int t = 0; t < time_stamps; t++) {
      //inicializacao da grid
      counter_aux = 0;

      for (int i = 0; i < ((matrix[id].col + 2) * (matrix[id].row + 2)); i++)
      {
        grid[i].x = 0;
        grid[i].y = 0;
        grid[i].m = 0;
      }
      //calculo centro de massa
      for (int i = 0; i < counter; i++) {
        if ((par[i].x >= (matrix[id].j /(double) ncside)) && (par[i].y >= (matrix[id].i /(double) ncside)) && (par[i].x < (matrix[id].j /(double) ncside + matrix[id].col /(double) ncside)) && (par[i].y < (matrix[id].i /(double) ncside + matrix[id].row /(double) ncside))) {
          aux_x = floor(par[i].x * ncside - matrix[id].j);
          aux_y = floor(par[i].y * ncside - matrix[id].i);

          grid[(aux_y + 1) * (matrix[id].col + 2) + 1 + aux_x].m += par[i].m;
          grid[(aux_y + 1) * (matrix[id].col + 2) + 1 + aux_x].x += par[i].x * par[i].m;
          grid[(aux_y + 1) * (matrix[id].col + 2) + 1 + aux_x].y += par[i].y * par[i].m;
        }
      }

      for (int i = 0; i < (matrix[id].col + 2) * (matrix[id].row + 2); i++)
      {
        if (grid[i].m != 0) {
          grid[i].x = grid[i].x / grid[i].m;
          grid[i].y = grid[i].y / grid[i].m;
        }
      }

      // envio e receção das vizinhanças

      /**********************************************/
      //enviar 1 vizinhanca - canto superior esquerdo
      if(dest[0]!= id){
        grid_send[0] = grid[get_index(matrix[id].col + 2, 1, 1)];
        MPI_Send(grid_send, sizeof(cell) * 1, MPI_BYTE, dest[0], 1, MPI_COMM_WORLD);

        // Recebe canto superior esquerdo e mete-o no inferior direito
        MPI_Recv(grid_send, sizeof(cell) * 1, MPI_BYTE, dest[7], 1, MPI_COMM_WORLD, &status);
        grid[(matrix[id].col + 2) * (matrix[id].row + 2) - 1] = grid_send[0];
      }
      else{
        grid[(matrix[id].col + 2) * (matrix[id].row + 2) - 1] = grid[get_index(matrix[id].col + 2, 1, 1)];
      }


      /**********************************************/
      //enviar 2 vizinhanca - linha de cima
      if(dest[1]!=id){
        for (int i = 1; i <= matrix[id].col; i++)
          grid_send[i-1] = grid[get_index(matrix[id].col + 2, i, 1)];
        MPI_Send(grid_send, matrix[id].col * sizeof(cell), MPI_BYTE, dest[1], 2, MPI_COMM_WORLD);

        // Recebe linha de cima e mete-a na linha de baixo
        MPI_Recv(grid_send, matrix[id].col * sizeof(cell), MPI_BYTE, dest[6], 2, MPI_COMM_WORLD, &status);
        for (int i = 1; i <= matrix[id].col; i++)
          grid[get_index(matrix[id].col + 2, i, matrix[id].row + 1)] = grid_send[i - 1];
      }
      else {
        for (int i = 1; i <= matrix[id].col; i++)
          grid[get_index(matrix[id].col + 2, i, matrix[id].row + 1)] = grid[get_index(matrix[id].col + 2, i, 1)];
      }

      /***********************************************/
      //enviar 3 vizinhanca - canto superior direito
      if(dest[2]!=id){
        grid_send[0] = grid[get_index(matrix[id].col + 2, matrix[id].row, 1)];
        MPI_Send(grid_send, sizeof(cell), MPI_BYTE, dest[2], 3, MPI_COMM_WORLD);

        // Recebe canto superior direito e mete no canto inferior esquerdo
        MPI_Recv(grid_send, sizeof(cell) * 1, MPI_BYTE, dest[5], 3, MPI_COMM_WORLD, &status);
        grid[(matrix[id].col + 2) * (matrix[id].row + 1)] = grid_send[0];
      }
      else {
        grid[(matrix[id].col + 2) * (matrix[id].row + 1)] = grid[get_index(matrix[id].col + 2, matrix[id].row, 1)];
      }

      /*****************************************************/
      //enviar 4 vizinhanca - coluna esquerda
      if(dest[4]!=id){
        for (int i = 1; i <= matrix[id].row; i++)
          grid_send[i - 1] = grid[get_index(matrix[id].col + 2, 1, i)];
        MPI_Send(grid_send, matrix[id].row * sizeof(cell), MPI_BYTE, dest[3], 4, MPI_COMM_WORLD);

        // Recebe coluna esquerda e mete na coluna direita
        MPI_Recv(grid_send, matrix[id].row * sizeof(cell), MPI_BYTE, dest[4], 4, MPI_COMM_WORLD, &status);
        for (int i = 1; i <= matrix[id].row; i++)
          grid[get_index(matrix[id].col + 2, matrix[id].col + 1, i)] = grid_send[i - 1];
      }
      else {
        for (int i = 1; i <= matrix[id].row; i++)
          grid[get_index(matrix[id].col + 2, matrix[id].col + 1, i)] = grid[get_index(matrix[id].col + 2, 1, i)];
      }
      /********************************************************/
      //enviar 5 vizinhanca - coluna direita
      if(dest[5]!=id){
        for (int i = 1; i <= matrix[id].row; i++)
          grid_send[i - 1] = grid[get_index(matrix[id].col + 2, matrix[id].col, i)];
        MPI_Send(grid_send, matrix[id].row * sizeof(cell), MPI_BYTE, dest[4], 5, MPI_COMM_WORLD);

        // Recebe coluna direita e mete na coluna esquerda
        MPI_Recv(grid_send, matrix[id].row * sizeof(cell), MPI_BYTE, dest[3], 5, MPI_COMM_WORLD, &status);
        for (int i = 1; i <= matrix[id].row; i++)
          grid[get_index(matrix[id].col + 2, 0, i)] = grid_send[i - 1];
        }
      else {
        for (int i = 1; i <= matrix[id].row; i++)
          grid[get_index(matrix[id].col + 2, 0, i)] = grid[get_index(matrix[id].col + 2, matrix[id].col, i)];
      }

      /***********************************************************/
      //enviar 6 vizinhanca - canto inferior esquerdo
      if(dest[5]!=id){
        grid_send[0] = grid[get_index(matrix[id].col + 2, 1, matrix[id].col)];
        MPI_Send(grid_send, sizeof(cell), MPI_BYTE, dest[5], 6, MPI_COMM_WORLD);

        // Recebe canto inferior esquerdo e guarda-o no canto superior direito
        MPI_Recv(grid_send, sizeof(cell) * 1, MPI_BYTE, dest[2], 6, MPI_COMM_WORLD, &status);
        grid[matrix[id].col + 1] = grid_send[0];
      }
      else {
        grid[matrix[id].col + 1] = grid[get_index(matrix[id].col + 2, 1, matrix[id].col)];     
      }

      /**********************************************************/
      //enviar 7 vizinhanca - linha de baixo
      if(dest[6]!=id){
        for (int i = 1; i <= matrix[id].col; i++)
          grid_send[i - 1] = grid[get_index(matrix[id].col + 2, i, matrix[id].row)];
        MPI_Send(grid_send, matrix[id].col * sizeof(cell), MPI_BYTE, dest[6], 7, MPI_COMM_WORLD);

        // Recebe linha de baixo e guarda-o na linha de cima
        MPI_Recv(grid_send, sizeof(cell)*matrix[id].col, MPI_BYTE, dest[1], 7, MPI_COMM_WORLD, &status);
        for (int i = 1; i <= matrix[id].col; i++)
          grid[get_index(matrix[id].col + 2, i, 0)] = grid_send[i - 1];
      }
      else {
        for (int i = 1; i <= matrix[id].col; i++)
          grid[get_index(matrix[id].col + 2, i, 0)] = grid[get_index(matrix[id].col + 2, i, matrix[id].row)];
      }

      /***********************************************************/
      //enviar 8 vizinhanca - canto inferior direito
      if(dest[7]!=id){
        grid_send[0] = grid[get_index(matrix[id].col + 2, matrix[id].col, matrix[id].row)];
        MPI_Send(grid_send, sizeof(cell), MPI_BYTE, dest[7], 8, MPI_COMM_WORLD);

        // Recebe canto inferior direito e guarda-o no canto superior esquerdo
        MPI_Recv(grid_send, sizeof(cell) * 1, MPI_BYTE, dest[0], 8, MPI_COMM_WORLD, &status);
        grid[0] = grid_send[0];
      }
      else {
        grid[0] = grid[get_index(matrix[id].col + 2, matrix[id].col, matrix[id].row)];
      }

      // calculo das forças das particulas
      //calculo de forças
      for (int part = 0; part < counter; part++) {
        if ((par[part].x >= (matrix[id].j /(double) ncside)) && (par[part].y >= (matrix[id].i /(double) ncside)) && (par[part].x < (matrix[id].j /(double) ncside + matrix[id].col /(double) ncside)) && (par[part].y < (matrix[id].i /(double) ncside + matrix[id].row /(double) ncside))) {

          Fax = 0;
          Fay = 0;

          for (int i = 0, linha = floor(par[part].y * ncside - matrix[id].i)-1; i < 3; i++, linha++) {

            for (int j = 0, coluna = floor(par[part].x * ncside - matrix[id].j)-1; j < 3; j++, coluna++) {
 
              d_x = grid[(matrix[id].col+2)+(matrix[id].col+2)*linha +1 + coluna].x - par[part].x;
              d_y = grid[(matrix[id].col+2)+(matrix[id].col+2)*linha +1 + coluna].y - par[part].y;

              if (d_x > 0.5) d_x = 1 - d_x;
              if (d_y > 0.5) d_y = 1 - d_y;
              if (d_x < -0.5) d_x = 1 + d_x;
              if (d_y < -0.5) d_y = 1 + d_y;

              d_2 = sqrt(d_x * d_x + d_y * d_y);
              //calculo da força resultante
              if (d_2 >= EPSLON) {
                Fax += G * (par[part].m * grid[(matrix[id].col+2)+(matrix[id].col+2)*linha +1 + coluna].m) * (d_x) / (d_2 * d_2 * d_2);
                Fay += G * (par[part].m * grid[(matrix[id].col+2)+(matrix[id].col+2)*linha +1 + coluna].m) * (d_y) / (d_2 * d_2 * d_2);
              }
            }
          }

          //calcular aceleraçao
          Fax = Fax / par[part].m;
          Fay = Fay / par[part].m;

          //calcular velocidade
          par[part].vx = par[part].vx + Fax;
          par[part].vy = par[part].vy + Fay;

          par[part].x = par[part].x + par[part].vx + Fax / 2;
          if (par[part].x >= 1)
            par[part].x = par[part].x - 1;
          if (par[part].x < 0)
            par[part].x = par[part].x + 1;

          par[part].y = par[part].y + par[part].vy + Fay / 2;
          if (par[part].y >= 1)
            par[part].y = par[part].y - 1;
          if (par[part].y < 0)
            par[part].y = par[part].y + 1;

          
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);

      for (int i=0; i<8; i++){
        esc_counter[i] = 0;
      }

      for (int part = 0; part < counter; part++) {
        if(par[part].send==true){
          continue;
        }
        //superior esquerdo
        if ((par[part].x < (matrix[id].j /(double) ncside)) && (par[part].y < (matrix[id].i /(double) ncside))) {
          par[part].send = true;
          esc_par[0][esc_counter[0]] = par[part];
          esc_counter[0]++;
        }
        //superior direito
        else if ((par[part].y < (matrix[id].i /(double) ncside)) && (par[part].x >= (matrix[id].j /(double) ncside + matrix[id].col / (double)ncside))) {
          par[part].send = true;
          esc_par[2][esc_counter[2]] = par[part];
          esc_counter[2]++;
        }
        //inferior esquerdo
        else if ((par[part].x < (matrix[id].j /(double) ncside)) && (par[part].y >= (matrix[id].i /(double) ncside + matrix[id].row / (double)ncside))) {
          par[part].send = true;
          esc_par[5][esc_counter[5]] = par[part];
          esc_counter[5]++;
        }
        //inferior direito
        else if ((par[part].x >= (matrix[id].j /(double) ncside+ matrix[id].col / (double)ncside)) && (par[part].y >= (matrix[id].i /(double) ncside + matrix[id].row / (double)ncside))) {
          par[part].send = true;
          esc_par[7][esc_counter[7]] = par[part];
          esc_counter[7]++;
        }
        //cima
        else if (par[part].y < (matrix[id].i /(double) ncside)) {
          par[part].send = true;
          esc_par[1][esc_counter[1]] = par[part];
          esc_counter[1]++;
        }
        //coluna esquerda
        else if (par[part].x < (matrix[id].j /(double) ncside)) {
          par[part].send = true;
          esc_par[3][esc_counter[3]] = par[part];
          esc_counter[3]++;
        }
        //coluna direita
        else if (par[part].x >= (matrix[id].j /(double) ncside + matrix[id].col / (double)ncside)) {
          par[part].send = true;
          esc_par[4][esc_counter[4]] = par[part];
          esc_counter[4]++;
        }
        //baixo
        else if (par[part].y >= (matrix[id].i /(double) ncside + matrix[id].row / (double)ncside)) {
          par[part].send = true;
          esc_par[6][esc_counter[6]] = par[part];
          esc_counter[6]++;
        }
      }
 
      for (int i = 0; i < 8; i++) {

        counter_aux = 0;
        MPI_Send(&esc_counter[i], 1, MPI_INT, dest[i], i, MPI_COMM_WORLD);
        MPI_Recv(&counter_aux, 1, MPI_INT, dest[7 - i], i, MPI_COMM_WORLD, &status);
        MPI_Barrier(MPI_COMM_WORLD);

        if (esc_counter[i]>0){
          MPI_Send(esc_par[i], esc_counter[i] * sizeof(particle_t), MPI_BYTE, dest[i], i, MPI_COMM_WORLD);
        }

        if(counter_aux > 0){
          MPI_Recv(send_par, counter_aux * sizeof(particle_t), MPI_BYTE, dest[7 - i], i, MPI_COMM_WORLD, &status);
          
          for(int j=0, h=0; h<counter_aux; j++){
            if(j<counter){
              if(par[j].send==true){
                par[j]=send_par[h];
                h++;
                par[j].send = false; 
              }
            }
            else{
              counter++;
              par[j]=send_par[h];
              h++;
              par[j].send = false;
            }
          }
        }
        MPI_Barrier(MPI_COMM_WORLD);
      }

      MPI_Barrier(MPI_COMM_WORLD);

    }

    //calculo centro de massa total
    for (int i = 0; i < ((matrix[id].col + 2) * (matrix[id].row + 2)); i++)
    {
      grid[i].x = 0;
      grid[i].y = 0;
      grid[i].m = 0;
    }

    for (int i = 0; i < counter; i++) {
      if ((par[i].x >= (matrix[id].j /(double) ncside)) && (par[i].y >= (matrix[id].i /(double) ncside)) && (par[i].x < (matrix[id].j /(double) ncside + matrix[id].col /(double) ncside)) && (par[i].y < (matrix[id].i /(double) ncside + matrix[id].row /(double) ncside))) {

        aux_x = floor(par[i].x * ncside - matrix[id].j);
        aux_y = floor(par[i].y * ncside - matrix[id].i);

        grid[(aux_y + 1) * (matrix[id].col + 2) + 1 + aux_x].m += par[i].m;
        grid[(aux_y + 1) * (matrix[id].col + 2) + 1 + aux_x].x += par[i].x * par[i].m;
        grid[(aux_y + 1) * (matrix[id].col + 2) + 1 + aux_x].y += par[i].y * par[i].m;
      }
    }

    grid[0].x = 0;
    grid[0].y = 0;
    grid[0].m = 0;

    for (int i = 1; i <= matrix[id].row; i++) {
      for (int j = 1; j <= matrix[id].col; j++) {
          grid[0].x += grid[i*(matrix[id].col + 2)+j].x;
          grid[0].y += grid[i*(matrix[id].col + 2)+j].y;
          grid[0].m += grid[i*(matrix[id].col + 2)+j].m;
      }
    }

    if(id == 0) {
      grid[1].x = 0;
      grid[1].y = 0;
      grid[1].m = 0;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Op final_cell;
    MPI_Op_create((MPI_User_function *) centro_calc, 1, &final_cell );
    MPI_Reduce(&grid[0], &grid[1], sizeof(cell), MPI_BYTE, final_cell, 0, MPI_COMM_WORLD);

    for(int i=0; i<counter; i++){
      if(par[i].tag==true){
        if ((par[i].x >= (matrix[id].j /(double) ncside)) && (par[i].y >= (matrix[id].i /(double) ncside)) && (par[i].x < (matrix[id].j / (double)ncside + matrix[id].col / (double) ncside)) && (par[i].y < (matrix[id].i /(double) ncside + matrix[id].row / (double) ncside))) {

          printf("%.2f %.2f\n",par[i].x, par[i].y);fflush(stdout);
          break;
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(id == 0) {
      grid[1].x /= grid[1].m;
      grid[1].y /= grid[1].m;

      printf("%.2f %.2f\n", grid[1].x, grid[1].y);fflush(stdout);
    }

    free(par);
    free(send_par);
    for(int i=0; i<8; i++)
      free(esc_par[i]);
    free(esc_par);
    free(grid);
    free(grid_send);
    free(matrix);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Op_free(&final_cell);
  }
  MPI_Finalize();

  exit(0);
}
