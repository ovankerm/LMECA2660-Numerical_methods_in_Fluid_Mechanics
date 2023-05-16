# ifndef STRUCT_H
# define STRUCT_H

#include "poisson.h"

# define GR 10000000000.0
# define SQRT_GR 100000.0
# define PR 4.0
# define PR_M1 0.25
# define L0 0.001
# define DT_DTAU 1000
# define OMEGA_S 0.1

typedef struct{
    int nx;
    int ny;

    double h;
    
    int nxu;
    int nyu;
    double *u_data;
    double **u;
    double **u_star;
    double **H_u;
    double **H_1_u;
    int nxu_low;
    int nxu_high;
    int nyu_low;
    int nyu_high;
    int n_u_points;
    int *u_i;
    int *u_j;

    int nxv;
    int nyv;
    double *v_data;
    double **v;
    double **v_star;
    double **H_v;
    double **H_1_v;
    int nxv_low;
    int nxv_high;
    int nyv_low;
    int nyv_high;
    int n_v_points;
    int *v_i;
    int *v_j;

    int nxP;
    int nyP;
    double *P_data;
    double **P;

    double *phi_data;
    double **phi;

    int nxT;
    int nyT;
    double *T_data;
    double **T;
    double **T_1;
    double **H_T;
    double **H_1_T;
    int nxT_low;
    int nxT_high;
    int nyT_low;
    int nyT_high;
    int n_T_points;
    int *T_i;
    int *T_j;
} data_sim;

typedef struct{
    data_sim *sim_data;

    double t;
    double dt;
    int iter;

    double Reh;
    double Rehw;

    Poisson_data *poiss_data;

    FILE *diag_file;
} problem_struct;

# endif