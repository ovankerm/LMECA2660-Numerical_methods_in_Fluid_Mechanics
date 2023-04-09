# ifndef STRUCT_H
# define STRUCT_H

#include "poisson.h"

# define GR 10000000000.0
# define SQRT_GR 0.00001
# define PR 4.0
# define PR_M1 0.25
# define L0 0.001

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

    int nxv;
    int nyv;
    double *v_data;
    double **v;
    double **v_star;
    double **H_v;
    double **H_1_v;

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
    double **H_T;
    double **H_1_T;
} data_sim;

typedef struct{
    data_sim *sim_data;

    double t;
    double dt;
    int iter;

    Poisson_data *poiss_data;
} problem_struct;

# endif