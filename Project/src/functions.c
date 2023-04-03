#include "functions.h"

int P_IND(int i, int j, int Ny){
    return (Ny-1) * i + j;
}

int T_IND(int i, int j, int Ny){
    return (Ny+1) * i + j;
}

int U_IND(int i, int j, int Ny){
    return (Ny+1) * i + j;
}

int V_IND(int i, int j, int Ny){
    return Ny * i + j;
}

problem_struct *create_problem(int Nx){
    problem_struct *problem = malloc(sizeof(problem_struct));

    problem->Gr = 1e10;
    problem->Pr = 4.0;
    problem->l0 = 1e-3;
    problem->dt = 1;

    problem->h = 2.0/(3.0 * (Nx-1));
    problem->Nx = Nx;
    problem->Ny = 1.0/problem->h + 1;
    problem->P = calloc((problem->Nx - 1) * (problem->Ny - 1), sizeof(double));
    problem->T = calloc((problem->Nx + 1) * (problem->Ny + 1), sizeof(double));
    problem->H_n_T = calloc((problem->Nx + 1) * (problem->Ny + 1), sizeof(double));
    problem->phi = calloc((problem->Nx - 1) * (problem->Ny - 1), sizeof(double));
    problem->u = calloc((problem->Nx) * (problem->Ny + 1), sizeof(double));
    problem->v = calloc((problem->Nx + 1) * (problem->Ny), sizeof(double));
    problem->u_star = calloc((problem->Nx) * (problem->Ny + 1), sizeof(double));
    problem->v_star = calloc((problem->Nx + 1) * (problem->Ny), sizeof(double));
    problem->H_n_u = calloc((problem->Nx) * (problem->Ny + 1), sizeof(double));
    problem->H_n_v = calloc((problem->Nx + 1) * (problem->Ny), sizeof(double));
    problem->H_n_1_u = calloc((problem->Nx) * (problem->Ny + 1), sizeof(double));
    problem->H_n_1_v = calloc((problem->Nx + 1) * (problem->Ny), sizeof(double));
    // int i;
    // for(i = 0; i < (problem->Nx-1) * (problem->Ny-1); i++){
    //     problem->phi[i] = i;
    // }

    return problem;
}

void free_problem(problem_struct *problem){
    free(problem->P);
    free(problem->T);
    free(problem->H_n_T);
    free(problem->phi);
    free(problem->u);
    free(problem->v);
    free(problem->u_star);
    free(problem->v_star);
    free(problem->H_n_u);
    free(problem->H_n_v);
    free(problem->H_n_1_u);
    free(problem->H_n_1_v);
    free(problem);
}

void print_mesh(problem_struct *problem){
    int i,j;
    printf("Nx : %d, Ny : %d, h : %f\n", problem->Nx, problem->Ny, problem->h);
    printf("PRESSURE\n");
    for(i = 1; i < problem->Nx; i++){
        for(j = 1; j < problem->Ny; j++){
            printf("T[%d, %d] = %f\n", i, j, problem->T[T_IND(i, j, problem->Ny)]);
        }
    }
    printf("U\n");
    for(i = 0; i < problem->Nx; i++){
        for(j = 0; j < problem->Ny - 1; j++){
            printf("U[%d, %d] = %f\n", i, j, problem->u[U_IND(i, j, problem->Ny)]);
        }
    }
}

void compute_BC(problem_struct *problem){
    int Nx = problem->Nx, Ny = problem->Ny;
    int i;
    for(i = 1; i < problem->Nx - 1; i++){
        problem->u[U_IND(i, 0, Ny)] = -1.0/5.0 * (problem->u[U_IND(i, 3, Ny)] - 5.0 * problem->u[U_IND(i, 2, Ny)] + 15.0 * problem->u[U_IND(i, 1, Ny)]);
        problem->u[U_IND(i, Ny, Ny)] = problem->u[U_IND(i, Ny-1, Ny)];
        problem->T[T_IND(i, 0, Ny)] = problem->T[T_IND(i, 1, Ny)] + problem->h;
        problem->T[T_IND(i, Ny, Ny)] = problem->T[T_IND(i, Ny - 1, Ny)] * (problem->l0/problem->h - 0.5)/(problem->l0/problem->h + 0.5);
    }
    problem->T[T_IND(i, 0, Ny)] = problem->T[T_IND(i, 1, Ny)] + problem->h;
    problem->T[T_IND(i, Ny, Ny)] = problem->T[T_IND(i, Ny - 1, Ny)] * (problem->l0/problem->h - 0.5)/(problem->l0/problem->h + 0.5);
    for(i = 1; i < problem->Ny - 1; i++){
        problem->v[V_IND(0, i, Ny)] = -1.0/5.0 * (problem->v[V_IND(3, i, Ny)] - 5.0 * problem->v[V_IND(2, i, Ny)] + 15.0 * problem->v[V_IND(1, i, Ny)]);
        problem->v[V_IND(Nx, i, Ny)] = -1.0/5.0 * (problem->v[V_IND(Nx - 3, i, Ny)] - 5.0 * problem->v[V_IND(Nx - 2, i, Ny)] + 15.0 * problem->v[V_IND(Nx - 1, i, Ny)]);
        problem->T[T_IND(0, i, Ny)] = problem->T[T_IND(0, i, Ny)];
        problem->T[T_IND(Nx, i, Ny)] = problem->T[T_IND(Nx - 1, i, Ny)];
    }
    problem->T[T_IND(0, i, Ny)] = problem->T[T_IND(0, i, Ny)];
    problem->T[T_IND(Nx, i, Ny)] = problem->T[T_IND(Nx - 1, i, Ny)];
}

void compute_v_star(problem_struct *problem){
    int i, j;
    int Nx = problem->Nx, Ny = problem->Ny;
    double h = problem->h;
    double *u = problem->u;
    double *v = problem->v;
    double *P = problem->P;
    double *T = problem->T;
    double *H_n_u = problem->H_n_u;
    double *H_n_v = problem->H_n_v;
    for(i = 1; i < Nx - 1; i++){
        for(j = 1; j < Ny; j++){
            problem->u_star[U_IND(i, j, Ny)] = u[U_IND(i, j, Ny)] + problem->dt * (H_n_u[U_IND(i, j, Ny)]
                                            - 1.0/h * (P[P_IND(i, j - 1, Ny)] - P[P_IND(i - 1, j - 1, Ny)])
                                            + pow(problem->Gr, -0.5)/(h*h) * (u[U_IND(i + 1, j, Ny)] + u[U_IND(i - 1, j, Ny)] + u[U_IND(i, j + 1, Ny)] + u[U_IND(i, j - 1, Ny)] - 4 * u[U_IND(i, j, Ny)]));
        }
    }

    for(i = 1; i < Nx; i++){
        for(j = 1; j < Ny - 1; j++){
            problem->v_star[V_IND(i, j, Ny)] = v[V_IND(i, j, Ny)] + problem->dt * (H_n_v[V_IND(i, j, Ny)]
                                            - 1.0/h * (P[P_IND(i - 1, j, Ny)] - P[P_IND(i - 1, j - 1, Ny)])
                                            + pow(problem->Gr, -0.5)/(h*h) * (v[V_IND(i + 1, j, Ny)] + v[V_IND(i - 1, j, Ny)] + v[V_IND(i, j + 1, Ny)] + v[V_IND(i, j - 1, Ny)] - 4 * v[V_IND(i, j, Ny)])
                                            - 0.5 * (T[T_IND(i - 1, j, Ny)] + T[T_IND(i - 1, j - 1, Ny)]));
        }
    }
}

void compute_v(problem_struct *problem){
    int i, j;
    int Nx = problem->Nx, Ny = problem->Ny;
    double h = problem->h;
    double *u_star = problem->u_star;
    double *v_star = problem->v_star;
    double *phi = problem->phi;
    for(i = 1; i < Nx - 1; i++){
        for(j = 1; j < Ny; j++){
            problem->u[U_IND(i, j, Ny)] = u_star[U_IND(i, j, Ny)] - problem->dt/h * (phi[P_IND(i, j - 1, Ny)] - phi[P_IND(i - 1, j - 1, Ny)]);
        }
    }

    for(i = 1; i < Nx; i++){
        for(j = 1; j < Ny - 1; j++){
            problem->v[V_IND(i, j, Ny)] = v_star[V_IND(i, j, Ny)] - problem->dt/h * (phi[P_IND(i - 1, j, Ny)] - phi[P_IND(i - 1, j - 1, Ny)]);
        }
    }
}

void compute_T(problem_struct *problem){
    int i, j;
    int Nx = problem->Nx, Ny = problem->Ny;
    double h = problem->h;
    double *T = problem->T;
    for(i = 1; i < Nx; i++){
        for(j = 1; j < Ny; j++){
            T[T_IND(i, j, Ny)] += problem->dt * (problem->H_n_T[T_IND(i, j, Ny)]
                                + pow(problem->Gr, -0.5)/(problem->Pr * h*h) * (T[T_IND(i + 1, j, Ny)] + T[T_IND(i - 1, j, Ny)] + T[T_IND(i, j + 1, Ny)] + T[T_IND(i, j - 1, Ny)] - 4 * T[T_IND(i, j, Ny)]));
        }
    }
}

void compute_H_n(problem_struct *problem){
    int i, j;
    int Nx = problem->Nx, Ny = problem->Ny;
    double h = problem->h;
    double *u = problem->u;
    double *v = problem->v;
    double *T = problem->T;
    for(i = 1; i < Nx - 1; i++){
        for(j = 1; j < Ny; j++){
            problem->H_n_u[U_IND(i, j, Ny)] = 1/(4 * h) * ((u[U_IND(i, j, Ny)] + u[U_IND(i + 1, j, Ny)]) * (u[U_IND(i + 1, j, Ny)] - u[U_IND(i, j, Ny)])
                                                + (u[U_IND(i, j, Ny)] + u[U_IND(i - 1, j, Ny)]) * (u[U_IND(i, j, Ny)] - u[U_IND(i - 1, j, Ny)])
                                                + (v[V_IND(i - 1, j + 1, Ny)] + v[V_IND(i, j + 1, Ny)]) * (u[U_IND(i, j + 1, Ny)] - u[U_IND(i, j, Ny)])
                                                + (v[V_IND(i, j, Ny)] + v[V_IND(i - 1, j, Ny)]) * (u[U_IND(i, j, Ny)] - u[U_IND(i - 1, j, Ny)]));
        }
        for(j = 1; j < Ny - 1; j++){
            problem->H_n_T[T_IND(i, j, Ny)] = 1/(4 * h) * ((u[U_IND(i, j + 1, Ny)] + u[U_IND(i + 1, j + 1, Ny)]) * (T[T_IND(i + 1, j, Ny)] - T[T_IND(i - 1, j, Ny)])
                                                + (v[V_IND(i + 1, j, Ny)] + v[V_IND(i + 1, j + 1, Ny)]) * (T[T_IND(i, j + 1, Ny)] - T[T_IND(i, j - 1, Ny)]));
        }
    }

    for(i = 1; i < Nx; i++){
        for(j = 1; j < Ny - 1; j++){
            problem->H_n_v[V_IND(i, j, Ny)] = 1/(4 * h) * ((u[U_IND(i, j, Ny)] + u[U_IND(i, j + 1, Ny)]) * (v[V_IND(i + 1, j, Ny)] - v[V_IND(i, j, Ny)])
                                                + (u[U_IND(i - 1, j + 1, Ny)] + u[U_IND(i - 1, j, Ny)]) * (v[V_IND(i, j, Ny)] - v[V_IND(i - 1, j, Ny)])
                                                + (v[V_IND(i, j + 1, Ny)] + v[V_IND(i, j, Ny)]) * (v[V_IND(i, j + 1, Ny)] - v[V_IND(i, j, Ny)])
                                                + (v[V_IND(i, j, Ny)] + v[V_IND(i, j - 1, Ny)]) * (v[V_IND(i, j, Ny)] - v[V_IND(i - 1, j, Ny)]));
        }
    }

}

void add_phi_P(problem_struct *problem){
    int i,j;
    int Nx = problem->Nx, Ny = problem->Ny;
    for(i = 0; i < Nx - 1; i++){
        for(j = 0; j < Ny - 1; j++){
            problem->P[P_IND(i, j, Ny)] += problem->phi[P_IND(i, j, Ny)];
        }
    }
}

