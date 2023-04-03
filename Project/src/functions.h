#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef functions_h
#define functions_h

typedef struct{
    double Gr;
    double Pr;
    double l0;

    double dt;

    int Nx;
    int Ny;
    double h;
    
    double *P;
    double *T;
    double *H_n_T;
    double *phi;
    double *u;
    double *v;
    double *u_star;
    double *v_star;
    double *H_n_u;
    double *H_n_v;
    double *H_n_1_u;
    double *H_n_1_v;
} problem_struct;


problem_struct *create_problem(int Nx);
void free_problem(problem_struct *problem);
void print_mesh(problem_struct *problem);

void compute_BC(problem_struct *problem);
void compute_v_star(problem_struct *problem);
void compute_v(problem_struct *problem);
void compute_T(problem_struct *problem);
void compute_H_n(problem_struct *problem);

void add_phi_P(problem_struct *problem);

#endif