#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "poisson.h"

#ifndef functions_h
#define functions_h

typedef struct{
    double Gr;
    double Pr;
    double l0;

    double dt;
    double t;

    int iter;

    int Nx;
    int Ny;
    double h;
    
    double *P;
    double *T;
    double *H_n_T;
    double *H_n_1_T;
    double *phi;
    double *u;
    double *v;
    double *u_star;
    double *v_star;
    double *H_n_u;
    double *H_n_v;
    double *H_n_1_u;
    double *H_n_1_v;

    Poisson_data *data;
} problem_struct;


problem_struct *create_problem(int Nx);
void free_problem(problem_struct *problem);

void print_mesh(problem_struct *problem);
void problem_to_file(problem_struct *problem);
void T_to_file(problem_struct *problem);
void U_to_file(problem_struct *problem);
void V_to_file(problem_struct *problem, const char* name);
void U_star_to_file(problem_struct *problem);
void V_star_to_file(problem_struct *problem);
void P_to_file(problem_struct *problem);

void compute_BC(problem_struct *problem);
void compute_v_star(problem_struct *problem);
void compute_v(problem_struct *problem);
void compute_T(problem_struct *problem);
void compute_H_n(problem_struct *problem);
void add_phi_P(problem_struct *problem);
void copy_H_n(problem_struct *problem);

void iterate(problem_struct *problem);


void test_poisson();

#endif
