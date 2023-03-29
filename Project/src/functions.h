#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef functions_h
#define functions_h

typedef struct{
    int Nx;
    int Ny;
    double h;
    
    double *P;
    double *v;
    double *u;
    double *T;
    double *phi;
}mesh_struct;

typedef struct{
    double Gr;
    double Pr;

    double dt;

    mesh_struct *mesh;
} problem_struct;


problem_struct *create_problem(int Nx);
void free_problem(problem_struct *problem);

mesh_struct *create_mesh(int Nx);
void free_mesh(mesh_struct *mesh);
void print_mesh(mesh_struct *mesh);

void add_phi_P(problem_struct *problem);

#endif