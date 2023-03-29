#include "functions.h"

problem_struct *create_problem(int Nx){
    problem_struct *problem = malloc(sizeof(problem_struct));
    problem->mesh = create_mesh(Nx);
    problem->Gr = 1e10;
    problem->Pr = 4.0;
    problem->dt = 0.01;
    return problem;
}

void free_problem(problem_struct *problem){
    free_mesh(problem->mesh);
    free(problem);
}

mesh_struct *create_mesh(int Nx){
    mesh_struct *mesh = malloc(sizeof(mesh_struct));
    mesh->h = 1.0/(Nx-1);
    mesh->Nx = Nx;
    mesh->Ny = 1.5/mesh->h + 1;
    mesh->P = calloc((mesh->Nx-1) * (mesh->Ny - 1), sizeof(double));
    mesh->T = calloc((mesh->Nx-1) * (mesh->Ny - 1), sizeof(double));
    mesh->u = calloc((mesh->Nx) * (mesh->Ny - 1), sizeof(double));
    mesh->v = calloc((mesh->Nx - 1) * (mesh->Ny), sizeof(double));
    mesh->phi = calloc((mesh->Nx) * (mesh->Ny), sizeof(double));
    int i;
    for(i = 0; i < mesh->Nx * mesh->Ny; i++){
        mesh->phi[i] = 1.0;
    }
    return mesh;
}

void free_mesh(mesh_struct *mesh){
    free(mesh->P);
    free(mesh->u);
    free(mesh->v);
    free(mesh->T);
    free(mesh);
}

void print_mesh(mesh_struct *mesh){
    int i,j;
    printf("Nx : %d, Ny : %d, h : %f\n", mesh->Nx, mesh->Ny, mesh->h);
    for(i = 0; i < mesh->Nx - 1; i++){
        for(j = 0; j < mesh->Ny - 1; j++){
            printf("P[%d, %d] = %f\n", i, j, mesh->P[(mesh->Nx-1) * i + j]);
        }
    }
}

void add_phi_P(problem_struct *problem){
    int i,j;
    mesh_struct *mesh = problem->mesh;
    for(i = 0; i < mesh->Nx - 1; i++){
        for(j = 0; j < mesh->Ny - 1; j++){
            mesh->P[(mesh->Ny-1) * i + j] += 1.0/4.0 * (mesh->phi[(mesh->Ny-1) * i + j] + mesh->phi[(mesh->Ny-1) * i + j + 1] + mesh->phi[(mesh->Ny-1) * (i+1) + j] + mesh->phi[(mesh->Ny-1) * (i+1) + j + 1]);
        }
    }
}