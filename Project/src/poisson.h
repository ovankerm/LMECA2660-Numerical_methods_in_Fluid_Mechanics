#ifndef _POISSON_H_
#define _POISSON_H_

/*To include in the file in which you will call initialize_poisson_solver and poisson_solver*/

#include <petsc.h>
#include <petscsys.h>

//Structure storing petsc vectors

typedef struct {

	Vec b;
	Vec x;
	Mat A;
	KSP sles;

} Poisson_data;

PetscErrorCode initialize_poisson_solver(Poisson_data* data, int Nx, int Ny);
void poisson_solver(Poisson_data *data, double *phi, double *u_star, double *v_star, double h, double dt, int Nx, int Ny);
void free_poisson_solver(Poisson_data* data);

#endif

