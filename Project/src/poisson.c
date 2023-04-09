#include <mpi.h>
#include "poisson.h"

int PHI_IND(int i, int j, int Ny){
    return (Ny-1) * i + j;
}

int U_INDEX(int i, int j, int Ny){
    return (Ny+1) * i + j;
}

int V_INDEX(int i, int j, int Ny){
    return Ny * i + j;
}

/*Called by poisson_solver at each time step*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*  OKKKK  -Impose zero mass flow here by changing value of U_star*/
/*  OKKKK  -Fill vector rhs*/
void computeRHS(double *rhs, PetscInt rowStart, PetscInt rowEnd, double *u_star, double *v_star, double h, double dt, int Nx, int Ny)
{
    int i, j;
    for(i = 1; i < Nx; i++){
        v_star[V_INDEX(i, 0, Ny)] = 0.0;
        v_star[V_INDEX(i, Ny - 1, Ny)] = 0.0;
    }
    
    for(i = 1; i < Ny; i++){
        u_star[U_INDEX(0, i, Ny)] = 0.0;
        u_star[U_INDEX(Nx - 1, i, Ny)] = 0.0;
    }

    int r = rowStart;
    for(r=rowStart; r<rowEnd ; r++){
        i = (int) r/(Ny - 1);
        j = r%(Ny - 1);
		rhs[r] = 1/(dt * h) * (u_star[U_INDEX(i+1, j+1, Ny)] - u_star[U_INDEX(i, j+1, Ny)] + v_star[V_INDEX(i+1, j+1, Ny)] - v_star[V_INDEX(i+1, j, Ny)]); /*WRITE HERE (nabla dot u_star)/dt at each mesh point r*/
        /*Do not forget that the solution for the Poisson equation is defined within a constant.
        One point from Phi must then be set to an abritrary constant.*/
    }


}

/*To call at each time step after computation of U_star. This function solves the poisson equation*/
/*and copies the solution of the equation into your vector Phi*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*   OKKKKK - Change the call to computeRHS as you have to modify its prototype too*/
/*   OKKKKK - Copy solution of the equation into your vector PHI*/
void poisson_solver(Poisson_data *data, double *phi, double *u_star, double *v_star, double h, double dt, int Nx, int Ny)
{

    /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
    int its;
    PetscInt rowStart, rowEnd;
    PetscScalar *rhs, *sol;

    KSP sles = data->sles;
    Vec b = data->b;
    Vec x = data->x;

    /* Fill the right-hand-side vector : b */
    VecGetOwnershipRange(b, &rowStart, &rowEnd);
    VecGetArray(b, &rhs);
    computeRHS(rhs, rowStart, rowEnd, u_star, v_star, h, dt, Nx, Ny); /*MODIFY THE PROTOTYPE HERE*/
    VecRestoreArray(b, &rhs);

    // VecView(b, PETSC_VIEWER_STDOUT_WORLD);


    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);
    KSPGetIterationNumber(sles, &its);
    PetscPrintf(PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n", its);

    VecGetArray(x, &sol);

    int r;
    for(r=rowStart; r<rowEnd; r++){
        phi[r] = sol[r];
    }

    VecRestoreArray(x, &sol);
    free(sol);
    free(rhs);
}

/*This function is called only once during the simulation, i.e. in initialize_poisson_solver.*/
/*In its current state, it inserts unity on the main diagonal.*/
/*More than probably, you should need to add arguments to the prototype ... .*/
/*Modification to do in this function : */
/*  OKKKK -Insert the correct factor in matrix A*/
/*

Compute the discretized laplacian of phi
we need to set a value of phi, ex, phi[0] = 0

*/
void computeLaplacianMatrix(Mat A, int rowStart, int rowEnd, int Nx, int Ny)
{
    int i,j;
    int r;
    for(j = 1; j < Ny - 2; j++){
        r = rowStart + PHI_IND(0, j, Ny);
        MatSetValue(A, r, r , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r+1 , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r-1 , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r+(Ny-1) , 1.0, INSERT_VALUES);
        r = rowStart + PHI_IND(Nx-2, j, Ny);
        MatSetValue(A, r, r , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r+1 , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r-1 , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r-(Ny-1) , 1.0, INSERT_VALUES);
    }
    for(i = 1; i < Nx - 2; i++){
        r = rowStart + PHI_IND(i, 0, Ny);
        MatSetValue(A, r, r , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r+1 , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r+(Ny-1) , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r-(Ny-1) , 1.0, INSERT_VALUES);
        r = rowStart + PHI_IND(i, Ny - 2, Ny);
        MatSetValue(A, r, r , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r-1 , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r+(Ny-1) , 1.0, INSERT_VALUES);
        // MatSetValue(A, r, r-(Ny-1) , 1.0, INSERT_VALUES);
        for(j = 1; j < Ny - 2; j++){
            r = rowStart + PHI_IND(i, j, Ny);
            MatSetValue(A, r, r , -4.0, INSERT_VALUES);
            MatSetValue(A, r, r+1 , 1.0, INSERT_VALUES);
            MatSetValue(A, r, r-1 , 1.0, INSERT_VALUES);
            MatSetValue(A, r, r+(Ny-1) , 1.0, INSERT_VALUES);
            MatSetValue(A, r, r-(Ny-1) , 1.0, INSERT_VALUES);
        }
    }

    r = rowStart + PHI_IND(0, 0, Ny);
    MatSetValue(A, r, r , 1.0, INSERT_VALUES);
    // MatSetValue(A, r, r+1 , 1.0, INSERT_VALUES);
    // MatSetValue(A, r, r+(Ny-1) , 1.0, INSERT_VALUES);

    r = rowStart + PHI_IND(0, Ny - 2, Ny);
    MatSetValue(A, r, r , 1.0, INSERT_VALUES);
    // MatSetValue(A, r, r-1 , 1.0, INSERT_VALUES);
    // MatSetValue(A, r, r+(Ny-1) , 1.0, INSERT_VALUES);

    r = rowStart + PHI_IND(Nx - 2, 0, Ny);
    MatSetValue(A, r, r , 1.0, INSERT_VALUES);
    // MatSetValue(A, r, r+1 , 1.0, INSERT_VALUES);
    // MatSetValue(A, r, r-(Ny-1) , 1.0, INSERT_VALUES);

    r = rowStart + PHI_IND(Nx - 2, Ny-2, Ny);
    MatSetValue(A, r, r , 1.0, INSERT_VALUES);
    // MatSetValue(A, r, r-1 , 1.0, INSERT_VALUES);
    // MatSetValue(A, r, r-(Ny-1) , 1.0, INSERT_VALUES);


    // int r;
    // MatSetValue(A, 0, 0 , 1.0, INSERT_VALUES);
    // for (r = rowStart + 1; r < rowStart + Ny; r++){
    //     MatSetValue(A, r, r, -4.0, INSERT_VALUES);
    //     MatSetValue(A, r, r+1 , 1.0, INSERT_VALUES);
    //     MatSetValue(A, r, r-1 , 1.0, INSERT_VALUES);
    //     MatSetValue(A, r, r+(Ny-1) , 1.0, INSERT_VALUES);
    // }
    // for (;r < rowEnd - Ny; r++){
    //     MatSetValue(A, r, r, -4.0, INSERT_VALUES);
    //     MatSetValue(A, r, r+1 , 1.0, INSERT_VALUES);
    //     MatSetValue(A, r, r-1 , 1.0, INSERT_VALUES);
    //     MatSetValue(A, r, r+(Ny-1) , 1.0, INSERT_VALUES);
    //     MatSetValue(A, r, r-(Ny-1) , 1.0, INSERT_VALUES);
    // }
    // MatSetValue(A, 1, 1 , -4.0, INSERT_VALUES);
    // MatSetValue(A, 1, 2 , 1.0, INSERT_VALUES);
    // MatSetValue(A, 1, 0 , 1.0, INSERT_VALUES);
    // for (r = rowStart + 2; r < rowEnd; r++){
    //     MatSetValue(A, r, r, -4.0, INSERT_VALUES);
    //     /*USING MATSETVALUE FUNCTION, INSERT THE GOOD FACTOR AT THE GOOD PLACE*/
    //     /*Be careful; the solution from the system solved is defined within a constant.
    //     One point from Phi must then be set to an abritrary constant.*/
    // }

}

/*To call during the initialization of your solver, before the begin of the time loop*/
/*Maybe you should need to add an argument to specify the number of unknows*/
/*Modification to do in this function :*/
/*  OKKKKK -Specify the number of unknows*/
/*  OKKKKK -Specify the number of non-zero diagonals in the sparse matrix*/
PetscErrorCode initialize_poisson_solver(Poisson_data* data, int Nx, int Ny)
{
    PetscInt rowStart; /*rowStart = 0*/
    PetscInt rowEnd; /*rowEnd = the number of unknows*/
    PetscErrorCode ierr;

	int nphi = (Nx - 1) * (Ny - 1); /*WRITE HERE THE NUMBER OF UNKNOWS*/

    /* Create the right-hand-side vector : b */
    VecCreate(PETSC_COMM_WORLD, &(data->b));
    VecSetSizes(data->b, PETSC_DECIDE, nphi);
    VecSetType(data->b,VECSTANDARD);

    /* Create the solution vector : x */
    VecCreate(PETSC_COMM_WORLD, &(data->x));
    VecSetSizes(data->x, PETSC_DECIDE, nphi);
    VecSetType(data->x,VECSTANDARD);

    /* Create and assemble the Laplacian matrix : A  */
    MatCreate(PETSC_COMM_WORLD, &(data->A));
    MatSetSizes(data->A, PETSC_DECIDE, PETSC_DECIDE, nphi , nphi);
    MatSetType(data->A, MATAIJ);
    MatSeqAIJSetPreallocation(data->A,2 * (Ny - 1) + 1, NULL); // /*SET HERE THE NUMBER OF NON-ZERO DIAGONALS*/
    MatGetOwnershipRange(data->A, &rowStart, &rowEnd);

    computeLaplacianMatrix(data->A, rowStart, rowEnd, Nx, Ny);
    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    // MatView(data->A, PETSC_VIEWER_STDOUT_WORLD);

    /* Create the Krylov context */
    KSPCreate(PETSC_COMM_WORLD, &(data->sles));
    KSPSetOperators(data->sles, data->A, data->A);
    KSPSetType(data->sles,KSPGMRES); //KSPGMRES seems the best, it will not be used if PC LU.
    PC prec;
    KSPGetPC(data->sles, &prec);
    PCSetType(prec,PCLU);
    KSPSetFromOptions(data->sles); // to uncomment if we want to specify the solver to use in command line. Ex: mpirun -ksp_type gmres -pc_type gamg
    KSPSetTolerances(data->sles, 1.e-12, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetReusePreconditioner(data->sles,PETSC_TRUE);
    KSPSetUseFischerGuess(data->sles,1,4);
    KSPGMRESSetPreAllocateVectors(data->sles);

    PetscPrintf(PETSC_COMM_WORLD, "Assembly of Mattrix and Vectors is done \n");

    return ierr;
}

/*To call after the simulation to free the vectors needed for Poisson equation*/
/*Modification to do : nothing */
void free_poisson_solver(Poisson_data* data){
    MatDestroy(&(data->A));
    VecDestroy(&(data->b));
    VecDestroy(&(data->x));
    KSPDestroy(&(data->sles));
}
