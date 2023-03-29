#include "src/poisson.h"

int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);

    /*WRITE YOUR PROJECT ...*/
    Poisson_data *data = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(data);

    VecView(data->b, PETSC_VIEWER_STDOUT_WORLD);
    VecSet(data->b, 1.0);
    VecView(data->b, PETSC_VIEWER_STDOUT_WORLD);

    free_poisson_solver(data);
    free(data);

    PetscFinalize();

}
