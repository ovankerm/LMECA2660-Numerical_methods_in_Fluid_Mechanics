#include "src/poisson.h"
#include "src/functions.h"

int main(int argc, char *argv[]){

    // PetscInitialize(&argc, &argv, 0, 0);

    problem_struct *problem = create_problem(3);

    int i;
    for(i = 0; i < 100; i++){
        add_phi_P(problem);
    }

    print_mesh(problem->mesh);
    free_problem(problem);

    // PetscFinalize();
}
