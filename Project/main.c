#include "src/poisson.h"
#include "src/functions.h"

int main(int argc, char *argv[]){

    // PetscInitialize(&argc, &argv, 0, 0);

    problem_struct *problem = create_problem(5);

    int i;
    // for(i = 0; i < 4; i++){
    //     add_phi_P(problem);
    // }
    compute_BC(problem);
    compute_H_n(problem);
    // compute_v_star(problem);
    // compute_v(problem);
    compute_T(problem);

    for(i = 0; i < (problem->Nx + 1) * (problem->Ny); i++){
        problem->v[i] = problem->v_star[i];
    }

    print_mesh(problem);

    free_problem(problem);

    // PetscFinalize();
}
