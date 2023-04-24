#include "src/functions.h"

int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);

    // Nx == 2^n + 1
    problem_struct *problem = create_problem(65, 1);

    int i, j;

    for(j = 0; j < 500; j++){
        for(i = 0; i < 800; i++){
            iterate(problem);
        }
        problem_to_file(problem);
    }

    free_problem(problem);

    PetscFinalize();
}
