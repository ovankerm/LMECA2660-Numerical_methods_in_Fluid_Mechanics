#include "src/functions.h"

int main(int argc, char *argv[]){
    PetscInitialize(&argc, &argv, 0, 0);
    
    problem_struct *problem = create_problem(351, 1);

    int i, j;
    int N_frame = 1/(problem->dt);

    for(j = 0; j < 1000; j++){
        for(i = 0; i < N_frame; i++){
            iterate(problem);
            iterate(problem);
            iterate(problem);
            iterate(problem);
            iterate(problem);
        }
        problem_to_file(problem);
        compute_diagnostics(problem);
        if(problem->Rehw > 40){
            printf("\033[1;31mRehw too big !!\033[0m\n");
            return 0;
        }

        if(problem->Reh > 25){
            printf("\033[1;31mReh too big !!\033[0m\n");
            return 0;
        }
    }

    free_problem(problem);

    PetscFinalize();
    return 0;
}
