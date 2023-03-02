#include "src/functions.h"

int main(int argc, char **argv){
    int i;
    problem *pE2 = initProblem(100, 16.0, E2);
    for(i = 0; i < 1000; i++){
        RK4Iteration(pE2);
    }
    problemToFile(pE2, "output/E2.txt");
    freeProblem(pE2);

    problem *pE4 = initProblem(100, 16.0, E4);
    for(i = 0; i < 1000; i++){
        RK4Iteration(pE4);
    }
    problemToFile(pE4, "output/E4.txt");
    freeProblem(pE4);

    problem *pI4 = initProblem(100, 16.0, I4);
    for(i = 0; i < 1000; i++){
        RK4Iteration(pI4);
    }
    problemToFile(pI4, "output/I4.txt");
    freeProblem(pI4);
    return EXIT_SUCCESS;
}