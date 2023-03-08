#include "src/functions.h"

int main(int argc, char **argv){
    int i;
    // problem *pE2 = initProblem(100, 16.0, E2);
    // for(i = 0; i < 1000; i++){
    //     RK4Iteration(pE2);
    // }
    // problemToFile(pE2, "output/E2.txt");
    // freeProblem(pE2);

    printf("NEW PROBLEM\n");

    myProblem *pE4 = initProblem(100, 16.0, ED);
    printf("after init : %f\n", pE4->U[0]);
    for(i = 0; i < 10000; i++){
        RK4Iteration(pE4);
        // printf("%d : %f\n", i, pE4->U[0]);
    }
    printf("after loop : %f\n", pE4->U[0]);
    problemToFile(pE4, "output/E4.txt");
    printf("after file : %f\n", pE4->U[0]);
    freeProblem(pE4);

    // problem *pI4 = initProblem(100, 16.0, I4);
    // for(i = 0; i < 1000; i++){
    //     RK4Iteration(pI4);
    // }
    // problemToFile(pI4, "output/I4.txt");
    // freeProblem(pI4);
    return EXIT_SUCCESS;
}