#include "src/functions.h"

int main(int argc, char **argv){
    int i;

    myProblem *pE2 = initProblem(100, 16.0, ED);
    for(i = 0; i < 100000; i++){
        if((i%100) == 0){
            problemToFile(pE2, "output/E2");
        }
        RK4Iteration(pE2);
    }
    problemToFile(pE2, "output/E2");
    freeProblem(pE2);

    // myProblem *pE4 = initProblem(100, 16.0, E4);
    // for(i = 0; i < 10000; i++){
    //     if((i%100) == 0){
    //         problemToFile(pE4, "output/E4");
    //     }
    //     RK4Iteration(pE4);
    // }
    // problemToFile(pE4, "output/E4");
    // freeProblem(pE4);

    // myProblem *pI4 = initProblem(100, 16.0, I4);
    // for(i = 0; i < 10000; i++){
    //     if((i%100) == 0){
    //         problemToFile(pI4, "output/I4");
    //     }
    //     RK4Iteration(pI4);
    // }
    // problemToFile(pI4, "output/I4");
    // freeProblem(pI4);

    // myProblem *pED = initProblem(100, 16.0, ED);
    // for(i = 0; i < 10000; i++){
    //     if((i%100) == 0){
    //         problemToFile(pED, "output/ED");
    //     }
    //     RK4Iteration(pED);
    // }
    // problemToFile(pED, "output/ED");
    // freeProblem(pED);

    return EXIT_SUCCESS;
}