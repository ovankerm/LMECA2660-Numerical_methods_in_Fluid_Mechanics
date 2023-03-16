#include "src/functions.h"

const int N = 32;

int main(int argc, char **argv){
    int i, j;
    problemStruct *pE2 = initProblem(N, I4, "output/E2", 0, 3.0/5.0, 0);
    computeDiagnostics(pE2);
    for(i = 0; i < 2; i++){
        for(j = 0; j < N; j++){
            RK4Iteration(pE2);
            computeDiagnostics(pE2);
        }
        problemToFile(pE2);
    }
    freeProblem(pE2);

    return EXIT_SUCCESS;
}