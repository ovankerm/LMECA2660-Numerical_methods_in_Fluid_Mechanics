#include "src/functions.h"

const int N = 128;

int main(int argc, char **argv){
    int i, j;
    problemStruct *pE2 = initProblem(N, E2, "output/E2_128points_wave_packet_mapping", 1, 3.0/5.0, 0);
    // computeDiagnostics(pE2);
    for(i = 0; i < 2; i++){
        for(j = 0; j < N/2; j++){
            RK4Iteration(pE2);
            // computeDiagnostics(pE2);
        }
        problemToFile(pE2);
    }
    freeProblem(pE2);

    problemStruct *pE4 = initProblem(N, E4, "output/E4_128points_wave_packet_mapping", 1, 3.0/5.0, 0);
    // computeDiagnostics(pE4);
    for(i = 0; i < 2; i++){
        for(j = 0; j < N/2; j++){
            RK4Iteration(pE4);
            // computeDiagnostics(pE4);
        }
        problemToFile(pE4);
    }
    freeProblem(pE4);

    problemStruct *pED = initProblem(N, ED, "output/ED_128points_wave_packet_mapping", 1, 3.0/5.0, 0);
    // computeDiagnostics(pED);
    for(i = 0; i < 2; i++){
        for(j = 0; j < N/2; j++){
            RK4Iteration(pED);
            // computeDiagnostics(pED);
        }
        problemToFile(pED);
    }
    freeProblem(pED);

    problemStruct *pI4 = initProblem(N, I4, "output/I4_128points_wave_packet_mapping", 1, 3.0/5.0, 0);
    // computeDiagnostics(pI4);
    for(i = 0; i < 2; i++){
        for(j = 0; j < N/2; j++){
            RK4Iteration(pI4);
            // computeDiagnostics(pI4);
        }
        problemToFile(pI4);
    }
    freeProblem(pI4);

    return EXIT_SUCCESS;
}