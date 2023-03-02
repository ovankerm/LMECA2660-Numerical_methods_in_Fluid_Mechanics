#include "src/functions.h"
#include "src/thomas.h"

int main(int argc, char **argv){
    problem *problem = initProblem(100, 10.5);
    int i;
    for(i = 0; i < 1000; i++){
        RK4Iteration(problem);
    }
    problemToFile(problem);
    freeProblem(problem);
    return EXIT_SUCCESS;
}