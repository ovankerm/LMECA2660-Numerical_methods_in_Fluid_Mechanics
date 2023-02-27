#include "src/functions.h"
#include "src/thomas.h"

int main(int argc, char **argv){
    problem *problem = initProblem(100, 10.5);
    problemToFile(problem);
    freeProblem(problem);

    return EXIT_SUCCESS;
}