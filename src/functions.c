#include "functions.h"

problem *initProblem(int N, double L){
    problem *problem = malloc(sizeof(problem));
    problem->N = N;
    problem->L = L;
    problem->h = L/((double) N);
    problem->sigma = L/16.0;
    problem->X = malloc(N * sizeof(double));
    problem->U = malloc(N * sizeof(double));

    int i;
    for(i = 0; i < problem->N; i++){
        problem->X[i] = -problem->L/2  + i * problem->h;
    }

    initialConditionGaussian(problem);

    return problem;
}

void freeProblem(problem *problem){
    free(problem->X);
    free(problem->U);
    free(problem);
}

void problemToFile(problem *problem){
    int i;
    FILE *file = fopen("output/test.txt", "w");
    if(file == NULL) return;
    fprintf(file, "L = %lf  N = %d  h = %lf\n", problem->L, problem->N, problem->h);
    for(i = 0; i < problem->N; i++){
        fprintf(file, "i = %d  Xi = %lf  Ui = %lf\n", i, problem->X[i], problem->U[i]);
    }
    fclose(file);
}

void initialConditionGaussian(problem *problem){
    int i;
    for(i = 0; i < problem->N; i++){
        problem->U[i] = exp(-problem->X[i] * problem->X[i]/(problem->sigma * problem->sigma));
    }

}