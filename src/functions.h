#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef functions_h
#define functions_h

typedef struct{
    int N;
    double h;
    double L;
    double sigma;
    double U_0;
    double *X;
    double *U;
} problem;


problem *initProblem(int N, double L);
void freeProblem(problem *problem);
void problemToFile(problem *problem);

void initialConditionGaussian(problem *problem);


#endif