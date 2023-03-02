#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef functions_h
#define functions_h

typedef struct{
    int N;
    double h;
    double L;
    double sigma;
    double U_0;
    double c;
    double t;
    double dt;
    double *X;
    double *U;
} problem;


problem *initProblem(int N, double L);
void freeProblem(problem *problem);
void problemToFile(problem *problem);

void computeDU(double *dU, double dt, double *Ul, double tl, double N, double h, double c);

void RK4Iteration(problem *problem);

void initialConditionGaussian(problem *problem);


#endif