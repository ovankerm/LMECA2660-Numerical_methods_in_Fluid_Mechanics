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
    double c;
    double t;
    double dt;
    void (*integrator)(double *, double *, double , double , double , int);
    double *X;
    double *U;

    double I;
    double E;
    double R;
} problem;

// Space integrators
void E2(double *U, double *dU, double h, double c, double dt, int N);
void E4(double *U, double *dU, double h, double c, double dt, int N);
void I4(double *U, double *dU, double h, double c, double dt, int N);

problem *initProblem(int N, double L, void (*integrator)(double *, double *, double , double , double , int));
void freeProblem(problem *problem);
void problemToFile(problem *problem, const char* filename);

void RK4Iteration(problem *problem);

void computeDiagnostics(problem *problem);

void initialConditionGaussian(problem *problem);
double exactGaussian(double x, double t, double sigma, double c);

#endif