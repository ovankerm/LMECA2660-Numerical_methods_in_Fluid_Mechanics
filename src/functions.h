#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef functions_h
#define functions_h

typedef struct{
    int N;
    double h;
    double sigma;
    double t;
    double dt;
    void (*integrator)(double *, double *, double , double , int);
    double *X;
    double *U;

    FILE *diagnostics_file;
} problemStruct;

// Space integrators
void E2(double *U, double *dU, double h, double dt, int N);
void E4(double *U, double *dU, double h, double dt, int N);
void I4(double *U, double *dU, double h, double dt, int N);
void ED(double *U, double *dU, double h, double dt, int N);

// Time integration
void RK4Iteration(problemStruct *problem);

problemStruct *initProblem(int N, void (*integrator)(double *, double *, double , double , int), const char *basename);
void freeProblem(problemStruct *problem);
void problemToFile(problemStruct *problem, const char *basename);

void computeDiagnostics(problemStruct *problem);

void initialConditionGaussian(problemStruct *problem);
double exactGaussian(double x, double t, double sigma);

#endif