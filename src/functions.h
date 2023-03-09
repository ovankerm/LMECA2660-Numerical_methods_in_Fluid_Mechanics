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
} myProblem;

// Space integrators
void E2(double *U, double *dU, double h, double c, double dt, int N);
void E4(double *U, double *dU, double h, double c, double dt, int N);
void I4(double *U, double *dU, double h, double c, double dt, int N);
void ED(double *U, double *dU, double h, double c, double dt, int N);

// Time integration
void RK4Iteration(myProblem *problem);

myProblem *initProblem(int N, double L, void (*integrator)(double *, double *, double , double , double , int));
void freeProblem(myProblem *problem);
void problemToFile(myProblem *problem, const char *basename);

void computeDiagnostics(myProblem *problem);

void initialConditionGaussian(myProblem *problem);
double exactGaussian(double x, double t, double sigma, double c);

#endif