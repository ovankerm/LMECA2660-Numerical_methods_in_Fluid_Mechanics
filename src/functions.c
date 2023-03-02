#include "functions.h"
#include "thomas.h"

int IND(i, N) {return i % N;}

const double RK4_alpha[4] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
const double RK4_gamma[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

void E2(double *U, double *dU, double h, double c, double dt, int N){
    int i;
    for(i = 0; i < N; i++){
        dU[IND(i, N)] = -c * dt * 1/(2 * h) * (U[IND(i+1, N)] - U[IND(i-1, N)]);
    }
}

void E4(double *U, double *dU, double h, double c, double dt, int N){
    int i;
    for(i = 0; i < N; i++){
        dU[IND(i, N)] = -c * dt * (4/(6 * h) * (U[IND(i+1, N)] - U[IND(i-1, N)]) - 1/(12 * h) * (U[IND(i+2, N)] - U[IND(i-2, N)]));
    }
}

void I4(double *U, double *dU, double h, double c, double dt, int N){
    double *q = malloc(N * sizeof(double));
    int i;
    for(i = 0; i < N; i++){
        q[IND(i, N)] = -c * dt * 3/(4 * h) * (U[IND(i+1, N)] - U[IND(i-1, N)]);
    }
    solve_period_3diag(N, 1.0, 0.25, 0.25, dU, q);
    free(q);
}

problem *initProblem(int N, double L, void (*integrator)(double *, double *, double , double , double , int)){
    problem *problem = malloc(sizeof(problem));
    problem->N = N;
    problem->L = L;
    problem->h = L/((double) N);
    problem->sigma = L/16.0;
    problem->c = 1;
    problem->t = 0.0;
    problem->dt = 0.001;
    problem->X = malloc(N * sizeof(double));
    problem->U = malloc(N * sizeof(double));

    problem->integrator = integrator;

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

void problemToFile(problem *problem, const char *filename){
    int i;
    FILE *file = fopen(filename, "w");
    if(file == NULL) return;
    fprintf(file, "L = %f  N = %d  h = %f\n", problem->L, problem->N, problem->h);
    for(i = 0; i < problem->N; i++){
        fprintf(file, "i = %d  Xi = %e Ui = %e\n", i, problem->X[i], problem->U[i]);
    }
    fclose(file);
}

void RK4Iteration(problem *problem){
    int i,j;
    double ts, tl;
    double *Ul = malloc(problem->N * sizeof(double));
    double *Us = malloc(problem->N * sizeof(double));
    double *dU = malloc(problem->N * sizeof(double));
    memcpy(Us, problem->U, problem->N * sizeof(double));
    ts = problem->t;
    for(i = 0; i < 4; i++){
        for(j = 0; j < problem->N; j++){
            Ul[j] = Us[j] + RK4_alpha[i] * dU[j];
        }
        tl = ts + RK4_alpha[i] * problem->dt;
        
        problem->integrator(Ul, dU, problem->h, problem->c, problem->dt, problem->N);

        for(j = 0; j < problem->N; j++){
            problem->U[j] += RK4_gamma[i] * dU[j];
        }
        problem->t += RK4_gamma[i] * problem->dt;
    }

    free(Ul);
    free(Us);
    free(dU);
}

void initialConditionGaussian(problem *problem){
    int i;
    for(i = 0; i < problem->N; i++){
        double value = exp(-problem->X[i] * problem->X[i]/(problem->sigma * problem->sigma));
        problem->U[i] = fabs(value) < 1e-12 ? 0 : value;
    }
}