#include "functions.h"

int IND(i, N) {return i % N;}

const double RK4_alpha[4] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
const double RK4_gamma[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

problem *initProblem(int N, double L){
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
    fprintf(file, "L = %12.12f  N = %d  h = %12.12f\n", problem->L, problem->N, problem->h);
    for(i = 0; i < problem->N; i++){
        fprintf(file, "i = %d  Xi = %e Ui = %e\n", i, problem->X[i], problem->U[i]);
    }
    fclose(file);
}

void computeDU(double *dU, double dt, double *Ul, double tl, double N, double h, double c){
    int i;
    for(i = 0; i < N; i++){
        dU[IND(i, N)] = -c * dt/(2 * h) * (Ul[IND(i+1, N)] - Ul[IND(i-1, N)]);
    }
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
        computeDU(dU, problem->dt, Ul, tl, problem->N, problem->h, problem->c);
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
        problem->U[i] = exp(-problem->X[i] * problem->X[i]/(problem->sigma * problem->sigma));
    }

}