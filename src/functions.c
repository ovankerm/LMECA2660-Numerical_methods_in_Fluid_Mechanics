#include "functions.h"
#include "thomas.h"

int IND(i, N) {return (i+N) % N;}

const double RK4_alpha[4] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
const double RK4_gamma[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

void E2(double *U, double *dU, double h, double dt, int N){
    int i;
    for(i = 0; i < N; i++){
        dU[IND(i, N)] = -dt * 1/(2 * h) * (U[IND(i+1, N)] - U[IND(i-1, N)]);
    }
}

void E4(double *U, double *dU, double h, double dt, int N){
    int i;
    for(i = 0; i < N; i++){
        dU[IND(i, N)] = -dt * (4/(6 * h) * (U[IND(i+1, N)] - U[IND(i-1, N)]) - 1/(12 * h) * (U[IND(i+2, N)] - U[IND(i-2, N)]));
    }
}

void I4(double *U, double *dU, double h, double dt, int N){
    double *q = malloc(N * sizeof(double));
    int i;
    for(i = 0; i < N; i++){
        q[IND(i, N)] = -dt * 3/(4 * h) * (U[IND(i+1, N)] - U[IND(i-1, N)]);
    }
    solve_period_3diag(N, 1.0, 0.25, 0.25, dU, q);
    free(q);
}

void ED(double *U, double *dU, double h, double dt, int N){
    int i;
    for(i = 0; i < N; i++){
        dU[IND(i, N)] = -dt * 1/(6 * h) * (U[IND(i-2, N)] - 6 * U[IND(i-1, N)] + 3 * U[IND(i, N)] + 2 * U[IND(i+1, N)]);
    }
}

problemStruct *initProblem(int N, void (*integrator)(double *, double *, double , double , int), const char *basename){
    problemStruct *problem = malloc(sizeof(problemStruct));
    problem->N = N;
    problem->h = 1/((double) N);
    problem->sigma = 1/16.0;
    problem->t = 0.0;
    problem->dt = 0.5 * problem->h;
    problem->X = malloc(N * sizeof(double));
    problem->U = malloc(N * sizeof(double));
    problem->integrator = integrator;

    const char *filename_struct = "%s-diagnostics.txt";
    char filename[256];
    sprintf(filename,filename_struct,basename);
    FILE* file = fopen(filename,"w");
    problem->diagnostics_file = file;

    int i;
    for(i = 0; i < problem->N; i++){
        problem->X[i] = -1.0/2.0  + i * problem->h;
    }

    initialConditionGaussian(problem);

    return problem;
}

void freeProblem(problemStruct *problem){
    fclose(problem->diagnostics_file);
    free(problem->X);
    free(problem->U);
    free(problem);
}

void problemToFile(problemStruct *problem, const char *baseResult){
    int i;
    const char *basename = "%s-%.4f.txt";
    char filename[256];
    sprintf(filename,basename,baseResult,problem->t);
    FILE* file = fopen(filename,"w");

    fprintf(file, "sigma = %f  N = %d  h = %f\n", problem->sigma, problem->N, problem->h);
    for(i = 0; i < problem->N; i++){
        fprintf(file, "i = %d  Xi = %e Ui = %e\n", i, problem->X[i], problem->U[i]);
    }
    fclose(file);
}

void RK4Iteration(problemStruct *problem){
    int i,j;
    double ts, tl;
    double *Ul = malloc(problem->N * sizeof(double));
    double *Us = malloc(problem->N * sizeof(double));
    double *dU = calloc(problem->N, sizeof(double));

    for(i = 0; i < problem->N; i++){
        Us[i] = problem->U[i];
    }
    ts = problem->t;

    for(i = 0; i < 4; i++){
        for(j = 0; j < problem->N; j++){
            Ul[j] = Us[j] + RK4_alpha[i] * dU[j];
        }
        tl = ts + RK4_alpha[i] * problem->dt;
        
        problem->integrator(Ul, dU, problem->h, problem->dt, problem->N);

        for(j = 0; j < problem->N; j++){
            problem->U[j] += RK4_gamma[i] * dU[j];
        }
        problem->t += RK4_gamma[i] * problem->dt;
    }

    free(Ul);
    free(Us);
    free(dU);
}

void computeDiagnostics(problemStruct *problem){
    int i;
    double I = 0;
    double E = 0;
    double R = 0;
    double U, Ur;
    double factor = problem->h / problem->sigma;
    for(i = 0; i < problem->N; i++){
        U = problem->U[i];
        Ur = exactGaussian(problem->X[i], problem->t, problem->sigma);
        I += U;
        E += U*U;
        R += (U - Ur) * (U - Ur);
    }
    printf("factor : %f\n", factor);
    I *= factor;
    E *= factor;
    R *= factor;
    fprintf(problem->diagnostics_file, "t = %f  I = %e E = %e R = %e\n", problem->t, I, E, R);
}

void initialConditionGaussian(problemStruct *problem){
    int i;
    for(i = 0; i < problem->N; i++){
        double value = exp(-problem->X[i] * problem->X[i]/(problem->sigma * problem->sigma));
        problem->U[i] = fabs(value) < 1e-12 ? 0 : value;
    }
}

double exactGaussian(double x, double t, double sigma){
    return exp(-(x - t) * (x - t)/(sigma * sigma));
}