#include "functions.h"

problem_struct *create_problem(int nx){
    problem_struct *problem = malloc(sizeof(problem_struct));
    
    problem->sim_data = init_data(nx);


    problem->poiss_data = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(problem->poiss_data, problem->sim_data->nx, problem->sim_data->ny);

    problem->dt = 0.05 * problem->sim_data->h * problem->sim_data->h/SQRT_GR;
    problem->t = 0.0;

    problem->iter = 0;

    return problem;
}

data_sim *init_data(int nx){
    int i;
    data_sim *sim_data = malloc(sizeof(data_sim));

    sim_data->nx = nx;
    sim_data->h = 2.0/(3.0 * (nx-1));

    int ny = 1.0/sim_data->h + 1;
    sim_data->ny = ny;

    sim_data->nxu = nx;
    sim_data->nyu = ny+1;
    sim_data->u_data = (double *) calloc(4 * nx * (ny+1), sizeof(double));
    sim_data->u = (double **) malloc(nx * sizeof(double *));
    sim_data->u_star = (double **) malloc(nx * sizeof(double *));
    sim_data->H_u = (double **) malloc(nx * sizeof(double *));
    sim_data->H_1_u = (double **) malloc(nx * sizeof(double *));
    for(i = 0; i < nx; i++){
        sim_data->u[i] = sim_data->u_data + i * (ny+1);
        sim_data->u_star[i] = sim_data->u_data + 1 * nx * (ny+1) + i * (ny+1);
        sim_data->H_u[i] = sim_data->u_data + 2 * nx * (ny+1) + i * (ny+1);
        sim_data->H_1_u[i] = sim_data->u_data + 3 * nx * (ny+1) + i * (ny+1);
    }

    sim_data->nxv = nx+1;
    sim_data->nyv = ny;
    sim_data->v_data = (double *) calloc(4 * (nx+1) * ny, sizeof(double));
    sim_data->v = (double **) malloc((nx+1) * sizeof(double *));
    sim_data->v_star = (double **) malloc((nx+1) * sizeof(double *));
    sim_data->H_v = (double **) malloc((nx+1) * sizeof(double *));
    sim_data->H_1_v = (double **) malloc((nx+1) * sizeof(double *));
    for(i = 0; i < nx+1; i++){
        sim_data->v[i] = sim_data->v_data + i * ny;
        sim_data->v_star[i] = sim_data->v_data + 1 * (nx+1) * ny + i * ny;
        sim_data->H_v[i] = sim_data->v_data + 2 * (nx+1) * ny + i * ny;
        sim_data->H_1_v[i] = sim_data->v_data + 3 * (nx+1) * ny + i * ny;
    }

    sim_data->nxT = nx+1;
    sim_data->nyT = ny+1;
    sim_data->T_data = (double *) calloc(3 * (nx+1) * (ny+1), sizeof(double));
    sim_data->T = (double **) malloc((nx+1) * sizeof(double *));
    sim_data->H_T = (double **) malloc((nx+1) * sizeof(double *));
    sim_data->H_1_T = (double **) malloc((nx+1) * sizeof(double *));
    for(i = 0; i < nx+1; i++){
        sim_data->T[i] = sim_data->T_data + i * (ny+1);
        sim_data->H_T[i] = sim_data->T_data + 1 * (nx+1) * (ny+1) + i * (ny+1);
        sim_data->H_1_T[i] = sim_data->T_data + 2 * (nx+1) * (ny+1) + i * (ny+1);
    }

    sim_data->nxP = nx-1;
    sim_data->nyP = ny-1;
    sim_data->P_data = (double *) calloc(3 * (nx-1) * (ny-1), sizeof(double));
    sim_data->P = (double **) malloc((nx-1) * sizeof(double *));
    sim_data->phi_data = (double *) calloc(3 * (nx-1) * (ny-1), sizeof(double));
    sim_data->phi = (double **) malloc((nx-1) * sizeof(double *));
    for(i = 0; i < nx-1; i++){
        sim_data->P[i] = sim_data->P_data + i * (ny-1);
        sim_data->phi[i] = sim_data->phi_data + i * (ny-1);
    }

    return sim_data;
}

void free_problem(problem_struct *problem){
    free_data_sim(problem->sim_data);
    free(problem);
}

void free_data_sim(data_sim *sim_data){
    free(sim_data->u_data);
    free(sim_data->u);
    free(sim_data->u_star);
    free(sim_data->H_u);
    free(sim_data->H_1_u);

    free(sim_data->v_data);
    free(sim_data->v);
    free(sim_data->v_star);
    free(sim_data->H_v);
    free(sim_data->H_1_v);

    free(sim_data->P_data);
    free(sim_data->P);

    free(sim_data->phi_data);
    free(sim_data->phi);

    free(sim_data->T_data);
    free(sim_data->T);
    free(sim_data->H_T);
    free(sim_data->H_1_T);

    free(sim_data);
}

void iterate(problem_struct *problem){
    problem->t += problem->dt;
    problem->iter++;
    data_sim *data = problem->sim_data;
    copy_H_n(data);
    compute_BC(data);
    compute_H_n(data);
    compute_v_star(data, problem->dt);
    poisson_solver(problem->poiss_data, data->nx, data->ny, data->u_star, data->v_star, problem->dt, data->h, data->phi_data);
    compute_v(data, problem->dt);
    add_phi_P(data);
    compute_T(data, problem->dt);
}

// void print_mesh(problem_struct *problem){
//     int i,j;
//     printf("Nx : %d, Ny : %d, h : %f\n", problem->Nx, problem->Ny, problem->h);
//     printf("PRESSURE\n");
//     for(i = 1; i < problem->Nx; i++){
//         for(j = 1; j < problem->Ny; j++){
//             printf("T[%d, %d] = %f\n", i, j, problem->T[T_IND(i, j, problem->Ny)]);
//         }
//     }
//     printf("U\n");
//     for(i = 0; i < problem->Nx; i++){
//         for(j = 0; j < problem->Ny - 1; j++){
//             printf("U[%d, %d] = %f\n", i, j, problem->u[U_IND(i, j, problem->Ny)]);
//         }
//     }
// }

void problem_to_file(problem_struct *problem){
    T_to_file(problem);
    U_to_file(problem);
    V_to_file(problem);
}

void T_to_file(problem_struct *problem){
    int i, j;
    const char *basename = "output/T-%d.txt";
    char filename[256];
    sprintf(filename,basename,problem->iter);
    data_sim *data = problem->sim_data;
    FILE* file = fopen(filename,"w");
    fprintf(file, "Nx : %d Ny : %d h : %e\n", data->nx, data->ny, data->h);
    for(i = 1; i < data->nx; i++){
        for(j = 1; j < data->ny; j++){
            fprintf(file, "i : %d j : %d T : %e\n", i, j, data->T[i][j]);
        }
    }
    fclose(file);
}

void U_to_file(problem_struct *problem){
    int i, j;
    const char *basename = "output/U-%d.txt";
    char filename[256];
    sprintf(filename,basename,problem->iter);
    data_sim *data = problem->sim_data;
    FILE* file = fopen(filename,"w");
    fprintf(file, "Nx : %d Ny : %d h : %e\n", data->nx, data->ny, data->h);
    for(i = 0; i < data->nx; i++){
        for(j = 1; j < data->ny; j++){
            fprintf(file, "i : %d j : %d U : %e\n", i, j, data->u[i][j]);
        }
    }
    fclose(file);
}

void V_to_file(problem_struct *problem){
    int i, j;
    const char *basename = "output/V-%d.txt";
    char filename[256];
    sprintf(filename,basename,problem->iter);
    data_sim *data = problem->sim_data;
    FILE* file = fopen(filename,"w");
    fprintf(file, "Nx : %d Ny : %d h : %e\n", data->nx, data->ny, data->h);
    for(i = 1; i < data->nx; i++){
        for(j = 0; j < data->ny; j++){
            fprintf(file, "i : %d j : %d V : %e\n", i, j, data->v[i][j]);
        }
    }
    fclose(file);
}

// void U_star_to_file(problem_struct *problem){
//     int i, j;
//     const char *basename = "output/U_star-%d.txt";
//     char filename[256];
//     sprintf(filename,basename,problem->iter);
//     FILE* file = fopen(filename,"w");
//     fprintf(file, "Nx : %d Ny : %d h : %e\n", problem->Nx, problem->Ny, problem->h);
//     for(i = 0; i < problem->Nx; i++){
//         for(j = 1; j < problem->Ny; j++){
//             fprintf(file, "i : %d j : %d U_star : %e\n", i, j, problem->u_star[U_IND(i, j, problem->Ny)]);
//         }
//     }
//     fclose(file);
// }

// void V_star_to_file(problem_struct *problem){
//     int i, j;
//     const char *basename = "output/V_star-%d.txt";
//     char filename[256];
//     sprintf(filename,basename,problem->iter);
//     FILE* file = fopen(filename,"w");
//     fprintf(file, "Nx : %d Ny : %d h : %e\n", problem->Nx, problem->Ny, problem->h);
//     for(i = 1; i < problem->Nx; i++){
//         for(j = 0; j < problem->Ny; j++){
//             fprintf(file, "i : %d j : %d V_star : %e\n", i, j, problem->v_star[V_IND(i, j, problem->Ny)]);
//         }
//     }
//     fclose(file);
// }

// void P_to_file(problem_struct *problem){
//     int i, j;
//     const char *basename = "output/P-%.8f.txt";
//     char filename[256];
//     sprintf(filename,basename,problem->t);
//     FILE* file = fopen(filename,"w");
//     fprintf(file, "Nx : %d Ny : %d h : %e\n", problem->Nx, problem->Ny, problem->h);
//     for(i = 0; i < problem->Nx - 1; i++){
//         for(j = 0; j < problem->Ny - 1; j++){
//             fprintf(file, "i : %d j : %d P : %e\n", i, j, problem->P[P_IND(i, j, problem->Ny)]);
//         }
//     }
//     fclose(file);
// }

// void print_phi(double *phi, int Nx, int Ny){
//     int i, j;
//     printf("PHI\n");
//     for(j = 0; j < Ny-1; j++){
//         for(i = 0; i < Nx-1; i++){
//             printf("%f\t", phi[P_IND(i, j, Ny)]);
//         }
//         printf("\n");
//     }
// }

// void U_star(double *u_star, int Nx, int Ny){
//     int i, j;
//     FILE* file = fopen("output/U_star.txt","w");
//     for(i = 0; i < Nx; i++){
//         for(j = 1; j < Ny; j++){
//             fprintf(file, "i : %d j : %d U_star : %e\n", i, j, u_star[U_IND(i, j, Ny)]);
//         }
//     }
//     fclose(file);
// }

// void V_star(double *v_star, int Nx, int Ny){
//     int i, j;
//     FILE* file = fopen("output/V_star.txt","w");
//     for(i = 1; i < Nx; i++){
//         for(j = 0; j < Ny; j++){
//             fprintf(file, "i : %d j : %d V_star : %e\n", i, j, v_star[V_IND(i, j, Ny)]);
//         }
//     }
//     fclose(file);
// }


// void test_poisson(){
//     int Nx = 8, Ny = 8, i, j;
//     double dt = 0.1;
//     double h = 0.1;
//     Poisson_data *data = malloc(sizeof(Poisson_data));
//     initialize_poisson_solver(data, Nx, Ny);
//     double *phi = calloc((Nx-1)*(Ny-1), sizeof(double));
//     double *u_star = calloc(Nx * (Ny+1), sizeof(double));
//     double *v_star = calloc(Nx * (Ny+1), sizeof(double));
//     poisson_solver(data, phi, u_star, v_star, 0.5, 0.1, Nx, Ny);

//     for(i = 1; i < Nx - 1; i++){
//         for(j = 1; j < Ny; j++){
//             u_star[U_IND(i, j, Ny)] += - dt/h * (phi[P_IND(i, j - 1, Ny)] - phi[P_IND(i - 1, j - 1, Ny)]);
//         }
//     }

//     for(i = 1; i < Nx; i++){
//         for(j = 1; j < Ny - 1; j++){
//             // printf("%f\n", dt/h * (phi[P_IND(i - 1, j, Ny)] - phi[P_IND(i - 1, j - 1, Ny)]));
//             v_star[V_IND(i, j, Ny)] += - dt/h * (phi[P_IND(i - 1, j, Ny)] - phi[P_IND(i - 1, j - 1, Ny)]);
//         }
//     }

//     V_star(v_star, Nx, Ny);
//     U_star(u_star, Nx, Ny);

//     print_phi(phi, Nx, Ny);
//     free(phi);
//     free(u_star);
//     free(v_star);
//     free(data);
// }


void compute_BC(data_sim *data){
    int i;
    double **u = data->u;
    double **T = data->T;
    double **v = data->v;

    // horizontal walls
    for(i = 1; i < data->nxu - 1; i++){
        u[i][0] = -1.0/5.0 * (u[i][3] - 5.0 * u[i][2] + 15.0 * u[i][1]);
        u[i][data->nyu-1] = u[i][data->nyu-2];
    }
    for(i = 1; i < data->nxT - 1; i++){
        T[i][0] = T[i][1] + data->h;
        T[i][data->nyT-1] = T[i][data->nyT-2] * (L0/data->h - 0.5)/(L0/data->h + 0.5);
    }

    // vertical walls
    for(i = 1; i < data->nyv - 1; i++){
        v[0][i] = -1.0/5.0 * (v[3][i] - 5.0 * v[2][i] + 15.0 * v[1][i]);
        v[data->nxv-1][i] = -1.0/5.0 * (v[data->nxv-4][i] - 5.0 * v[data->nxv-3][i] + 15.0 * v[data->nxv-2][i]);
    }
    for(i = 1; i < data->nyT - 1; i++){
        T[0][i] = T[1][i];
        T[data->nxT-1][i] = T[data->nxT-2][i];
    }
}

void add_phi_P(data_sim *data){
    int i,j;
    for(i = 0; i < data->nxP; i++){
        for(j = 0; j < data->nyP; j++){
            data->P[i][j] += data->phi[i][j];
        }
    }
}

void compute_T(data_sim *data, double dt){
    int i, j;
    double **T = data->T;
    for(i = 1; i < data->nxT - 1; i++){
        for(j = 1; j < data->nyT - 1; j++){
            T[i][j] += dt * (-0.5 * (3 * data->H_T[i][j] - data->H_1_T[i][j])
                            + SQRT_GR * PR_M1 * (T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4.0 * T[i][j])/(data->h * data->h));
        }
    }
}

void compute_v(data_sim *data, double dt){
    int i, j;
    double **u = data->u;
    double **v = data->v;
    double **u_star = data->u_star;
    double **v_star = data->v_star;
    double **phi = data->phi;
    for(i = 1; i < data->nxu - 1; i++){
        for(j = 1; j < data->nyu - 1; j++){
            u[i][j] = u_star[i][j] - dt/data->h * (phi[i][j-1] - phi[i-1][j-1]);
        }
    }

    for(i = 1; i < data->nxv - 1; i++){
        for(j = 1; j < data->nyv - 1; j++){
            v[i][j] = v_star[i][j] - dt/data->h * (phi[i-1][j] - phi[i-1][j-1]);
        }
    }
}

void compute_v_star(data_sim *data, double dt){
    int i, j;
    double **u = data->u;
    double **v = data->v;
    double **P = data->P;
    double **T = data->T;
    double **H_u = data->H_u;
    double **H_v = data->H_v;
    double **H_1_u = data->H_1_u;
    double **H_1_v = data->H_1_v;
    for(i = 1; i < data->nxu - 1; i++){
        for(j = 1; j < data->nyu - 1; j++){
            data->u_star[i][j] = u[i][j] + dt * (-0.5 * (3 * H_u[i][j] - H_1_u[i][j])
                                            - 1.0/data->h * (P[i][j-1] - P[i-1][j-1])
                                            + SQRT_GR * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4.0 * u[i][j])/(data->h * data->h));
        }
    }

    for(i = 1; i < data->nxv - 1; i++){
        for(j = 1; j < data->nyv - 1; j++){
            data->v_star[i][j] = v[i][j] + dt * (-0.5 * (3 * H_v[i][j] - H_1_v[i][j])
                                            - 1.0/data->h * (P[i-1][j] - P[i-1][j-1])
                                            + SQRT_GR * (v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1] - 4.0 * v[i][j])/(data->h * data->h)
                                            + 0.5 * (T[i][j] + T[i][j+1]));
        }
    }
}

void compute_H_n(data_sim *data){
    int i, j;
    double **u = data->u;
    double **v = data->v;
    double **T = data->T;
    for(i = 1; i < data->nxu - 1; i++){
        for(j = 1; j < data->nyu - 1; j++){
            data->H_u[i][j] = 1/(4 * data->h) * ((u[i][j] + u[i+1][j]) * (u[i+1][j] - u[i][j])
                                                + (u[i][j] + u[i-1][j]) * (u[i][j] - u[i-1][j])
                                                + (v[i][j] + v[i+1][j]) * (u[i][j+1] - u[i][j])
                                                + (v[i][j-1] + v[i+1][j-1]) * (u[i][j] - u[i][j-1]));
        }
    }

    for(i = 1; i < data->nxv - 1; i++){
        for(j = 1; j < data->nyv - 1; j++){
            data->H_v[i][j] = 1/(4 * data->h) * ((u[i-1][j] + u[i-1][j+1]) * (v[i][j] - v[i-1][j])
                                                + (u[i][j] + u[i][j+1]) * (v[i+1][j] - v[i][j])
                                                + (v[i][j+1] + v[i][j]) * (v[i][j+1] - v[i][j])
                                                + (v[i][j] + v[i][j-1]) * (v[i][j] - v[i-1][j]));
        }
    }

    for(i = 1; i < data->nxT - 1; i++){
        for(j = 1; j < data->nyT; j++){
            data->H_T[i][j] = 1/(4 * data->h) * (u[i][j] * (T[i+1][j] - T[i][j]) + u[i-1][j] * (T[i][j] - T[i-1][j])
                                                + v[i][j] * (T[i][j+1] - T[i][j]) + v[i][j-1] * (T[i][j] - T[i][j-1]));
        }
    }
}

void copy_H_n(data_sim *data){
    memcpy(data->H_1_u[0], data->H_u[0], data->nxu * data->nyu * sizeof(double));
    memcpy(data->H_1_v[0], data->H_v[0], data->nxv * data->nyv * sizeof(double));
    memcpy(data->H_1_T[0], data->H_T[0], data->nxT * data->nyT * sizeof(double));
}
