#include "functions.h"

double u_x(int i, double h){
    return i*h;
}

double u_y(int j, double h){
    return -h/2 + j*h;
}

double v_x(int i, double h){
    return -h/2 + i*h;
}

double v_y(int j, double h){
    return j*h;
}

double T_x(int i, double h){
    return -h/2 + i*h;
}

double T_y(int j, double h){
    return -h/2 + j*h;
}

problem_struct *create_problem(int nx, int mixer){
    /*
    mixer = 1 if mixer or = 0 if no mixer
    nx : number of cells in the horizontal direction + 1
    */
    problem_struct *problem = malloc(sizeof(problem_struct));
    
    problem->sim_data = init_data(nx, mixer);

    problem->poiss_data = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(problem->poiss_data, problem->sim_data->nx, problem->sim_data->ny, problem->sim_data->h);

    problem->dt = fmin(0.5 * problem->sim_data->h * problem->sim_data->h * SQRT_GR, 0.5 * problem->sim_data->h/0.25);
    problem->t = 0.0;

    problem->iter = 0;

    problem->Reh = 0.0;
    problem->Rehw = 0.0;

    FILE* file = fopen("output/diagnostics.txt","w");
    problem->diag_file = file;

    return problem;
}

data_sim *init_data(int nx, int mixer){
    int i;
    data_sim *sim_data = malloc(sizeof(data_sim));

    sim_data->nx = nx;
    sim_data->h = 2.0/(3.0 * (nx-1));

    int ny = 1.0/sim_data->h + 1.0;
    sim_data->ny = ny;

    sim_data->nxu = nx;
    sim_data->nyu = ny+1;
    sim_data->u_data = (double *) calloc(4 * nx * (ny+1), sizeof(double));
    sim_data->u = (double **) malloc(nx * sizeof(double *));
    sim_data->u_star = (double **) malloc(nx * sizeof(double *));
    sim_data->H_u = (double **) malloc(nx * sizeof(double *));
    sim_data->H_1_u = (double **) malloc(nx * sizeof(double *));
    sim_data->nxu_low = (int) (2.0/(15.0 * sim_data->h));
    sim_data->nxu_high = (int) (8.0/(15.0 * sim_data->h)) + 1;
    sim_data->nyu_low = (int) (1.0/2.0 + 2.0/(15.0 * sim_data->h));
    sim_data->nyu_high = (int) (1.0/2.0 + 8.0/(15.0 * sim_data->h)) + 1;
    sim_data->n_u_points = 0;
    sim_data->u_i = (int *) calloc((sim_data->nyu_high - sim_data->nyu_low + 1) * (sim_data->nxu_high - sim_data->nxu_low + 1), sizeof(int));
    sim_data->u_j = (int *) calloc((sim_data->nyu_high - sim_data->nyu_low + 1) * (sim_data->nxu_high - sim_data->nxu_low + 1), sizeof(int));
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
    sim_data->nxv_low = (int) (1.0/2.0 + 2.0/(15.0 * sim_data->h));
    sim_data->nxv_high = (int) (1.0/2.0 + 8.0/(15.0 * sim_data->h)) + 1;
    sim_data->nyv_low = (int) (2.0/(15.0 * sim_data->h));
    sim_data->nyv_high = (int) (8.0/(15.0 * sim_data->h)) + 1;
    sim_data->n_v_points = 0;
    sim_data->v_i = (int *) calloc((sim_data->nyv_high - sim_data->nyv_low + 1) * (sim_data->nxv_high - sim_data->nxv_low + 1), sizeof(int));
    sim_data->v_j = (int *) calloc((sim_data->nyv_high - sim_data->nyv_low + 1) * (sim_data->nxv_high - sim_data->nxv_low + 1), sizeof(int));
    for(i = 0; i < nx+1; i++){
        sim_data->v[i] = sim_data->v_data + i * ny;
        sim_data->v_star[i] = sim_data->v_data + 1 * (nx+1) * ny + i * ny;
        sim_data->H_v[i] = sim_data->v_data + 2 * (nx+1) * ny + i * ny;
        sim_data->H_1_v[i] = sim_data->v_data + 3 * (nx+1) * ny + i * ny;
    }

    sim_data->nxT = nx+1;
    sim_data->nyT = ny+1;
    sim_data->T_data = (double *) calloc(4 * (nx+1) * (ny+1), sizeof(double));
    sim_data->T = (double **) malloc((nx+1) * sizeof(double *));
    sim_data->H_T = (double **) malloc((nx+1) * sizeof(double *));
    sim_data->H_1_T = (double **) malloc((nx+1) * sizeof(double *));
    sim_data->T_1 = (double **) malloc((nx+1) * sizeof(double *));
    sim_data->nxT_low = (int) (1.0/2.0 + 2.0/(15.0 * sim_data->h));
    sim_data->nxT_high = (int) (1.0/2.0 + 8.0/(15.0 * sim_data->h)) + 1;
    sim_data->nyT_low = (int) (1.0/2.0 + 2.0/(15.0 * sim_data->h));
    sim_data->nyT_high = (int) (1.0/2.0 + 8.0/(15.0 * sim_data->h)) + 1;
    sim_data->n_T_points = 0;
    sim_data->T_i = (int *) calloc((sim_data->nyT_high - sim_data->nyT_low + 1) * (sim_data->nxT_high - sim_data->nxT_low + 1), sizeof(int));
    sim_data->T_j = (int *) calloc((sim_data->nyT_high - sim_data->nyT_low + 1) * (sim_data->nxT_high - sim_data->nxT_low + 1), sizeof(int));
    for(i = 0; i < nx+1; i++){
        sim_data->T[i] = sim_data->T_data + i * (ny+1);
        sim_data->H_T[i] = sim_data->T_data + 1 * (nx+1) * (ny+1) + i * (ny+1);
        sim_data->H_1_T[i] = sim_data->T_data + 2 * (nx+1) * (ny+1) + i * (ny+1);
        sim_data->T_1[i] = sim_data->T_data + 3 * (nx+1) * (ny+1) + i * (ny+1);
    }

    sim_data->nxP = nx-1;
    sim_data->nyP = ny-1;
    sim_data->P_data = (double *) calloc((nx-1) * (ny-1), sizeof(double));
    sim_data->P = (double **) malloc((nx-1) * sizeof(double *));
    sim_data->phi_data = (double *) calloc((nx-1) * (ny-1), sizeof(double));
    sim_data->phi = (double **) malloc((nx-1) * sizeof(double *));
    for(i = 0; i < nx-1; i++){
        sim_data->P[i] = sim_data->P_data + i * (ny-1);
        sim_data->phi[i] = sim_data->phi_data + i * (ny-1);
    }

    if(mixer) get_indices(sim_data);

    return sim_data;
}

void free_problem(problem_struct *problem){
    free_data_sim(problem->sim_data);
    fclose(problem->diag_file);
    free(problem);
}

void free_data_sim(data_sim *sim_data){
    free(sim_data->u_data);
    free(sim_data->u);
    free(sim_data->u_star);
    free(sim_data->H_u);
    free(sim_data->H_1_u);
    free(sim_data->u_i);
    free(sim_data->u_j);

    free(sim_data->v_data);
    free(sim_data->v);
    free(sim_data->v_star);
    free(sim_data->H_v);
    free(sim_data->H_1_v);
    free(sim_data->v_i);
    free(sim_data->v_j);

    free(sim_data->P_data);
    free(sim_data->P);

    free(sim_data->phi_data);
    free(sim_data->phi);

    free(sim_data->T_data);
    free(sim_data->T);
    free(sim_data->T_1);
    free(sim_data->H_T);
    free(sim_data->H_1_T);
    free(sim_data->T_i);
    free(sim_data->T_j);

    free(sim_data);
}

void get_indices(data_sim *data){
    int i,j;
    
    for(i = data->nxu_low; i <= data->nxu_high; i++){
        for(j = data->nyu_low; j <= data->nyu_high; j++){
            if((u_x(i, data->h) - 1.0/3.0) * (u_x(i, data->h) - 1.0/3.0) + (u_y(j, data->h) - 1.0/3.0) * (u_y(j, data->h) - 1.0/3.0) <= 1.0/25.0){
                data->u_i[data->n_u_points] = i;
                data->u_j[data->n_u_points] = j;
                data->n_u_points++;
            }
        }
    }

    for(i = data->nxv_low; i <= data->nxv_high; i++){
        for(j = data->nyv_low; j <= data->nyv_high; j++){
            if((v_x(i, data->h) - 1.0/3.0) * (v_x(i, data->h) - 1.0/3.0) + (v_y(j, data->h) - 1.0/3.0) * (v_y(j, data->h) - 1.0/3.0) <= 1.0/25.0){
                data->v_i[data->n_v_points] = i;
                data->v_j[data->n_v_points] = j;
                data->n_v_points++;
            }
        }
    }

    for(i = data->nxT_low; i <= data->nxT_high; i++){
        for(j = data->nyT_low; j <= data->nyT_high; j++){
            if((T_x(i, data->h) - 1.0/3.0) * (T_x(i, data->h) - 1.0/3.0) + (T_y(j, data->h) - 1.0/3.0) * (T_y(j, data->h) - 1.0/3.0) <= 1.0/25.0){
                data->T_i[data->n_T_points] = i;
                data->T_j[data->n_T_points] = j;
                data->n_T_points++;
            }
        }
    }
}

void iterate(problem_struct *problem){
    problem->t += problem->dt;
    problem->iter++;
    data_sim *data = problem->sim_data;
    copy_H_n(data);
    compute_BC(data);
    compute_H_n(data);
    compute_v_star(data, problem->dt, problem->t);
    poisson_solver(problem->poiss_data, data->nx, data->ny, data->u_star, data->v_star, problem->dt, data->phi_data);
    compute_v(data, problem->dt);
    add_phi_P(data);
    compute_T(data, problem->dt, problem->t);
}

void problem_to_file(problem_struct *problem){
    printf("------ %d ------\n", problem->iter);
    T_to_file(problem);
    U_to_file(problem);
    omega_to_file(problem);
    v_to_file(problem);
}

void T_to_file(problem_struct *problem){
    int i, j;
    const char *basename = "output/T-%d.txt";
    char filename[256];
    sprintf(filename,basename,problem->iter);
    data_sim *data = problem->sim_data;
    FILE* file = fopen(filename,"w");
    fprintf(file, "Nx : %d Ny : %d h : %e t : %e\n", data->nx, data->ny, data->h, problem->t);
    for(i = 1; i < data->nxT-1; i++){
        for(j = 1; j < data->nyT-1; j++){
            fprintf(file, "%d %d %e\n", i, j, data->T[i][j]);
        }
    }
    fclose(file);
}

void v_to_file(problem_struct *problem){
    int i, j;
    const char *basename = "output/v-%d.txt";
    char filename[256];
    sprintf(filename,basename,problem->iter);
    data_sim *data = problem->sim_data;
    FILE* file = fopen(filename,"w");
    fprintf(file, "Nx : %d Ny : %d h : %e t : %e\n", data->nx, data->ny, data->h, problem->t);
    for(i = 0; i < data->nxv; i++){
        for(j = 0; j < data->nyv; j++){
            fprintf(file, "%d %d %e\n", i, j, data->v[i][j]);
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
    fprintf(file, "Nx : %d Ny : %d h : %e t : %e\n", data->nx, data->ny, data->h, problem->t);
    for(i = 0; i < data->nx; i++){
        for(j = 0; j < data->ny; j++){
            fprintf(file, "%d %d %e\n", i, j, sqrt((data->v[i+1][j] + data->v[i][j])/2 * (data->v[i+1][j] + data->v[i][j])/2
                                                                + (data->u[i][j+1] + data->u[i][j])/2 * (data->u[i][j+1] + data->u[i][j])/2));
        }
    }
    fclose(file);
}

void omega_to_file(problem_struct *problem){
    int i, j;
    const char *basename = "output/w-%d.txt";
    char filename[256];
    sprintf(filename,basename,problem->iter);
    data_sim *data = problem->sim_data;
    FILE* file = fopen(filename,"w");
    fprintf(file, "Nx : %d Ny : %d h : %e t : %e\n", data->nx, data->ny, data->h, problem->t);
    for(i = 0; i < data->nx; i++){
        for(j = 0; j < data->ny; j++){
            fprintf(file, "%d %d %e\n", i, j, 1/data->h * (data->v[i+1][j] - data->v[i][j] - data->u[i][j+1] + data->u[i][j]));
        }
    }
    fclose(file);
}

void compute_Rehw(problem_struct *problem){
    double Rehw_max = 0, Rehw;
    int i,j;
    data_sim *data = problem->sim_data;
    for(i = 0; i < data->nx; i++){
        for(j = 0; j < data->ny; j++){
            Rehw = data->h * fabs(data->v[i+1][j] - data->v[i][j] - data->u[i][j+1] + data->u[i][j]) * SQRT_GR;
            Rehw_max = fmax(Rehw, Rehw_max);
        }
    }
    problem->Rehw = Rehw_max;
}

void compute_Reh(problem_struct *problem){
    double Reh_max = 0, Reh;
    int i,j;
    data_sim *data = problem->sim_data;
    for(i = 0; i < data->nx; i++){
        for(j = 0; j < data->ny; j++){
            Reh = data->h * (fabs((data->v[i+1][j] + data->v[i][j])/2) + fabs((data->u[i][j+1] + data->u[i][j])/2)) * SQRT_GR;
            Reh_max = fmax(Reh, Reh_max);
        }
    }
    problem->Reh = Reh_max;
}


void compute_BC(data_sim *data){
    int i;
    double **u = data->u;
    double **T = data->T;
    double **v = data->v;

    // horizontal walls
    for(i = 1; i < data->nxu - 1; i++){
        u[i][0] = -0.2 * (u[i][3] - 5.0 * u[i][2] + 15.0 * u[i][1]);
        u[i][data->nyu-1] = u[i][data->nyu-2];
    }
    for(i = 1; i < data->nxT - 1; i++){
        T[i][0] = T[i][1] + data->h;
        T[i][data->nyT-1] = T[i][data->nyT-2] * (L0/data->h - 0.5)/(L0/data->h + 0.5);
    }

    // vertical walls
    for(i = 1; i < data->nyv - 1; i++){
        v[0][i] = -0.2 * (v[3][i] - 5.0 * v[2][i] + 15.0 * v[1][i]);
        v[data->nxv-1][i] = -0.2 * (v[data->nxv-4][i] - 5.0 * v[data->nxv-3][i] + 15.0 * v[data->nxv-2][i]);
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

void compute_T(data_sim *data, double dt, double t){
    int i, j;
    double trash;
    double **T = data->T;
    double **T_1 = data->T_1;
    for(i = 0; i < data->nxT; i++){
        for(j = 0; j < data->nyT; j++){
            T_1[i][j] = T[i][j];
        }
    }
    for(i = 1; i < data->nxT - 1; i++){
        for(j = 1; j < data->nyT - 1; j++){
            T[i][j] += dt * (-0.5 * (3.0 * data->H_T[i][j] - data->H_1_T[i][j])
                            + PR_M1/SQRT_GR * (T_1[i+1][j] + T_1[i-1][j] + T_1[i][j+1] + T_1[i][j-1] - 4.0 * T_1[i][j])/(data->h * data->h));
        }
    }

    double Ts = compute_Ts(data);

    for(i = 0; i < data->n_T_points; i++){
        if(in_mixer(T_x(data->T_i[i], data->h), T_y(data->T_j[i], data->h), t, &trash, &trash)){
            data->T[data->T_i[i]][data->T_j[i]] /= (1.0 + DT_DTAU);
            data->T[data->T_i[i]][data->T_j[i]] += Ts;
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

void compute_v_star(data_sim *data, double dt, double t){
    int i, j;
    double us, vs;
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
            data->u_star[i][j] = u[i][j] + dt * (-0.5 * (3.0 * H_u[i][j] - H_1_u[i][j])
                                            - 1.0/data->h * (P[i][j-1] - P[i-1][j-1])
                                            + 1.0/SQRT_GR * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4.0 * u[i][j])/(data->h * data->h));
        }
    }

    for(i = 0; i < data->n_u_points; i++){
        if(in_mixer(u_x(data->u_i[i], data->h), u_y(data->u_j[i], data->h), t, &us, &vs)){
            data->u_star[data->u_i[i]][data->u_j[i]] /= (1.0 + DT_DTAU);
            data->u_star[data->u_i[i]][data->u_j[i]] += us;
        }
    }

    for(i = 1; i < data->nxv - 1; i++){
        for(j = 1; j < data->nyv - 1; j++){
            data->v_star[i][j] = v[i][j] + dt * (-0.5 * (3.0 * H_v[i][j] - H_1_v[i][j])
                                            - 1.0/data->h * (P[i-1][j] - P[i-1][j-1])
                                            + 1.0/SQRT_GR * (v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1] - 4.0 * v[i][j])/(data->h * data->h)
                                            + 0.5 * (T[i][j] + T[i][j+1]));
        }
    }

    for(i = 0; i < data->n_v_points; i++){
        if(in_mixer(v_x(data->v_i[i], data->h), v_y(data->v_j[i], data->h), t, &us, &vs)){
            data->v_star[data->v_i[i]][data->v_j[i]] /= (1.0 + DT_DTAU);
            data->v_star[data->v_i[i]][data->v_j[i]] += vs;
        }
    }
}

void compute_H_n(data_sim *data){
    int i, j;
    double **u = data->u;
    double **v = data->v;
    double **T = data->T;
    double f = 1.0/(4.0 * data->h);
    for(i = 1; i < data->nxu - 1; i++){
        for(j = 1; j < data->nyu - 1; j++){
            data->H_u[i][j] = 1.0 * f * (
                                        + (u[i+1][j] + u[i][j]) * (u[i+1][j] + u[i][j])
                                        - (u[i][j] + u[i-1][j]) * (u[i][j] + u[i-1][j])
                                        + (u[i][j+1] + u[i][j]) * (v[i][j] + v[i+1][j])
                                        - (u[i][j-1] + u[i][j]) * (v[i][j-1] + v[i+1][j-1])
                                        );
        }
    }

    for(i = 1; i < data->nxv - 1; i++){
        for(j = 1; j < data->nyv - 1; j++){
            data->H_v[i][j] = 1.0 * f * (
                                        + (u[i][j] + u[i][j+1]) * (v[i+1][j] + v[i][j])
                                        - (u[i-1][j] + u[i-1][j+1]) * (v[i][j] + v[i-1][j])
                                        + (v[i][j+1] + v[i][j]) * (v[i][j] + v[i][j+1])
                                        - (v[i][j-1] + v[i][j]) * (v[i][j-1] + v[i][j])
                                        );
        }
    }

    for(i = 1; i < data->nxT - 1; i++){
        for(j = 1; j < data->nyT - 1; j++){
            data->H_T[i][j] = 2.0 * f * (
                                u[i][j] * (T[i+1][j] - T[i][j])
                                + u[i-1][j] * (T[i][j] - T[i-1][j])
                                + v[i][j] * (T[i][j+1] - T[i][j])
                                + v[i][j-1] * (T[i][j] - T[i][j-1])
                                );
        }
    }
}

void copy_H_n(data_sim *data){
    int i,j;
    for(i = 1; i < data->nxu - 1; i++){
        for(j = 1; j < data->nyu - 1; j++){
            data->H_1_u[i][j] = data->H_u[i][j];
        }
    }

    for(i = 1; i < data->nxv - 1; i++){
        for(j = 1; j < data->nyv - 1; j++){
            data->H_1_v[i][j] = data->H_v[i][j];
        }
    }

    for(i = 1; i < data->nxT - 1; i++){
        for(j = 1; j < data->nyT - 1; j++){
            data->H_1_T[i][j] = data->H_T[i][j];
        }
    }
}

int in_mixer(double x, double y, double t, double *us, double *vs){
    double X = x - 1.0/3.0;
    double Y = y - 1.0/3.0;
    double r = sqrt(X * X + Y * Y);
    double theta;
    if(fabs(X) < 1e-12 && fabs(Y) < 1e-12){
        *us = 0;
        *vs = 0;
        return 1;
    }
    if(fabs(X) < 1e-12) theta = Y > 0 ? M_PI/2 : -M_PI/2;
    else theta = X > 0 ? atan(Y/X) : atan(Y/X) + M_PI;

    *us = -r * OMEGA_S * sin(theta);
    *vs = r * OMEGA_S * cos(theta);
    
    return r <= 1.0/25.0 || r <= 1.0/5.0 * cos(3.0 * (theta - OMEGA_S * t));
}

double compute_Ts(data_sim *data){
    double Ts = 0.0;
    int i;
    if(data->n_T_points != 0){
        for(i = 0; i < data->n_T_points; i++){
            Ts += data->T[data->T_i[i]][data->T_j[i]];
        }
        Ts /= data->n_T_points;
    }
    return Ts;
}

void compute_diagnostics(problem_struct *problem){
    double Tavg = compute_Tavg(problem->sim_data);
    double sigT = compute_sigT(problem->sim_data, Tavg);
    double qe = compute_qe(problem->sim_data);
    double Ts = compute_Ts(problem->sim_data);
    compute_Reh(problem);
    compute_Rehw(problem);
    fprintf(problem->diag_file, "iter : %d Tavg : %e sigT : %e qe : %e Ts : %e Rehw : %e Reh : %e\n", problem->iter, Tavg, sigT, qe, Ts, problem->Rehw, problem->Reh);
}

double compute_Tavg(data_sim *data){
    double Tavg = 0.0;
    int i,j;
    for(i = 1; i < data->nxT - 1; i++){
        for(j = 1; j < data->nyT; j++){
            Tavg += data->T[i][j]/((data->nxT - 2) * (data->nyT - 2) - data->n_T_points);
        }
    }
    for(i = 0; i < data->n_T_points; i++){
        Tavg -= data->T[data->T_i[i]][data->T_j[i]]/((data->nxT - 2) * (data->nyT - 2) - data->n_T_points);
    }
    return Tavg;
}

double compute_sigT(data_sim *data, double Tavg){
    double sigT = 0.0;
    int i,j;
    for(i = 1; i < data->nxT - 1; i++){
        for(j = 1; j < data->nyT; j++){
            sigT += (data->T[i][j] - Tavg) * (data->T[i][j] - Tavg)/((data->nxT - 2) * (data->nyT - 2) - data->n_T_points);
        }
    }
    for(i = 0; i < data->n_T_points; i++){
        sigT -= (data->T[data->T_i[i]][data->T_j[i]] - Tavg) * (data->T[data->T_i[i]][data->T_j[i]] - Tavg)/((data->nxT - 2) * (data->nyT - 2) - data->n_T_points);
    }
    return sigT;
}

double compute_qe(data_sim *data){
    int i;
    double qe = 0;
    for(i = 1; i < data->nxT - 1; i++){
        qe += data->T[i][data->nyT - 1] - data->T[i][data->nyT - 2];
    }
    return qe * -3.0/2.0;
}