#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include "fmm.h"


void AllocateData(Data *data, int N, int level) {
    data->x = NULL;
    data->y = NULL;
    data->q = NULL;
    data->u = NULL;
    data->x = (double *) malloc(sizeof(double)*N);
    data->y = (double *) malloc(sizeof(double)*N);
    data->u = (double *) malloc(sizeof(double)*N);
    data->q = (double *) malloc(sizeof(double)*N);
    data->icell = (int *) malloc(sizeof(double)*N);
    data->gridsize = pow(2, level);
    data->Ncells = data->gridsize*data->gridsize;
    data->N = N;
    data->level = level;
}

void InitTestData(Data *data) {
    int iX[2];
    srand(0);
    for(int i = 0; i < data->N; i++) {
        data->x[i] = drand48();
        data->y[i] = drand48();
        iX[0] = data->x[i] * data->gridsize;
        iX[1] = data->y[i] * data->gridsize;
        data->icell[i] = getIndex(iX, data->level);
        data->u[i] = 0;
        data->q[i] = 1;
    }
}


void PrintData(Data *data) {
    int i;
    for(i = 0; i < data->N; i++) {
        fprintf(stdout, "%d %lf %lf %lf %lf\n", data->icell[i], data->x[i], data->y[i], data->q[i], data->u[i]);
        fflush(stdout);
    }
    fflush(stdout);
}

void DeallocateData(Data *data) {
    free(data->q);
    free(data->u);
    free(data->x);
    free(data->y);
    free(data->M);
    free(data->L);
    free(data->cells); 
    free(data->levelOffset);
    free(data->offset);
    free(data->icell);
}

void ParticleToMultipole(Data *data) {
    // Particle-to-Multipole code
    // e.g. getting the multipole of the lowest level cell
    // due to the particles located within that cell.
    // Convert particle postions into the
    // multipole expansion of the lowest level cell
    // Here we can just multiple by 4 store the result
    // in an integer variable because this allows us to
    // calculate the index of the cell on level 2.
    // e.g. if x = 0.2, ix = 0
    //         x = 0.3, ix = 1
    //         x = 0.5, ix = 2
    //         x = 0.8, ix = 3
    // The zeroth order multipole expansion is just
    // the "charge" of the source particle.
    int i, j;
    for (i = 0; i < data->Nleaf; i++) {
        for(j = data->offset[i]; j < data->offset[i+1]; j++) {
            data->M[data->cells[data->level*data->Ncells + i]+data->levelOffset[data->level]] += data->q[j];
        }
    }
}

void MultipoleToMultipole(Data *data) {
    int i, j, l;
    for(l = data->level; l > 2; l--) {
        //nx = pow(2, data->level);
        //Ncells = nx * nx;
        for(i = 0; i < data->Nleaf; i++) {
            j = data->cells[l*data->Ncells + i];
            data->M[j/4 + data->levelOffset[l-1]] += data->M[j + data->levelOffset[l]];
        }
    }
    // for(int i = 0; i < data->Ncells+data->levelOffset[data->level]; i++) {
    //     printf("%d %lf\n", i, data->M[i]);
    // }
}

void MultipoleToLocal(Data *data) {
    // Multipole-to-local code
    // These are calculated from the multipoles of each
    // "far" cell. In this code we classify near as
    // by being in a neighbouring cell (and hence the 
    // non neighbours are far). This has its
    // advantages and disadvantages. More sophisticated
    // is to use an adaptive tree structure and a
    // multipole acceptance criteria.
    
    // Loop over source boxes
    int iX[2], jX[2], i, j, l, nx, Ncells;
    double dx, dy, r;
    for(l = 2; l <= data->level; l++) {
        nx = pow(2, l);
        Ncells = nx * nx;
        for(i = 0; i < Ncells; i++) {
                // Loop over the target boxes
            getIX(iX, i);
            for(j = 0; j < Ncells; j++) {
                getIX(jX, j);
                        // If the box is "far"
                        // add the contribution to the local expansion.
                if (abs(iX[0]/2 - jX[0]/2) <= 1 && abs(iX[1]/2 - jX[1]/2) <=1) {
                    if (abs(iX[0] - jX[0]) > 1 || abs(iX[1] - jX[1]) > 1) {
                            // Can do this because particles uniformly
                            // distributed between 0 and 1, only a level 2
                            // tree.
                        dx = (double) (iX[0] - jX[0]) / (double) nx; // Have to cast otherwise int/int.
                        dy = (double) (iX[1] - jX[1]) / (double) nx;
                        r = sqrt(dx*dx + dy*dy);
                        i = getIndex(iX, data->level);
                        j = getIndex(jX, data->level);
                        data->L[i + data->levelOffset[l]] += data->M[j + data->levelOffset[l]] / r;
                    }
                }   
            }
        }
    }
}

void LocalToLocal(Data *data) {
    int l, i, j;
    for(l = 3; l <= data->level; l++) {
        // nx = pow(2, data->level);
        // Ncells = nx * nx;
        for(i = 0; i < data->Nleaf; i++) {
            j = data->cells[l*data->Ncells + i];
            data->L[j + data->levelOffset[l]] += data->L[j/4 + data->levelOffset[l-1]];
        }
    }
}

void LocalToParticle(Data *data) {
    // Local-to-particle code
    // Calculate the effect of the local expansion of a cell
    // on the particles within that cell.
    // Loop over each particle
    int i, j;
    for(i = 0; i < data->Nleaf; i++) {
        for(j = data->offset[i]; j < data->offset[i+1]; j++) {
            // Add the local expansion of that cell to the
            // potential at that particle.
            data->u[j] += data->L[data->cells[data->level*data->Ncells + i]+data->levelOffset[data->level]];
        }
    }
}


void ParticleToParticle(Data *data) {
// Particle-to-Particle code
    // Now need to add contributions directly from particles stored
    // in the nearest neighbour cells.
    int ic, jc; // indices for cells.
    double dx, dy, r;
    int i, iX[2], j, jX[2];
    for(ic = 0; ic < data->Nleaf; ic++) {
        getIX(iX, data->cells[data->level*data->Ncells + ic]);
        for(jc = 0; jc < data->Nleaf; jc++) {
            getIX(jX, data->cells[data->level*data->Ncells + jc]);
            if (abs(iX[0] - jX[0]) <= 1 && abs(iX[1] - jX[1]) <= 1) {
                for(i=data->offset[ic]; i < data->offset[ic+1]; i++) {
                    for(j=data->offset[jc]; j < data->offset[jc+1]; j++) {
                        dx = data->x[i] - data->x[j];
                        dy = data->y[i] - data->y[j];
                        r = sqrt(dx*dx + dy*dy);
                        // Check for r = 0.
                        if (r != 0) {
                            data->u[i] += data->q[j] / r;
                        }
                    }
                }
            }

        }
    }
}
    


void FMMComputePotential(Data *data) {
    double time[9];
    time[0] = omp_get_wtime();
    ClearPotential(data);
    time[1] = omp_get_wtime();
    ParticleSortTreeSetup(data);
    time[2] = omp_get_wtime();
    ParticleToMultipole(data);
    time[3] = omp_get_wtime();
    MultipoleToMultipole(data);
    time[4] = omp_get_wtime();
    MultipoleToLocal(data);
    time[5] = omp_get_wtime();
    LocalToLocal(data);
    time[6] = omp_get_wtime();
    LocalToParticle(data);
    time[7] = omp_get_wtime();
    ParticleToParticle(data);
    time[8] = omp_get_wtime();
    printf("TreeCreation\t%.15lf\nP2M\t%.15lf\nM2M\t%.15lf\nM2L\t%.15lf\nL2L\t%.15lf\nL2P\t%.15lf\nP2P\t%.15lf\n", time[2]-time[1], time[3]-time[2], time[4]-time[3], time[5]-time[4], time[6]-time[5], time[7]-time[6], time[8]-time[7]);
}

void DirectComputePotential(Data *data, double *u) {
    int i, j;
    double dx, dy;
    for(i = 0; i < data->N; i++) {
        u[i] = 0.0;
        for(j = 0; j < data->N; j++) {
            if(i != j) {
                // add contribution to ui temporary variable
                dx = data->x[i] - data->x[j];
                dy = data->y[i] - data->y[j];
                u[i] += data->q[i] / sqrt(dx*dx + dy*dy);
            }
        }
    }
}

void ClearPotential(Data *data) {
    for(int i = 0; i < data->N; i++) {
        data->u[i] = 0.0;
    }
    for(int i = 0; i < data->Ncells + data->levelOffset[data->level]; i++) {
        data->L[i] = 0.0;
        data->M[i] = 0.0;
    }
}

int getIndex(int iX[2], int level) {
    int l, index = 0;
    int jX[2] = {iX[0], iX[1]};
    for(l = 0; l < level; l++) {
        index += (jX[1] & 1) << (2*l); // least significant bit to the power of 2*l
        // least significan bit because!
        // 1100101 & 0000001 = 0000001
        jX[1] >>= 1;
        index += (jX[0] & 1) << (2*l+1);
        jX[0] >>= 1;
    }
    return index;
}

void getIX(int iX[2], int index) {
    int l = 0;
    iX[0] = iX[1] = 0;
    while(index > 0) {
        iX[1] += (index & 1) << l;
        index >>= 1;
        iX[0] += (index & 1) << l;
        index >>= 1;
        l++;
    }
}


void ParticleSortTreeSetup(Data *data) {
    data->cells = (int *) malloc(sizeof(int)*(data->level + 1)*pow(4, data->level));
    data->offset = (int *) malloc(sizeof(int)*pow(4, data->level));
    int i, iX[2];
    int iMax = data->icell[0];
    for (i=1; i<data->N; i++) {
        if (iMax < data->icell[i]) {
            iMax = data->icell[i];
        }
    }
    iMax++;
    int *bucket = malloc(iMax*sizeof(int));
    double *x2, *y2, *q2;
    x2 = (double *) malloc(sizeof(double) * data->N);
    y2 = (double *) malloc(sizeof(double) * data->N);
    q2 = (double *) malloc(sizeof(double) * data->N);
    for (i=0; i<iMax; i++) {
        bucket[i] = 0;
    }
    for (i=0; i<data->N; i++) {
        bucket[data->icell[i]]++;
    }
    for (i=1; i<iMax; i++) {
        bucket[i] += bucket[i-1];
    }
    for (i=data->N-1; i>=0; i--) {
        bucket[data->icell[i]]--;
        int inew = bucket[data->icell[i]];
        x2[inew] = data->x[i];
        y2[inew] = data->y[i];
        q2[inew] = data->q[i];
    }

    for (i=0; i<data->N; i++) {
        data->x[i] = x2[i];
        data->y[i] = y2[i];
        data->q[i] = q2[i]; 
    }
    data->Nleaf = 0;
    int ic = -1;
    for(i = 0; i < data->N; i++) {
        iX[0] = data->x[i] * data->gridsize;
        iX[1] = data->y[i] * data->gridsize;
        data->icell[i] = getIndex(iX, data->level);
        if (ic != data->icell[i]) {
            ic = data->icell[i];
            data->cells[data->level*data->Ncells + data->Nleaf] = ic;
            // Calculate offset of particles in each cell.
            data->offset[data->Nleaf] = i;
            data->Nleaf++;
        }
    }
    data->offset[data->Nleaf] = data->N;



    data->levelOffset = (int *) malloc(sizeof(int)*(data->level+1));
    data->levelOffset[0] = 0;
    for(int l = 0; l < data->level; l++) {
        data->levelOffset[l+1] = data->levelOffset[l] + pow(4, l); // Need to change 4 here to 8 for 3-D
    }

    data->M = (double *) malloc(sizeof(double)*(data->Ncells + data->levelOffset[data->level]));
    data->L = (double *) malloc(sizeof(double)*(data->Ncells + data->levelOffset[data->level]));



    for(int i = 0; i < data->Ncells + data->levelOffset[data->level]; i++) {
            data->M[i] = 0;
            data->L[i] = 0;
    }

















            
    free(bucket);
    free(x2);
    free(y2);
    free(q2);
}
