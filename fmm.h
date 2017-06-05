#ifndef __FMM_H__
#define __FMM_H__
#include<stdlib.h>
#include<omp.h>
typedef struct Data {
    double *x, *y, *u, *q, *M, *L;
    int N, level, gridsize, Ncells, Nleaf, *levelOffset, *icell, *cells, *offset;

} Data;

void AllocateData(Data *data, int N, int level);
void InitTestData(Data *data);
void PrintData(Data *data);
void DeallocateData(Data *data);
void ParticleToMultipole(Data *data);
void MultipoleToLocal(Data *data);
void LocalToParticle(Data *data);
void ParticleToParticle(Data *data);
void FMMComputePotential(Data *data);
void DirectComputePotential(Data *data, double *u);
void ClearPotential(Data *data);
int getIndex(int iX[2], int level);
void getIX(int iX[2], int index);
void ParticleSortTreeSetup(Data *data);

#endif
