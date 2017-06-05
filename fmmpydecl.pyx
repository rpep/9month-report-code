cdef extern from "fmm.h":
    ctypedef struct Data:
        double *u
        double *x
        double *y
        double *M
        double *L
        double *q
        int *icell
        int *N
        int *level
        int *gridsize

    int AllocateData(Data *data, int N, int level)
    void InitTestData(Data *data)
    void PrintData(Data *data)
    void DeallocateData(Data *data)
    void ParticleToMultipole(Data *data)
    void MultipoleToMultipole(Data *data)
    void MultipoleToLocal(Data *data)
    void LocalToLocal(Data *data)
    void ClearPotential(Data *data)
    void LocalToParticle(Data *data)
    void ParticleToParticle(Data *data)
    void FMMComputePotential(Data *data)
    void DirectComputePotential(Data *data, double *u)
    void ParticleSortTreeSetup(Data *data)