from libc.stdlib cimport malloc, free
import sys
from timeit import default_timer as timer
import sys

include "fmmpydecl.pyx"

def printf(*args):
    print(args)
    sys.stdout.flush()




cdef int N = 10000

cdef double *u = <double*>malloc(N*sizeof(double))

cdef Data data

AllocateData(&data, N, 2)
InitTestData(&data)

start = timer()

ParticleSortTreeSetup(&data)
FMMComputePotential(&data)

middle = timer()
DirectComputePotential(&data, u)
printf('Direct Done')
end = timer()
printf('FMM {} Direct {}'.format(middle - start, end-middle))
printf('{}'.format((end-middle)/(middle-start)))

#for i in range(N):
#    print('{} {} {} {}'.format(i, data.icell[i], u[i], data.u[i]))


print 'N\t{}\nFMM\t{}\nDirect\t{}\nspeedup\t{}\n'.format(N,middle - start, end-middle, (end-middle)/(middle-start))

sys.stdout.flush()




DeallocateData(&data)
free(u)
