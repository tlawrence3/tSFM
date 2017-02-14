#cython: boundscheck=False, cdivision=True, wraparound=False, nonecheck=False
from libc.math cimport pow, log, exp, log2
from libc.stdlib cimport malloc, free

cdef double Ramanujan_logfact(double n, double logpi):
    cdef double lf
    if (n == 0 or n == 1):
        return 0.0
    lf = (n * log(n)) - n
    lf += log(n + (4 * n**2) + (8 * n**3))/6
    lf += logpi
    return lf

def calc_exact(double n, list p, int numclasses):
    cdef double exact, log_probability, entropy
    cdef int j, y, r, i
    cdef double *counts
    cdef double *pC
    cdef double PI = 3.14159265358979323
    cdef double HALF_LOG_PI = log(PI)/2.0
    pC = <double *>malloc(numclasses * sizeof(double))
    counts = <double *>malloc(numclasses * sizeof(double))
    if counts is NULL or pC is NULL:
        raise MemoryError()
    
    counts[0] = n
    for i in range(1,numclasses):
        counts[i] = 0
    for i in range(numclasses):
        pC[i] = p[i]

    r = 0 ##first index for which counts >= 1
    exact = 0.0
    while True:
        log_probability = Ramanujan_logfact(n, HALF_LOG_PI)
        entropy = 0.0
        for j in range(numclasses):
            log_probability -= Ramanujan_logfact(counts[j], HALF_LOG_PI)
            log_probability += counts[j] * log(pC[j])
            if (counts[j] != 0):
                entropy -= (counts[j]/n) * (log2(counts[j]/n))
        exact += (exp(log_probability) * entropy)
    
        if (counts[0] != 0):
            counts[0] -= 1
            r = 1
        else:
            if (r == numclasses - 1):
                break
            else:
                counts[0] = counts[r] - 1
                counts[r] = 0
                r += 1
        counts[r] = counts[r] + 1
    
    free(counts)
    free(pC)
    return (n, exact) 
