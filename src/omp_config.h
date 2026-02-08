#ifndef OMP_CONFIG_H
#define OMP_CONFIG_H

#ifdef _OPENMP
    #include <omp.h>
#else
    inline int omp_get_max_threads() { return 1; }
    inline int omp_get_thread_num() { return 0; }
    inline int omp_get_num_threads() { return 1; }
    inline void omp_set_num_threads(int) {}
#endif

#endif // OMP_CONFIG_H
