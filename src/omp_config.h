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

/**
 * @brief Configure OpenMP settings for the application
 *
 * Disables nested parallelism to prevent thread oversubscription
 * when outer parallel loops (e.g., T/P scanning) call functions
 * that contain their own parallel regions (e.g., calcthermo).
 *
 * Should be called once at program startup.
 */
inline void configure_openmp() {
#ifdef _OPENMP
    #if _OPENMP >= 201811  // OpenMP 5.0+
        omp_set_max_active_levels(1);
    #else
        omp_set_nested(0);
    #endif
#endif
}

#endif // OMP_CONFIG_H
