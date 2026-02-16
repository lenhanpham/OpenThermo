/**
 * @file omp_config.h
 * @brief OpenMP configuration, thread detection, and parallelization strategy
 * @author Romain Despoullains
 * @author Le Nhan Pham
 * @date 2026
 *
 * This header provides the OpenMP threading infrastructure for OpenThermo:
 * - Physical CPU core detection (Linux, macOS, Windows)
 * - HPC job scheduler CPU allocation detection (SLURM, PBS, SGE, LSF)
 * - Thread count validation and clamping with user-friendly notifications
 * - Automatic parallelization strategy selection (outer T/P scan vs inner
 *   vibrational loop) based on workload characteristics
 * - OpenMP runtime configuration with nested parallelism control
 *
 * All functions are header-only (inline) to avoid link-time overhead.
 * When OpenMP is not available, stub functions are provided so the rest
 * of the codebase compiles and runs single-threaded without #ifdef guards.
 */

#ifndef OMP_CONFIG_H
#define OMP_CONFIG_H

#include <thread>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifdef _OPENMP
    #include <omp.h>
#else
    inline int  omp_get_max_threads()       { return 1; }
    inline int  omp_get_thread_num()        { return 0; }
    inline int  omp_get_num_threads()       { return 1; }
    inline void omp_set_num_threads(int)    {}
    inline int  omp_get_level()             { return 0; }
#endif

#ifdef _WIN32
    #ifndef WIN32_LEAN_AND_MEAN
        #define WIN32_LEAN_AND_MEAN
    #endif
    #include <windows.h>
#endif

/**
 * @brief Check HPC job scheduler environment for allocated CPU count
 *
 * Checks common scheduler env vars in priority order:
 *   SLURM_CPUS_PER_TASK, PBS_NUM_PPN, PBS_NP, NSLOTS, LSB_DJOB_NUMPROC
 * These represent the allocated resources, which are the correct ceiling
 * on HPC systems (not the total hardware cores).
 *
 * @return Allocated CPU count, or 0 if no scheduler detected
 */
inline int detect_scheduler_cpus()
{
    // SLURM (most common HPC scheduler)
    const char* env = std::getenv("SLURM_CPUS_PER_TASK");
    if (env)
    {
        int cpus = std::atoi(env);
        if (cpus > 0) return cpus;
    }

    // PBS/Torque
    env = std::getenv("PBS_NUM_PPN");
    if (env)
    {
        int cpus = std::atoi(env);
        if (cpus > 0) return cpus;
    }
    env = std::getenv("PBS_NP");
    if (env)
    {
        int cpus = std::atoi(env);
        if (cpus > 0) return cpus;
    }

    // SGE / Univa Grid Engine
    env = std::getenv("NSLOTS");
    if (env)
    {
        int cpus = std::atoi(env);
        if (cpus > 0) return cpus;
    }

    // IBM LSF
    env = std::getenv("LSB_DJOB_NUMPROC");
    if (env)
    {
        int cpus = std::atoi(env);
        if (cpus > 0) return cpus;
    }

    return 0;  // No scheduler detected
}

/**
 * @brief Detect physical CPU cores (best effort)
 *
 * Uses platform-specific methods to detect physical (not logical/HT) cores.
 * Falls back to std::thread::hardware_concurrency() / 2.
 *
 * @return Physical core count (minimum 1)
 */
inline int detect_physical_cores()
{
#if defined(__linux__)
    // Linux: parse lscpu for physical cores
    FILE* fp = popen("lscpu -p 2>/dev/null | grep -v '^#' | awk -F',' '{print $2}' | sort -u | wc -l", "r");
    if (fp)
    {
        int cores = 0;
        if (fscanf(fp, "%d", &cores) == 1 && cores > 0)
        {
            pclose(fp);
            return cores;
        }
        pclose(fp);
    }
#elif defined(__APPLE__)
    // macOS: sysctl physicalcpu
    FILE* fp = popen("sysctl -n hw.physicalcpu 2>/dev/null", "r");
    if (fp)
    {
        int cores = 0;
        if (fscanf(fp, "%d", &cores) == 1 && cores > 0)
        {
            pclose(fp);
            return cores;
        }
        pclose(fp);
    }
#elif defined(_WIN32)
    // Windows: GetLogicalProcessorInformation
    DWORD returnLength = 0;
    GetLogicalProcessorInformation(nullptr, &returnLength);
    if (returnLength > 0)
    {
        std::vector<SYSTEM_LOGICAL_PROCESSOR_INFORMATION> buffer(
            returnLength / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION));
        if (GetLogicalProcessorInformation(buffer.data(), &returnLength))
        {
            int cores = 0;
            DWORD count = returnLength / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
            for (DWORD i = 0; i < count; ++i)
            {
                if (buffer[i].Relationship == RelationProcessorCore)
                {
                    ++cores;
                }
            }
            if (cores > 0) return cores;
        }
    }
#endif

    // Fallback: assume hyperthreading, use logical / 2
    int logical = static_cast<int>(std::thread::hardware_concurrency());
    if (logical > 1)
    {
        return logical / 2;
    }
    return 1;  // Conservative minimum
}

/**
 * @brief Calculate default thread count (half physical cores, min 1)
 * @param physical_cores Number of physical cores detected
 * @return Default thread count
 */
inline int get_default_threads(int physical_cores)
{
    return std::max(1, physical_cores / 2);
}

/**
 * @brief Validate and clamp the requested thread count
 *
 * Logic:
 *   - requested <= 0 : auto mode, use half of physical cores
 *   - requested > physical_cores : clamp to default (half), warn user
 *   - otherwise : respect the explicit user request
 *
 * @param requested  User-requested thread count (0 = auto)
 * @param physical_cores Detected physical cores
 * @param actual     [out] Actual thread count to use
 * @param user_override [out] Whether the user explicitly requested threads
 * @return Notification/warning message for the user
 */
inline std::string validate_thread_count(int requested, int physical_cores,
                                         int scheduler_cpus,
                                         int& actual, bool& user_override)
{
    user_override = (requested > 0);

    // Determine the effective core ceiling and default
    // If a scheduler allocated CPUs, use that as the ceiling instead of hardware
    bool has_scheduler = (scheduler_cpus > 0);
    int  effective_cores = has_scheduler ? scheduler_cpus : physical_cores;
    int  default_threads = get_default_threads(effective_cores);

    std::string core_source = has_scheduler
        ? std::to_string(scheduler_cpus) + " scheduler-allocated CPUs"
        : std::to_string(physical_cores) + " physical cores";

    if (!user_override)
    {
        // Auto mode: use conservative default
        actual = default_threads;
        return "OpenMP threads: " + std::to_string(actual) +
               " (default: half of " + core_source +
               "). Use -omp-threads N to override.";
    }

    if (requested > effective_cores)
    {
        // Clamp to default — safety guardrail
        actual = default_threads;
        return "Warning: Requested " + std::to_string(requested) +
               " threads exceeds " + core_source +
               ". Using default (" + std::to_string(actual) +
               "). Use -omp-threads <= " + std::to_string(effective_cores) +
               " to use more threads.";
    }

    // Valid explicit request — respect it
    actual = requested;
    return "OpenMP threads: " + std::to_string(actual) + " (user-specified)";
}

/**
 * @brief Configure OpenMP runtime with validated thread count
 *
 * Disables nested parallelism and sets the thread count.
 * Should be called once at program startup after validation.
 *
 * @param num_threads Number of threads to use (must be >= 1)
 */
inline void configure_openmp(int num_threads)
{
#ifdef _OPENMP
    // Disable nested parallelism (single active level only)
    #if _OPENMP >= 201811  // OpenMP 5.0+
        omp_set_max_active_levels(1);
    #else
        omp_set_nested(0);
    #endif

    // Set thread count
    if (num_threads > 0)
    {
        omp_set_num_threads(num_threads);
    }
#else
    (void)num_threads;
#endif
}

/**
 * @brief OpenMP parallelization strategy
 *
 * Outer: parallelize the T/P scan loop (main.cpp), vibrational loops run serially.
 * Inner: T/P scan runs serially, vibrational loop (calc.cpp) is parallelized.
 */
enum class OMPStrategy
{
    Outer = 0,  ///< Parallelize T/P scan, serial vibrational
    Inner = 1   ///< Serial T/P scan, parallel vibrational
};

/**
 * @brief Auto-select the optimal parallelization strategy
 *
 * @param total_points Number of T/P grid points to calculate
 * @param nfreq        Number of vibrational frequencies
 * @param nthreads     Available OpenMP threads
 * @return Selected strategy
 */
inline OMPStrategy select_strategy(int total_points, int nfreq, int nthreads)
{
    if (total_points >= nthreads && total_points > 1)
    {
        // Many T/P points: parallelize outer loop (simpler, better cache)
        return OMPStrategy::Outer;
    }
    else if (nfreq > 50)
    {
        // Few T/P points but many frequencies: parallelize inner vibrational loop
        return OMPStrategy::Inner;
    }
    else
    {
        // Neither benefits strongly; default to outer (simpler)
        return OMPStrategy::Outer;
    }
}

/**
 * @brief Get a human-readable description of the selected strategy
 */
inline std::string strategy_description(OMPStrategy strategy, int total_points, int nfreq)
{
    if (strategy == OMPStrategy::Outer)
    {
        return "Strategy: outer (parallel T/P scan, " + std::to_string(total_points) +
               " points, " + std::to_string(nfreq) + " frequencies)";
    }
    else
    {
        return "Strategy: inner (serial T/P scan, parallel vibrational, " +
               std::to_string(nfreq) + " frequencies)";
    }
}

#endif // OMP_CONFIG_H
