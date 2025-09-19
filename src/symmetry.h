/**
 * @file symmetry.h
 * @brief Header for molecular symmetry detection and analysis
 * @author Le Nhan Pham
 * @date 2025
 *
 * This header file declares classes and functions for detecting molecular
 * point groups, analyzing symmetry elements, and performing symmetry-based
 * calculations. It includes implementations of the SYVA algorithm for
 * automatic point group determination.
 */

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "defvar.h"
#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @brief Namespace containing symmetry detection and analysis functions
 *
 * This namespace encapsulates all functions and classes related to molecular
 * symmetry analysis, point group determination, and symmetry element detection.
 */
namespace symmetry

{
    void sym_elements(int                                     natoms,
                        const std::vector<int>&                 nat,
                        const std::vector<std::vector<double>>& coord,
                        const std::vector<std::string>&         symb,
                        double                                  delta,
                        int&                                    ng,
                        int&                                    ni,
                        int&                                    nsg,
                        int&                                    ncr,
                        int&                                    nsr,
                        int&                                    np,
                        std::vector<std::vector<double>>&       symn,
                        std::vector<std::vector<int>>&          nsym,
                        int                                     nout,
                        int&                                    nprm,
                        std::vector<std::vector<int>>&          nper,
                        int&                                    nseq,
                        std::vector<int>&                       nccl,
                        std::vector<std::vector<int>>&          nscl);

    void symclass(int natoms, int nprm, const std::vector<std::vector<int>>& nper,
                  int& nseq, std::vector<int>& nccl, std::vector<std::vector<int>>& nscl,
                  const std::vector<int>& nat, const std::vector<std::string>& symb, int nout);



    void symm_rotate(int natoms, const std::vector<int>& nat,
                    const std::vector<std::vector<double>>& coord,
                    const std::array<double, 3>& axis, double sina, double cosa,
                    double delta, int& nc, std::vector<int>& ntrans, double& delta3);

    void symm_reflect(int natoms, const std::vector<int>& nat,
                     const std::vector<std::vector<double>>& coord,
                     const std::array<double, 3>& normal, const std::array<double, 3>& point,
                     double delta, int& nc, std::vector<int>& ntrans, double& delta3);

    void symm_srotate(int natoms, const std::vector<int>& nat,
                     const std::vector<std::vector<double>>& coord,
                    const std::array<double, 3>& axis, double sina, double cosa,
                     double delta, int& nc, std::vector<int>& ntrans, double& delta3);

    double symm_dot(const double* a, const double* b, int n);
    void symm_crossp(const std::array<double, 3>& x, const std::array<double, 3>& y, std::array<double, 3>& z);
    int symm_igcd(int a, int b);

    void symm_check(int natoms, double delta, const std::vector<int>& nat,
                   const std::vector<std::vector<double>>& coord,
                   const std::vector<std::vector<double>>& cord,
                   int& nc, std::vector<int>& ntrans, double& delta3);

    void symm_cmass(int natoms, const std::vector<int>& nat, const std::vector<double>& wt,
                   const std::vector<std::vector<double>>& coord, double& wmol,
                   double& cmx, double& cmy, double& cmz);

    void symm_cshift(int natoms, std::vector<std::vector<double>>& coord, const std::vector<double>& pc);

    void add_Cn(int& nrot, std::vector<std::array<double, 3>>& rotn, std::vector<double>& rota,
                const std::array<double, 3>& axis, const std::array<double, 3>& point, double angle, double delta);
    void add_SG(int& nsg, std::vector<std::array<double, 3>>& sigman,
                const std::array<double, 3>& normal, double delta);
    void add_perm(int natoms, const std::vector<int>& ntrans, int& nprm, 
             std::vector<std::vector<int>>& nper);















    // Additional helper functions that would be called by sym_elements






    // SymmetryData struct
    struct SymmetryData {
        static const int max_pgs = 57;
        static const int max_subgroups = 406;
        static std::array<std::string, max_pgs> pgsymb;
        static std::array<std::array<int, 2>, max_pgs> nsgb;
        static std::array<int, max_subgroups> nsgr;
        static const std::array<std::string, max_pgs> sg;
        static const std::array<std::array<int, 6>, max_pgs> ng;
        static const std::array<std::string, max_pgs> cg;
    };

    struct Vector3D {
        double x, y, z; /**< Cartesian coordinates */

        /** Default constructor */
        Vector3D() = default;

        /** Constructor with coordinates */
        Vector3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

        /** Constructor from std::array */
        Vector3D(const std::array<double, 3>& arr) : x(arr[0]), y(arr[1]), z(arr[2]) {}

        /** Get pointer to data for external library compatibility */
        double* data() { return &x; }

        /** Get const pointer to data */
        const double* data() const { return &x; }
    };


    void symm_inversion(int                                     natoms,
                        const std::vector<int>&                 nat,
                        const std::vector<std::vector<double>>& coord,
                        double                                  delta,
                        int&                                    nc,
                        std::vector<int>&                       ntrans,
                        double                                  delta3);

    /**
     * @brief Class for symmetry detection and point group analysis
     */
    class SymmetryDetector {
    public:
        // Data members
        int ncenter = 0;                        /**< Number of atoms */
        std::vector<Atom> a;                    /**< Array of atoms */
        std::vector<int> a_index;               /**< Atomic indices */
        std::string PGlabelinit = "?";          /**< Initial point group label */
        std::string PGlabel = "?";              /**< Detected point group label */
        int rotsym = 1;                         /**< Rotational symmetry number */

        // Methods
        void detectPG(int ishow = 0);           /**< Detect point group */
        void PGlabel2rotsym();                  /**< Convert point group label to rotational symmetry */
    };

}  // namespace symmetry


#endif  // SYMMETRY_H