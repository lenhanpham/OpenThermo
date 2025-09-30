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

#include "chemsys.h"
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
    /**
     * @brief Detects symmetry elements and operations in a molecule
     *
     * This function analyzes the molecular structure to identify symmetry elements
     * such as rotation axes, reflection planes, and inversion centers. It determines
     * the point group by finding all symmetry operations and their associated elements.
     *
     * @param natoms Number of atoms in the molecule
     * @param nat Vector of atomic numbers for each atom
     * @param coord 3D coordinates of atoms (vector of vectors, size 3 x natoms)
     * @param symb Vector of element symbols (unused in current implementation)
     * @param delta Tolerance for symmetry detection (distance threshold)
     * @param ng Output: Total number of symmetry operations (group order)
     * @param ni Output: Number of inversion operations (0 or 1)
     * @param nsg Output: Number of symmetry planes (mirror planes)
     * @param ncr Output: Number of proper rotation axes
     * @param nsr Output: Number of improper rotation axes
     * @param np Output: Order of the principal rotation axis
     * @param symn Output: Symmetry element normals/directions (3 x max_sym)
     * @param nsym Output: Symmetry operation details (max_sym x 5)
     * @param nout Output control level (0: silent, 1: basic, 2: detailed)
     * @param nprm Output: Number of permutations found
     * @param nper Output: Permutation matrix (natoms x max_perm)
     * @param nseq Output: Number of equivalence classes
     * @param nccl Output: Size of each equivalence class
     * @param nscl Output: Members of each equivalence class
     */

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

    /**
     * @brief Classify atoms into equivalence classes based on symmetry permutations
     *
     * This function groups atoms into equivalence classes where atoms in the same class
     * are related by symmetry operations. It uses the permutation matrix to determine
     * which atoms are equivalent under the symmetry group.
     *
     * @param natoms Number of atoms in the molecule
     * @param nprm Number of symmetry permutations
     * @param nper Permutation matrix (natoms x nprm) where nper[i][j] is the atom mapped to i by permutation j
     * @param nseq Output: Number of equivalence classes found
     * @param nccl Output: Vector containing the size of each equivalence class
     * @param nscl Output: 2D vector where nscl[k][c] contains the atoms in class c (0-based indexing)
     * @param nat Vector of atomic numbers (currently unused)
     * @param symb Vector of element symbols (currently unused)
     * @param nout Output verbosity level (currently unused)
     */

    void symclass(int                                  natoms,
                  int                                  nprm,
                  const std::vector<std::vector<int>>& nper,
                  int&                                 nseq,
                  std::vector<int>&                    nccl,
                  std::vector<std::vector<int>>&       nscl,
                  const std::vector<int>&              nat,
                  const std::vector<std::string>&      symb,
                  int                                  nout);


    /**
     * @brief Performs a proper rotation around an axis
     *
     * This function applies a rotation transformation to the molecule around the specified axis
     * by the given angle and checks which atoms are mapped to equivalent positions.
     *
     * @param natoms Number of atoms in the molecule
     * @param nat Vector of atomic numbers for each atom
     * @param coord 3D coordinates of atoms (vector of vectors, size 3 x natoms)
     * @param axis Unit vector defining the rotation axis
     * @param sina Sine of the rotation angle
     * @param cosa Cosine of the rotation angle
     * @param delta Tolerance for symmetry detection (distance threshold)
     * @param nc Output: Number of atoms that are properly mapped by this rotation
     * @param ntrans Output: Permutation array showing atom mappings (size natoms)
     * @param delta3 Output: Maximum deviation from perfect symmetry
     */

    void symm_rotate(int                                     natoms,
                     const std::vector<int>&                 nat,
                     const std::vector<std::vector<double>>& coord,
                     const std::array<double, 3>&            axis,
                     double                                  sina,
                     double                                  cosa,
                     double                                  delta,
                     int&                                    nc,
                     std::vector<int>&                       ntrans,
                     double&                                 delta3);

    /**
     * @brief Performs a reflection through a plane
     *
     * This function applies a reflection transformation through the specified plane
     * and determines which atoms are mapped to equivalent positions.
     *
     * @param natoms Number of atoms in the molecule
     * @param nat Vector of atomic numbers for each atom
     * @param coord 3D coordinates of atoms (vector of vectors, size 3 x natoms)
     * @param normal Unit vector normal to the reflection plane
     * @param point A point on the reflection plane
     * @param delta Tolerance for symmetry detection (distance threshold)
     * @param nc Output: Number of atoms that are properly mapped by this reflection
     * @param ntrans Output: Permutation array showing atom mappings (size natoms)
     * @param delta3 Output: Maximum deviation from perfect symmetry
     */

    void symm_reflect(int                                     natoms,
                      const std::vector<int>&                 nat,
                      const std::vector<std::vector<double>>& coord,
                      const std::array<double, 3>&            normal,
                      const std::array<double, 3>&            point,
                      double                                  delta,
                      int&                                    nc,
                      std::vector<int>&                       ntrans,
                      double&                                 delta3);

    /**
     * @brief Performs an improper rotation (rotation followed by reflection)
     *
     * This function applies an improper rotation transformation (rotation around an axis
     * followed by reflection through a plane perpendicular to that axis) and checks
     * which atoms are mapped to equivalent positions.
     *
     * @param natoms Number of atoms in the molecule
     * @param nat Vector of atomic numbers for each atom
     * @param coord 3D coordinates of atoms (vector of vectors, size 3 x natoms)
     * @param axis Unit vector defining the rotation axis
     * @param sina Sine of the rotation angle
     * @param cosa Cosine of the rotation angle
     * @param delta Tolerance for symmetry detection (distance threshold)
     * @param nc Output: Number of atoms that are properly mapped by this improper rotation
     * @param ntrans Output: Permutation array showing atom mappings (size natoms)
     * @param delta3 Output: Maximum deviation from perfect symmetry
     */

    void symm_srotate(int                                     natoms,
                      const std::vector<int>&                 nat,
                      const std::vector<std::vector<double>>& coord,
                      const std::array<double, 3>&            axis,
                      double                                  sina,
                      double                                  cosa,
                      double                                  delta,
                      int&                                    nc,
                      std::vector<int>&                       ntrans,
                      double&                                 delta3);

    /**
     * @brief Computes the dot product of two vectors
     *
     * Calculates the scalar (dot) product of two n-dimensional vectors.
     *
     * @param a First vector
     * @param b Second vector
     * @param n Dimension of the vectors
     * @return Dot product value
     */

    auto symm_dot(const double* a, const double* b, int n) -> double;
    /**
     * @brief Computes the cross product of two 3D vectors
     *
     * Calculates the vector (cross) product of two 3D vectors x and y, storing the result in z.
     *
     * @param x First 3D vector
     * @param y Second 3D vector
     * @param z Output: Cross product result
     */

    void symm_crossp(const std::array<double, 3>& x, const std::array<double, 3>& y, std::array<double, 3>& z);
    /**
     * @brief Computes the greatest common divisor of two integers
     *
     * Uses the Euclidean algorithm to find the greatest common divisor (GCD) of two integers.
     *
     * @param a First integer
     * @param b Second integer
     * @return Greatest common divisor of a and b
     */

    auto symm_igcd(int a, int b) -> int;

    /**
     * @brief Checks if transformed coordinates match original within tolerance
     *
     * Compares the transformed coordinates (cord) with the original coordinates (coord)
     * to determine which atoms are mapped to equivalent positions within the given tolerance.
     *
     * @param natoms Number of atoms in the molecule
     * @param delta Tolerance for symmetry detection (distance threshold)
     * @param nat Vector of atomic numbers for each atom
     * @param coord Original 3D coordinates of atoms (vector of vectors, size 3 x natoms)
     * @param cord Transformed 3D coordinates of atoms (vector of vectors, size 3 x natoms)
     * @param nc Output: Number of atoms that are properly mapped
     * @param ntrans Output: Permutation array showing atom mappings (size natoms)
     * @param delta3 Output: Maximum deviation from perfect symmetry
     */

    void symm_check(int                                     natoms,
                    double                                  delta,
                    const std::vector<int>&                 nat,
                    const std::vector<std::vector<double>>& coord,
                    const std::vector<std::vector<double>>& cord,
                    int&                                    nc,
                    std::vector<int>&                       ntrans,
                    double&                                 delta3);

    /**
     * @brief Calculates the center of mass of a molecule
     *
     * Computes the center of mass coordinates and total molecular weight using atomic masses.
     *
     * @param natoms Number of atoms in the molecule
     * @param nat Vector of atomic numbers for each atom
     * @param wt Vector of atomic weights for each atom
     * @param coord 3D coordinates of atoms (vector of vectors, size 3 x natoms)
     * @param wmol Output: Total molecular weight
     * @param cmx Output: X-coordinate of center of mass
     * @param cmy Output: Y-coordinate of center of mass
     * @param cmz Output: Z-coordinate of center of mass
     */

    void symm_cmass(int                                     natoms,
                    const std::vector<int>&                 nat,
                    const std::vector<double>&              wt,
                    const std::vector<std::vector<double>>& coord,
                    double&                                 wmol,
                    double&                                 cmx,
                    double&                                 cmy,
                    double&                                 cmz);

    /**
     * @brief Shifts coordinates to center them at the origin
     *
     * Translates all atomic coordinates by subtracting the center coordinates (pc).
     *
     * @param natoms Number of atoms in the molecule
     * @param coord Input/Output: 3D coordinates of atoms (modified in place)
     * @param pc Center coordinates to subtract (size 3)
     */

    void symm_cshift(int natoms, std::vector<std::vector<double>>& coord, const std::vector<double>& pc);

    /**
     * @brief Adds a proper rotation axis to the symmetry element list
     *
     * Records a Cn rotation axis if it doesn't already exist in the list,
     * or updates the order if a higher-order rotation is found for the same axis.
     *
     * @param nrot Input/Output: Current number of rotation axes (incremented if added)
     * @param rotn Input/Output: List of rotation axis vectors
     * @param rota Input/Output: List of rotation angles corresponding to axes
     * @param axis Unit vector defining the rotation axis
     * @param point A point on the rotation axis (unused in current implementation)
     * @param angle Rotation angle in radians
     * @param delta Tolerance for comparing axis directions
     */

    void add_Cn(int&                                nrot,
                std::vector<std::array<double, 3>>& rotn,
                std::vector<double>&                rota,
                const std::array<double, 3>&        axis,
                const std::array<double, 3>&        point,
                double                              angle,
                double                              delta);
    /**
     * @brief Adds a symmetry plane (mirror plane) to the list
     *
     * Records a mirror plane if it doesn't already exist in the list of symmetry planes.
     *
     * @param nsg Input/Output: Current number of symmetry planes (incremented if added)
     * @param sigman Input/Output: List of plane normal vectors
     * @param normal Unit vector normal to the mirror plane
     * @param delta Tolerance for comparing plane normals
     */

    void
    add_SG(int& nsg, std::vector<std::array<double, 3>>& sigman, const std::array<double, 3>& normal, double delta);
    /**
     * @brief Adds a permutation to the symmetry operation list
     *
     * Records a new permutation if it doesn't already exist in the list of symmetry permutations.
     *
     * @param natoms Number of atoms in the molecule
     * @param ntrans Permutation array (atom mapping)
     * @param nprm Input/Output: Current number of permutations (incremented if added)
     * @param nper Input/Output: 2D array of permutations (natoms x max_permutations)
     */

    void add_perm(int natoms, const std::vector<int>& ntrans, int& nprm, std::vector<std::vector<int>>& nper);

    /**
     * @brief Determines the point group of a molecule using symmetry analysis
     *
     * This function performs a complete point group determination by analyzing the molecular
     * structure for symmetry elements. It calculates the center of mass, shifts coordinates
     * to the origin, identifies all symmetry operations, determines equivalence classes,
     * and assigns the appropriate point group label.
     *
     * @param natoms Number of atoms in the molecule
     * @param nat Vector of atomic numbers for each atom
     * @param coord Input/Output: 3D coordinates of atoms (modified to center at origin)
     * @param delta Tolerance for symmetry detection (distance threshold)
     * @param PGlab Output: Detected point group label (e.g., "C1", "C2v", "D3h", etc.)
     */

    void PG_determ(int                               natoms,
                   const std::vector<int>&           nat,
                   std::vector<std::vector<double>>& coord,
                   double                            delta,
                   std::string&                      PGlab);

    /**
     * @brief Determines the point group label from symmetry element counts
     *
     * This function matches the provided symmetry element counts (group order, inversions,
     * planes, rotations, etc.) against a database of known point groups to determine
     * the appropriate point group label. It also provides detailed output about the
     * point group properties when requested.
     *
     * @param ngp Total number of symmetry operations (group order)
     * @param ni Number of inversion operations (0 or 1)
     * @param nsg Number of symmetry planes (mirror planes)
     * @param ncr Number of proper rotation axes
     * @param nsr Number of improper rotation axes
     * @param np Order of the principal rotation axis
     * @param nout Output verbosity level (0: silent, >=1: print diagnostics)
     * @return Point group label string (e.g., "C2v", "D3h", etc.)
     */

    auto symm_point_group(int ngp, int ni, int nsg, int ncr, int nsr, int np, int nout) -> std::string;


    // Additional helper functions that would be called by sym_elements


    /**
     * @brief Static data structure containing point group database and symmetry information
     *
     * This struct provides access to comprehensive point group data including
     * symmetry operation counts, subgroup relationships, and point group classifications.
     * All members are static constants or arrays initialized at compile time.
     */
    struct SymmetryData
    {
        static const int max_pgs       = 57;  /**< Maximum number of supported point groups */
        static const int max_subgroups = 406; /**< Maximum number of subgroup relationships */

        static std::array<std::string, max_pgs>        pgsymb; /**< Point group symbols (extended form) */
        static std::array<std::array<int, 2>, max_pgs> nsgb;   /**< Subgroup boundary indices for each point group */
        static std::array<int, max_subgroups>          nsgr;   /**< Subgroup relationship indices */
        static const std::array<std::string, max_pgs>  sg;     /**< Point group symbols (short form) */
        static const std::array<std::array<int, 6>, max_pgs>
            ng; /**< Symmetry operation counts: [order, inversion, planes, proper_rotations, improper_rotations,
                   principal_axis] */
        static const std::array<std::string, max_pgs> cg; /**< Point group descriptions and generator information */
    };

    struct Vector3D
    {
        double x, y, z; /**< Cartesian coordinates */

        /** Default constructor */
        Vector3D() = default;

        /** Constructor with coordinates */
        Vector3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

        /** Constructor from std::array */
        Vector3D(const std::array<double, 3>& arr) : x(arr[0]), y(arr[1]), z(arr[2]) {}

        /** Get pointer to data for external library compatibility */
        [[nodiscard]] auto data() -> double*
        {
            return &x;
        }

        /** Get const pointer to data */
        [[nodiscard]] auto data() const -> const double*
        {
            return &x;
        }
    };


    /**
     * @brief Performs an inversion through the origin
     *
     * Applies the inversion operation (x,y,z) -> (-x,-y,-z) and checks which atoms
     * are mapped to equivalent positions.
     *
     * @param natoms Number of atoms in the molecule
     * @param nat Vector of atomic numbers for each atom
     * @param coord 3D coordinates of atoms (vector of vectors, size 3 x natoms)
     * @param delta Tolerance for symmetry detection (distance threshold)
     * @param nc Output: Number of atoms that are properly mapped by inversion
     * @param ntrans Output: Permutation array showing atom mappings (size natoms)
     * @param delta3 Output: Maximum deviation from perfect symmetry
     */

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
    class SymmetryDetector
    {
    public:
        // Data members
        int               ncenter = 0;       /**< Number of atoms */
        std::vector<Atom> a;                 /**< Array of atoms */
        std::vector<int>  a_index;           /**< Atomic indices */
        std::string       PGnameinit = "?"; /**< Initial point group label */
        std::string       PGname     = "?"; /**< Detected point group label */
        int               rotsym      = 1;   /**< Rotational symmetry number */

        // Methods
        /**
         * @brief Detects the point group of the molecule
         *
         * Analyzes the atomic coordinates and types stored in the class members
         * to determine the molecular point group using systematic symmetry element detection.
         *
         * @param ishow Output verbosity level (0: silent, 1: show results)
         */
        void detectPG(int ishow = 0);

        /**
         * @brief Converts point group name to rotational symmetry number
         *
         * Based on the detected point group name (PGname), determines the rotational
         * symmetry number (rotsym) used in spectroscopic calculations.
         */
        void PGname2rotsym();
    };

}  // namespace symmetry


#endif  // SYMMETRY_H