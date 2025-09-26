/**
 * @file symmetry.cpp
 * @brief Implementation of molecular symmetry detection and analysis
 * @author Le Nhan Pham
 * @date 2025
 *
 * This file contains the implementation of classes and functions for detecting
 * molecular point groups, analyzing symmetry elements, and performing
 * symmetry-based calculations. 
 * algorithm for automatic point group determination.
 */


#include "symmetry.h"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <array>
#include <string>
#include <algorithm>
#include <numeric>
#include <sstream>



namespace symmetry {

// Forward declarations
auto initialize_nsgb() -> std::array<std::array<int, 2>, 57>;
auto initialize_nsymop() -> std::array<std::array<std::array<int, 55>, 4>, 14>;

// Data common block equivalent
struct DataCommon {
    const std::vector<double> wt;      // Atomic weights
    const std::vector<std::string> symb; // Element symbols
};





// Template function to convert Fortran column-major to C++ row-major
template<size_t Rows, size_t Cols>
auto convert_fortran_to_cpp(const std::array<int, Rows * Cols>& flat) -> std::array<std::array<int, Cols>, Rows> {
    std::array<std::array<int, Cols>, Rows> result{};
    for (size_t i = 0; i < Rows; ++i) {
        for (size_t j = 0; j < Cols; ++j) {
            result[i][j] = flat[i * Cols + j];
        }
    }
    return result;
}



// NIR 1D Fortran data (column-major order) - CORRECTED
const std::array<int, 110> nir_1d = {
    1, 1, 2, 2, 4, 2, 6, 2, 8, 2,    // C1,Cs,Ci,C2,C3,C4,C5,C6,C7,C8
    10, 3, 13, 3, 16, 4, 20, 4, 24, 5, // D2,D3,D4,D5,D6,D7,D8,C2v,C3v,C4v
    29, 4, 33, 3, 36, 5, 41, 4, 45, 6, // C5v,C6v,C7v,C8v,C2h,C3h,C4h,C5h,C6h
    51, 5, 56, 7, 63, 4, 67, 3, 70, 5, // C7h,C8h,D2h,D3h,D4h,D5h,D6h,C2v,C3v
    75, 4, 79, 6, 85, 5, 90, 7, 97, 4, // D2d,D3d,D4d,D5d,D6d,D7d,D8d,S4,S6,S8
    101, 4, 105, 6, 111, 6, 117, 8, 125, 8, // T,Th,Td,O,Oh,I,Ih,Civ,Dih
    133, 10, 143, 8, 151, 6, 157, 10, 167, 8, // Additional entries
    175, 12, 187, 10, 197, 14, 211, 5,216, 6, // More entries
    222, 7, 229, 8, 237, 9, 246, 10, 256, 11, // Final entries
    267, 3, 270, 4, 274, 5, 279, 3, 282, 6, // Corrected C2v: start=66, count=4
    288, 5, 293, 5, 298, 10, 308, 5, 313, 10 // Additional corrections
};

// 1. Number of irreducible representations (nir)
const std::array<std::array<int, 2>, 55> nir_ = convert_fortran_to_cpp<55, 2>(nir_1d);



///** @brief Atomic weights for all elements (1-90) in atomic units */
//const std::array<double, 90> atomic_weights_ = {{
//      1.00783,   4.00260,   6.94000,   9.01218, 10.81000,   12.00000,  14.00307,  15.99491,
//     18.99840,  20.17900,  22.98977,  24.30500, 26.98154,   28.08550,  30.97376,  32.06000,
//     35.45300,  39.94800,  39.09830,  40.08000, 44.95590,   47.90000,  50.94150,  51.99600,
//     54.93800,  55.84700,  58.93320,  58.71000, 63.54600,   65.38000,  69.73500,  72.59000,
//     74.92160,  78.96000,  79.90400,  83.80000, 85.46780,   87.62000,  88.90590,  91.22000,
//     92.90640,  95.94000,  98.90620, 101.07000, 102.90550, 106.40000, 107.86800, 112.41000,
//    114.82000, 118.69000, 121.75000, 127.60000, 126.90450, 131.30000, 132.90540, 137.33000,
//      0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,  0.000000,   0.00000,
//      0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000, 178.49000,
//    180.94790, 183.85000, 186.20700, 190.20000, 192.22000, 194.96480, 196.96650, 200.59000,
//    204.37000, 207.20000, 208.98040,   0.00000,   0.00000, 0.00000,     0.00000,   0.00000,
//      0.00000,   0.00000
//}};

/** @brief Updated atomic weights for all elements (1-90) in atomic units with some from IUPAC Commission on Isotopic Abundances and Atomic Weights (CIAAW)*/
constexpr std::array<double, 90> atomic_weights_ = {{
       1.00783,    4.00260,    6.94000,    9.01218,   10.81000,   12.00000,    14.00307,   15.99491,
      18.99840,   20.17900,   22.98977,   24.30500,   26.98154,   28.08550,    30.97376,   32.06000,
      35.45300,   39.94800,   39.09830,   40.08000,   44.95590,   47.90000,    50.94150,   51.99600,
      54.93800,   55.84700,   58.93320,   58.71000,   63.54600,   65.38000,    69.73500,   72.59000,
      74.92160,   78.96000,   79.90400,   83.80000,   85.46780,   87.62000,    88.90590,   91.22000,
      92.90640,   95.94000,   98.90620,  101.07000,  102.90550,  106.40000,   107.86800,  112.41000,
     114.82000,  118.69000,  121.75000,  127.60000,  126.90450,  131.30000,   132.90540,  137.33000,
    138.905470, 140.116000, 140.907660, 144.242000,   0.000000,  150.360000, 151.964000, 157.249000,
    158.925354, 162.500000, 164.930329, 167.259000, 168.934219,  173.045000, 174.966690,  178.49000,
     180.94790,  183.85000,  186.20700,  190.20000,  192.22000,   194.96480,  196.96650,  200.59000,
     204.37000,  207.20000,  208.98040,    0.00000,    0.00000,     0.00000,    0.00000,    0.00000,
       0.00000,  232.03770
}};


/** @brief Element symbols for all elements (1-90) */
const std::array<std::string, 90> element_symbols_ = {{
     "H", "He", "Li", "Be",  "B",  "C",  "N",  "O",  "F", "Ne",
    "Na", "Mg", "Al", "Si",  "P",  "S", "Cl", "Ar", " K", "Ca",
    "Sc", "Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",  "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te",  "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", " W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "  "
}};

// Initialize data_common with the arrays
DataCommon data_common = {
    std::vector<double>(atomic_weights_.begin(), atomic_weights_.end()),  // wt
    std::vector<std::string>(element_symbols_.begin(), element_symbols_.end())  // symb
};

/**
 * Calculates the centre of mass of a molecule
 */
void symm_cmass(int natoms, const std::vector<int>& nat, const std::vector<double>& wt,
                const std::vector<std::vector<double>>& coord, double& wmol,
                double& cmx, double& cmy, double& cmz) {
    double sumwx = 0.0;
    double sumwy = 0.0;
    double sumwz = 0.0;
    wmol = 0.0;
    
    for (int i = 0; i < natoms; i++) {
        int nati = nat[i];
        wmol = wmol + wt[nati - 1]; // Adjust for 0-based indexing
        sumwx = sumwx + wt[nati - 1] * coord[0][i];
        sumwy = sumwy + wt[nati - 1] * coord[1][i];
        sumwz = sumwz + wt[nati - 1] * coord[2][i];
    }
    
    cmx = sumwx / wmol;
    cmy = sumwy / wmol;
    cmz = sumwz / wmol;
}

/**
 * Shifts the centre of mass of a molecule to the origin
 */
void symm_cshift(int natoms, std::vector<std::vector<double>>& coord, const std::vector<double>& pc) {
    for (int i = 0; i < natoms; i++) {
        coord[0][i] = coord[0][i] - pc[0];
        coord[1][i] = coord[1][i] - pc[1];
        coord[2][i] = coord[2][i] - pc[2];
    }
}








/**
 * Helper function to trim whitespace from string (equivalent to Fortran's adjustl)
 */
auto trim(const std::string& str) -> std::string {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

/**
 * Check if pg1 is a subgroup of pg2
 */
auto issubgroup(const std::string& pg1, const std::string& pg2) -> bool {
    int npg;
    
    // Find pg2 in pgsymb array
    for (npg = 0; npg < 57; npg++) {
        if (trim(pg2) == trim(SymmetryData::pgsymb[npg])) break;
    }

    if (npg >= 57) return false; // pg2 not found

    // Check if pg1 is in the subgroup list for pg2
    for (int i = SymmetryData::nsgb[npg][0]; i < SymmetryData::nsgb[npg][0] + SymmetryData::nsgb[npg][1]; i++) {
        if (trim(SymmetryData::pgsymb[SymmetryData::nsgr[i] - 1]) == trim(pg1)) { // Adjust for 0-based indexing
            return true;
        }
    }
    
    return false;
}





// Implementation of symm_check function
void symm_check(int natoms, double delta, const std::vector<int>& nat,
                const std::vector<std::vector<double>>& coord,
                const std::vector<std::vector<double>>& cord,
                int& nc, std::vector<int>& ntrans, double& delta3) {
    nc = 0;
    delta3 = 0.0;

    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < natoms; ++j) {
            if (nat[i] != nat[j]) continue;

            std::array<double, 3> diff = {
                coord[0][i] - cord[0][j],
                coord[1][i] - cord[1][j],
                coord[2][i] - cord[2][j]
            };

            double vn = std::sqrt(symm_dot(diff.data(), diff.data(), 3));
            if (vn <= delta) {
                nc++;
                ntrans[i] = j;
                if (vn > delta3) delta3 = vn;
                break; // Found match, move to next i
            }
        }
    }
}

// Implementation of add_perm function
void add_perm(int natoms, const std::vector<int>& ntrans, int& nprm,
              std::vector<std::vector<int>>& nper) {
    // Check if this permutation already exists
    for (int i = 0; i < nprm; ++i) {
        bool match = true;
        for (int j = 0; j < natoms; ++j) {
            if (ntrans[j] != nper[j][i]) {
                match = false;
                break;
            }
        }
        if (match) return; // Already exists
    }

    // Add new permutation
    if (nprm >= 250) {
        throw std::runtime_error("ERROR: You need to enlarge the second dimension of 'nper' array and recompile the code");
    }

    for (int j = 0; j < natoms; ++j) {
        nper[j][nprm] = ntrans[j];
    }
    nprm++;
}

// Implementation of add_SG function
void add_SG(int& nsg, std::vector<std::array<double, 3>>& sigman,
            const std::array<double, 3>& v, double delta) {
    // Suppress unused parameter warning
    (void)v;
    // Check if this plane already exists
    for (int k = 0; k < nsg; ++k) {
        std::array<double, 3> existing = sigman[k];
        double vk = symm_dot(v.data(), existing.data(), 3);
        if (std::abs(vk) >= (1.0 - delta) && std::abs(vk) <= (1.0 + delta)) {
            return; // Already exists
        }
    }

    // Add new plane
    if (nsg >= 150) {
        throw std::runtime_error("ERROR: Too many symmetry operations. Try a lower tolerance.");
    }

    sigman[nsg][0] = v[0];
    sigman[nsg][1] = v[1];
    sigman[nsg][2] = v[2];
    nsg++;
}

// Implementation of add_Cn function
void add_Cn(int& nrot, std::vector<std::array<double, 3>>& rotn,
            std::vector<double>& rota, const std::array<double, 3>& v,
            const std::array<double, 3>& point, double alpha, double delta) {
    // Suppress unused parameter warning
    (void)point;

    // Check if this axis already exists
    for (int k = 0; k < nrot; ++k) {
        std::array<double, 3> existing = rotn[k];
        double vk = symm_dot(v.data(), existing.data(), 3);
        if (std::abs(vk) >= (1.0 - delta) && std::abs(vk) <= (1.0 + delta)) {
            if (rota[k] > alpha) rota[k] = alpha;
            return; // Already exists
        }
    }

    // Add new axis
    if (nrot >= 150) {
        throw std::runtime_error("ERROR: Too many symmetry operations. Try a lower tolerance.");
    }

    rotn[nrot][0] = v[0];
    rotn[nrot][1] = v[1];
    rotn[nrot][2] = v[2];
    rota[nrot] = alpha;
    nrot++;
}

/**
 * Scalar (dot or inner) product of two vectors: <x|y>=x.y
 */
auto symm_dot(const double* a, const double* b, int n) -> double {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += a[i] * b[i];
    }
    return result;
}



/**
 * Vector (cross or outer) product of two vectors: z=[x,y]
 */
void symm_crossp(const std::array<double, 3>& x, const std::array<double, 3>& y, std::array<double, 3>& z) {
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = -x[0] * y[2] + x[2] * y[0];
    z[2] = x[0] * y[1] - x[1] * y[0];
}

/**
 * Greatest common divisor
 */
auto symm_igcd(int a, int b) -> int {
    while (b != 0) {
        int t = b;
        b = a % b;
        a = t;
    }
    return a;
}

/**
 * Performs an inversion to the origin
 */
void symm_inversion(int natoms, const std::vector<int>& nat, 
                   const std::vector<std::vector<double>>& coord,
                   double delta, int& nc, std::vector<int>& ntrans, double delta3) {
    std::vector<std::vector<double>> cord(3, std::vector<double>(natoms));
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < natoms; j++) {
            cord[i][j] = -coord[i][j];
        }
    }
    
    symm_check(natoms, delta, nat, coord, cord, nc, ntrans, delta3);
    
    if (nc != natoms) {
        for (int i = 0; i < natoms; i++) {
            ntrans[i] = i + 1; // Fortran uses 1-based indexing
        }
    }
}

/**
 * Performs a rotation around axis v1
 */
void symm_rotate(int natoms, const std::vector<int>& nat,
                const std::vector<std::vector<double>>& coord,
                const std::array<double, 3>& v1, double sina, double cosa,
                double delta, int& nc, std::vector<int>& ntrans, double& delta3) {
    std::vector<std::vector<double>> cord(3, std::vector<double>(natoms));
    std::array<double, 3> v2, v3, p0, p;
    
    for (int j = 0; j < natoms; j++) {
        p0[0] = coord[0][j];
        p0[1] = coord[1][j];
        p0[2] = coord[2][j];
        
        symm_crossp(v1, p0, v2);
        symm_crossp(v1, v2, v3);
        
        p[0] = p0[0] + sina * v2[0] + (1.0 - cosa) * v3[0];
        p[1] = p0[1] + sina * v2[1] + (1.0 - cosa) * v3[1];
        p[2] = p0[2] + sina * v2[2] + (1.0 - cosa) * v3[2];
        
        cord[0][j] = p[0];
        cord[1][j] = p[1];
        cord[2][j] = p[2];
    }
    
    symm_check(natoms, delta, nat, coord, cord, nc, ntrans, delta3);
}

/**
 * Performs a reflection to the plane v
 */
void symm_reflect(int natoms, const std::vector<int>& nat,
                  const std::vector<std::vector<double>>& coord,
                  const std::array<double, 3>& v, const std::array<double, 3>& p0,
                  double delta, int& nc, std::vector<int>& ntrans, double& delta3) {
    std::vector<std::vector<double>> cord(3, std::vector<double>(natoms));
    std::array<double, 3> p, p_minus_p0;
    
    for (int i = 0; i < natoms; i++) {
        p[0] = coord[0][i];
        p[1] = coord[1][i];
        p[2] = coord[2][i];
        
        // Calculate p - p0
        p_minus_p0[0] = p[0] - p0[0];
        p_minus_p0[1] = p[1] - p0[1];
        p_minus_p0[2] = p[2] - p0[2];
        
        double vk = -symm_dot(v.data(), p_minus_p0.data(), 3);
        
        cord[0][i] = coord[0][i] + 2.0 * vk * v[0];
        cord[1][i] = coord[1][i] + 2.0 * vk * v[1];
        cord[2][i] = coord[2][i] + 2.0 * vk * v[2];
    }
    
    symm_check(natoms, delta, nat, coord, cord, nc, ntrans, delta3);
}

/**
 * Improper rotation around axis v1
 */
void symm_srotate(int natoms, const std::vector<int>& nat,
                  const std::vector<std::vector<double>>& coord,
                  const std::array<double, 3>& v1, double sina, double cosa,
                  double delta, int& nc, std::vector<int>& ntrans, double& delta3) {
    std::vector<std::vector<double>> cord(3, std::vector<double>(natoms));
    std::vector<std::vector<double>> cc(3, std::vector<double>(natoms));
    std::array<double, 3> v2, v3, p0, p;
    
    // First step: perform rotation
    for (int l = 0; l < natoms; l++) {
        p0[0] = coord[0][l];
        p0[1] = coord[1][l];
        p0[2] = coord[2][l];
        
        symm_crossp(v1, p0, v2);
        symm_crossp(v1, v2, v3);
        
        p[0] = p0[0] + sina * v2[0] + (1.0 - cosa) * v3[0];
        p[1] = p0[1] + sina * v2[1] + (1.0 - cosa) * v3[1];
        p[2] = p0[2] + sina * v2[2] + (1.0 - cosa) * v3[2];
        
        cord[0][l] = p[0];
        cord[1][l] = p[1];
        cord[2][l] = p[2];
    }
    
    // Second step: perform reflection through plane perpendicular to v1
    for (int j = 0; j < natoms; j++) {
        p0[0] = cord[0][j];
        p0[1] = cord[1][j];
        p0[2] = cord[2][j];
        
        double vk = -symm_dot(v1.data(), p0.data(), 3);
        
        cc[0][j] = cord[0][j] + 2.0 * vk * v1[0];
        cc[1][j] = cord[1][j] + 2.0 * vk * v1[1];
        cc[2][j] = cord[2][j] + 2.0 * vk * v1[2];
    }
    
    symm_check(natoms, delta, nat, coord, cc, nc, ntrans, delta3);
}



/**
 * @brief Determines the point group based on symmetry operation counts
 * @param ngp Total number of symmetry operations (group order)
 * @param ni Number of inversion operations (0 or 1)
 * @param nsg Number of mirror planes
 * @param ncr Number of proper rotations
 * @param nsr Number of improper rotations
 * @param np Order of principal rotation axis
 * @param PGlab Reference to string where point group symbol will be stored (e.g., "C1", "D2h")
 * @param nout Output control (0: silent, >=1: print diagnostics)
 */
auto symm_point_group(int ngp, int ni, int nsg, int ncr, int nsr, int np, int nout) -> std::string {
    // Compute principal axis order (ng[i][5], 0-based) for each point group
    std::array<int, SymmetryData::max_pgs> principal_axis{};
    for (int i = 0; i < SymmetryData::max_pgs; ++i) {
        const std::string& s = SymmetryData::sg[i];
        if (s == "Civ" || s == "Dih") {
            principal_axis[i] = -1;
        } else if (s[0] == 'C' || s[0] == 'D') {
            try {
                principal_axis[i] = s[1] >= '1' && s[1] <= '9' ? std::stoi(s.substr(1, 1)) : 0;
            } catch (const std::exception&) {
                principal_axis[i] = 0;
            }
        } else if (s[0] == 'S') {
            try {
                principal_axis[i] = std::stoi(s.substr(1, 1)) / 2;
            } catch (const std::exception&) {
                principal_axis[i] = 0;
            }
        } else if (s[0] == 'T') {
            principal_axis[i] = 3;
        } else if (s[0] == 'O') {
            principal_axis[i] = 4;
        } else if (s[0] == 'I') {
            principal_axis[i] = 5;
        } else {
            principal_axis[i] = 0;
        }
    }

    // Match point group
    std::string pgrp = "   ";
    bool yesno = false;
    for (int i = 0; i < SymmetryData::max_pgs; ++i) {
        if (ngp == SymmetryData::ng[i][0] &&
            ni == SymmetryData::ng[i][1] &&
            nsg == SymmetryData::ng[i][2] &&
            ncr == SymmetryData::ng[i][3] &&
            nsr == SymmetryData::ng[i][4] &&
            np == principal_axis[i]) {
            yesno = true;
            pgrp = SymmetryData::sg[i];
            if (nout >= 1) {
                std::cout << "\n-- POINT GROUP --\n";
                std::cout << "\n-- The structure should belong to the " << pgrp << " point group.\n";
                std::cout << "\n   " << pgrp << " = " << SymmetryData::cg[i] << "\n";
                std::cout << "\n       g       E       i      SG      Cn      Sn\n";
                std::cout << "   -----------------------------------------------\n";
                if (SymmetryData::ng[i][0] == -1) {
                    std::cout << "         inf" << std::setw(8) << 1 << std::setw(8) << SymmetryData::ng[i][1]
                              << std::setw(8) << "inf" << std::setw(8) << "inf";
                    if (SymmetryData::ng[i][4] == -1) {
                        std::cout << std::setw(8) << "inf\n";
                    } else {
                        std::cout << std::setw(8) << SymmetryData::ng[i][4] << "\n";
                    }
                } else {
                    std::cout << std::setw(10) << SymmetryData::ng[i][0] << std::setw(8) << 1;
                    for (int j = 1; j < 5; ++j) {
                        std::cout << std::setw(8) << SymmetryData::ng[i][j];
                    }
                    std::cout << "\n";
                }
                // Polarity and chirality
                if ((i >= 3 && i <= 16) || i == 48 || i == 51 || i == 53) {
                    if (i >= 3 && i <= 9) {
                        std::cout << "\n    The molecule is polar and chiral.\n";
                    } else {
                        std::cout << "\n    The molecule is not polar and chiral.\n";
                    }
                } else if (i == 1 || (i >= 17 && i <= 23) || i == 55) {
                    std::cout << "\n    The molecule is polar and not chiral.\n";
                } else {
                    std::cout << "\n    The molecule is not polar and not chiral.\n";
                }
            }
        }
    }

    // Handle special cases
    if (!yesno) {
        if (ngp == 1) {
            pgrp = "C1";
            if (nout >= 1) {
                std::cout << "\n-- POINT GROUP --\n";
                std::cout << "\n-- The structure should belong to the C1 point group.\n";
                std::cout << "\n    The molecule is polar and chiral.\n";
            }
        }
    }
    return pgrp;
}



/**
 * Get point group
 */
void PG_determ(int natoms, const std::vector<int>& nat, 
               std::vector<std::vector<double>>& coord, double delta, std::string& PGlab) {
    const int nmax = 150;  // maximum number of symmetric operations
    
    std::vector<double> pc(3);
    std::vector<std::vector<double>> symn(3, std::vector<double>(nmax));
    std::vector<std::vector<int>> nper(natoms, std::vector<int>(250));
    std::vector<std::vector<int>> nscl(natoms, std::vector<int>(natoms));
    std::vector<int> nccl(natoms);
    std::vector<std::vector<int>> nsym(nmax, std::vector<int>(5));
    
    if (natoms == 1) {
        PGlab = "C1";
        return;
    }
    
    int ncr = 0;
    int nsr = 0;
    int nsg = 0;
    int nout = 0;  // Suppress almost all output 
    int ng, ni, np, nprm, nseq;
    double wmol, cmx, cmy, cmz;
    
    // Calculation of the COM (centre of mass) of the molecule
    std::vector<double> wt_vec(data_common.wt.begin(), data_common.wt.end());
    symm_cmass(natoms, nat, wt_vec, coord, wmol, cmx, cmy, cmz);
    pc[0] = cmx;
    pc[1] = cmy;
    pc[2] = cmz;
    
    // Shift the origin of the Cartesian system to COM
    symm_cshift(natoms, coord, pc);
    
    // Convert symb array from data_common to vector for compatibility
    std::vector<std::string> symb_vec(90);
    for (int i = 0; i < 90; i++) {
        symb_vec[i] = data_common.symb[i];
    }
    
    // Find symmetry operations
    sym_elements(natoms, nat, coord, symb_vec, delta, ng, ni, nsg, ncr, nsr, np,
                symn, nsym, nout, nprm, nper, nseq, nccl, nscl);
    
    // Determine the equivalence classes defined by the symmetry operations
    symclass(natoms, nprm, nper, nseq, nccl, nscl, nat, symb_vec, nout);
    
    // Determine point group and framework group
    PGlab = symm_point_group(ng, ni, nsg, ncr, nsr, np, nout);
}

/**
 * Detect point group 
 * ishow=1: Output prompt =0: Do not print
 */
void SymmetryDetector::detectPG(int ishow) {
    std::vector<std::vector<double>> tmpmat(3, std::vector<double>(this->ncenter));
    std::string PGlabel3; // Used for input and output from PG_eqvatm

    if (this->PGlabelinit == "?") { // Not directly specified
        if (ishow == 1) std::cout << "Identifying point group..." << "\n";

        for (int i = 0; i < this->ncenter; i++) {
            tmpmat[0][i] = this->a[i].x;
            tmpmat[1][i] = this->a[i].y;
            tmpmat[2][i] = this->a[i].z;
        }

        std::vector<int> nat(this->ncenter);
        for (int i = 0; i < this->ncenter; i++) {
            nat[i] = this->a[i].index;
        }

        // This tolerance is suitable for most systems
        PG_determ(this->ncenter, nat, tmpmat, 0.01, PGlabel3);

        if (PGlabel3 == " " || PGlabel3.empty()) {
            for (int i = 1; i <= 20; i++) {
                PG_determ(this->ncenter, nat, tmpmat, i * 0.005, PGlabel3);
                if (PGlabel3 != " " && !PGlabel3.empty()) break;
            }
        }

        if (PGlabel3 == " " || PGlabel3.empty()) {
            std::cout << "Warning: Failed to identify point group; C1 will be used " << "\n";
            this->PGlabel = "C1  ";
        } else {
            if (ishow == 1) std::cout << "Point group has been successfully identified" << "\n";
            this->PGlabel = "    "; // Initialize with spaces
            this->PGlabel.replace(0, 3, PGlabel3);
        }
    } else {
        this->PGlabel = this->PGlabelinit;
    }

    this->PGlabel2rotsym();
}

/**
 * Convert point group label (PGlabel) to rotational symmetry number (rotsym)
 */
void SymmetryDetector::PGlabel2rotsym() {
    std::string PGlabel_trimmed = trim(this->PGlabel);
    int ie = PGlabel_trimmed.length();

    if (PGlabel_trimmed == "C1" || PGlabel_trimmed == "Ci" ||
        PGlabel_trimmed == "Cs" || PGlabel_trimmed == "Civ") {
        this->rotsym = 1;
    } else if (PGlabel_trimmed == "Dih") {
        this->rotsym = 2;
    } else if (PGlabel_trimmed == "Ih") {
        this->rotsym = 60;
    } else if (PGlabel_trimmed[0] == 'S') {
        std::string num_str = PGlabel_trimmed.substr(1);
        this->rotsym = std::stoi(num_str) / 2;
    } else if (PGlabel_trimmed == "T" || PGlabel_trimmed == "Td" || PGlabel_trimmed == "Th") {
        this->rotsym = 12; // Rotsym of Th is in line with link 717 of Gaussian
    } else if (PGlabel_trimmed == "Oh") {
        this->rotsym = 24;
    } else if (PGlabel_trimmed[0] == 'C') {
        try {
            this->rotsym = std::stoi(PGlabel_trimmed.substr(1));
        } catch (const std::exception&) {
            // The label is e.g. C2v, try without the last character
            this->rotsym = std::stoi(PGlabel_trimmed.substr(1, ie - 2));
        }
    } else if (PGlabel_trimmed[0] == 'D') {
        try {
            this->rotsym = std::stoi(PGlabel_trimmed.substr(1));
        } catch (const std::exception&) {
            // The label is e.g. D2h, try without the last character
            this->rotsym = std::stoi(PGlabel_trimmed.substr(1, ie - 2));
        }
        this->rotsym = this->rotsym * 2;
    } else {
        // Although We can identify 'O' and 'I', rotsym is not available
        // For simplicity, we'll assume rotsym = 1 for unknown cases
        std::cout << "Warning: Rotational symmetry number cannot be identified for this point group "
                  << PGlabel_trimmed << "\n";
        std::cout << "Assuming rotational symmetry number to be 1" << "\n";
        this->rotsym = 1;
    }
}



/** @brief 2. Subgroup boundaries (nsgb) */
/** @brief Static member initialization for nsgb_ array */
//const std::array<std::array<int, 2>, 57> nsgb_ = {{
//    {1,1}, {2,2}, {4,2}, {6,2}, {8,2},
//    {10,3}, {13,2}, {15,4}, {19,2}, {21,4},
//    {25,3}, {28,4}, {32,5}, {37,4}, {41,7},
//    {48,4}, {52,7}, {59,4}, {63,4}, {67,6},
//    {73,4}, {77,8}, {85,4}, {89,8}, {97,5},
//    {102,4}, {106,8}, {114,4}, {118,10}, {128,4},
//    {132,10}, {142,7}, {149,9}, {158,15}, {173,9},
//    {182,20}, {202,9}, {211,22}, {233,7}, {240,9},
//    {249,10}, {259,8}, {267,13}, {280,8}, {288,12},
//    {300,3}, {303,4}, {307,4}, {311,5}, {316,12},
//    {328,11}, {339,9}, {348,25}, {373,9}, {382,21},
//    {403,2}, {405,2}
//}};

// 2. Subgroup boundaries (nsgb) 
const std::array<int, 114> nsgb_1d = {
      1,   1,     2,   2,     4,   2,     6,   2,     8,   2,
     10,   3,    13,   2,    15,   4,    19,   2,    21,   4,
     25,   3,    28,   4,    32,   5,    37,   4,    41,   7,
     48,   4,    52,   7,    59,   4,    63,   4,    67,   6,
     73,   4,    77,   8,    85,   4,    89,   8,    97,   5,
    102,   4,   106,   8,   114,   4,   118,  10,   128,   4,
    132,  10,   142,   7,   149,   9,   158,  15,   173,   9,
    182,  20,   202,   9,   211,  22,   233,   7,   240,   9,
    249,  10,   259,   8,   267,  13,   280,   8,   288,  12,
    300,   3,   303,   4,   307,   4,   311,   5,   316,  12,
    328,  11,   339,   9,   348,  25,   373,   9,   382,  21,
    403,   2,   405,   2
};

const std::array<std::array<int, 2>, 57> nsgb_ = convert_fortran_to_cpp<57, 2>(nsgb_1d);





// 5. Irreducible representation symbols (irsymb)
const std::array<std::string, 322> irsymb_ = {{
          "A",   "A'",  "A\"",   "Ag",   "Au",    "A",     "B",    "A",
          "E",   " A",    "B",    "E",    "A",   "E1",    "E2",    "A",
          "B",   "E1",   "E2",    "A",   "E1",   "E2",    "E3",    "A",
          "B",   "E1",   "E2",   "E3",   " A",   "B1",    "B2",   "B3",
         "A1",   "A2",   " E",   "A1",   "A2",   "B1",    "B2",    "E",
         "A1",   "A2",   "E1",   "E2",   "A1",   "A2",    "B1",   "B2",
         "E1",   "E2",   "A1",   "A2",   "E1",   "E2",    "E3",   "A1",
         "A2",   "B1",   "B2",   "E1",   "E2",   "E3",    "A1",   "A2",
         "B1",   "B2",   "A1",   "A2",   " E",   "A1",    "A2",   "B1",
         "B2",    "E",   "A1",   "A2",   "E1",   "E2",    "A1",   "A2",
         "B1",   "B2",   "E1",   "E2",   "A1",   "A2",    "E1",   "E2",
         "E3",   "A1",   "A2",   "B1",   "B2",   "E1",    "E2",   "E3",
         "Ag",   "Bg",   "Au",   "Bu",   "A'",  "A\"",    "E'",  "E\"",
         "Ag",   "Bg",   "Eg",   "Au",   "Bu",   "Eu",    "A'",  "A\"",
        "E1'", "E1\"",  "E2'", "E2\"",   "Ag",   "Bg",   "E1g",  "E2g",
         "Au",   "Bu",  "E1u",  "E2u",   "A'",  "A\"",   "E1'", "E1\"",
        "E2'", "E2\"",  "E3'", "E3\"",   "Ag",   "Bg",   "E1g",  "E2g",
        "E3g",   "Au",   "Bu",  "E1u",  "E2u",  "E3u",    "Ag",  "B1g",
        "B2g",  "B3g",   "Au",  "B1u",  "B2u",  "B3u",   "A1'", "A1\"",
        "A2'", "A2\"",   "E'",  "E\"",  "A1g",  "A2g",   "B1g",  "B2g",
         "Eg",  "A1u",  "A2u",  "B1u",  "B2u",   "Eu",   "A1'", "A1\"",
        "A2'", "A2\"",  "E1'", "E1\"",  "E2'", "E2\"",   "A1g",  "A2g",
        "B1g",  "B2g",  "E1g",  "E2g",  "A1u",  "A2u",   "B1u",  "B2u",
        "E1u",  "E2u",  "A1'", "A1\"",  "A2'", "A2\"",   "E1'", "E1\"",
        "E2'", "E2\"",  "E3'", "E3\"",  "A1g",  "A2g",   "B1g",  "B2g",
        "E1g",  "E2g",  "E3g",  "A1u",  "A2u",  "B1u",   "B2u",  "E1u",
        "E2u",  "E3u",   "A1",   "A2",   "B1",   "B2",     "E",  "A1g",
        "A2g",   "Eg",  "A1u",  "A2u",   "Eu",   "A1",    "A2",   "B1",
        " B2",   "E1",   "E2",   "E3",  "A1g",  "A2g",   "E1g",  "E2g",
        "A1u",  "A2u",  "E1u",  "E2u",   "A1",   "A2",    "B1",   "B2",
         "E1",   "E2",   "E3",   "E4",   "E5",  "A1g",   "A2g",  "E1g",
        "E2g",  "E3g",  "A1u",  "A2u",  "E1u",  "E2u",   "E3u",   "A1",
         "A2",   "B1",  "B2",    "E1",   "E2",   "E3",    "E4",   "E5",
         "E6",   "E7",   "A",     "B",    "E",   "Ag",    "Eg",   "Au",
         "Eu",    "A",   "B",    "E1",   "E2",   "E3",    " A",    "E",
          "T",   "Ag",  "Eg",    "Tg",   "Au",   "Eu",    "Tu",   "A1",
         "A2",    "E",  "T1",    "T2",   "A1",   "A2",     "E",   "T1",
         "T2",  "A1g", "A2g",    "Eg",  "T1g",  "T2g",   "A1u",   "2u",
         "Eu",  "T1u", "T2u",     "A",   "T1",   "T2",     "G",    "H",
         "Ag",  "T1g", "T2g",    "Gg",   "Hg",   "Au",   "T1u",  "T2u",
         "Gu",  "Hu"
}};




// 6. Rotational harmonics (nrotharm)
const std::array<int, 966> nrotharm_1d = {
    3, 3, 5,   1, 2, 3,   2, 1, 2,   3, 0, 5,   0, 3, 0,   1, 1, 3,
    2, 2, 2,   1, 1, 1,   2, 2, 4,   1, 1, 1,   0, 0, 2,   2, 2, 2,
    1, 1, 1,   2, 2, 2,   0, 0, 2,   1, 1, 1,   0, 0, 0,   2, 2, 2,
    0, 0, 2,   1, 1, 1,   2, 2, 2,   0, 0, 2,   0, 0, 0,   1, 1, 1,
    0, 0, 0,   2, 2, 2,   0, 0, 2,   0, 0, 0,   0, 0, 2,   1, 1, 1,
    1, 1, 1,   1, 1, 1,   0, 0, 1,   1, 1, 0,   2, 2, 4,   0, 0, 1,
    1, 1, 0,   0, 0, 1,   0, 0, 1,   2, 2, 2,   0, 0, 1,   1, 1, 0,
    2, 2, 2,   0, 0, 2,   0, 0, 1,   1, 1, 0,   0, 0, 0,   0, 0, 0,
    2, 2, 2,   0, 0, 2,   0, 0, 1,   1, 1, 0,   2, 2, 2,   0, 0, 2,
    0, 0, 0,   0, 0, 1,   1, 1, 0,   0, 0, 0,   0, 0, 0,   2, 2, 2,
    0, 0, 2,   0, 0, 0,   0, 1, 2,   1, 0, 1,   1, 1, 1,   1, 1, 1,
    0, 1, 1,   1, 0, 0,   2, 2, 4,   0, 1, 1,   1, 0, 0,   0, 0, 1,
    0, 0, 1,   2, 2, 2,   0, 1, 1,   1, 0, 0,   2, 2, 2,   0, 0, 2,
    0, 1, 1,   1, 0, 0,   0, 0, 0,   0, 0, 0,   2, 2, 2,   0, 0, 2,
    0, 1, 1,   1, 0, 0,   2, 2, 2,   0, 0, 2,   0, 0, 0,   0, 1, 1,
    1, 0, 0,   0, 0, 0,   0, 0, 0,   2, 2, 2,   0, 0, 2,   0, 0, 0,
    1, 0, 3,   2, 0, 2,   0, 1, 0,   0, 2, 0,   1, 0, 1,   0, 1, 0,
    0, 2, 2,   2, 0, 2,   1, 0, 1,   0, 0, 2,   2, 0, 2,   0, 1, 0,
    0, 0, 0,   0, 2, 0,   1, 0, 1,   0, 1, 0,   0, 2, 0,   2, 0, 2,
    0, 0, 2,   0, 0, 0,   1, 0, 1,   0, 0, 0,   2, 0, 2,   0, 0, 2,
    0, 1, 0,   0, 0, 0,   0, 2, 0,   0, 0, 0,   1, 0, 1,   0, 1, 0,
    0, 2, 0,   2, 0, 2,   0, 0, 2,   0, 0, 0,   0, 0, 0,   0, 0, 0,
    1, 0, 1,   0, 0, 0,   2, 0, 2,   0, 0, 2,   0, 0, 0,   0, 1, 0,
    0, 0, 0,   0, 2, 0,   0, 0, 0,   0, 0, 0,   0, 0, 2,   1, 0, 1,
    1, 0, 1,   1, 0, 1,   0, 0, 0,   0, 1, 0,   0, 1, 0,   0, 1, 0,
    0, 0, 1,   0, 0, 0,   1, 0, 0,   0, 1, 0,   0, 2, 2,   2, 0, 2,
    0, 0, 1,   1, 0, 0,   0, 0, 1,   0, 0, 1,   2, 0, 2,   0, 0, 0,
    0, 1, 0,   0, 0, 0,   0, 0, 0,   0, 2, 0,   0, 0, 1,   0, 0, 0,
    1, 0, 0,   0, 1, 0,   0, 2, 0,   2, 0, 2,   0, 0, 2,   0, 0, 0,
    0, 0, 1,   1, 0, 0,   0, 0, 0,   0, 0, 0,   2, 0, 2,   0, 0, 2,
    0, 0, 0,   0, 1, 0,   0, 0, 0,   0, 0, 0,   0, 2, 0,   0, 0, 0,
    0, 0, 1,   0, 0, 0,   1, 0, 0,   0, 1, 0,   0, 2, 0,   2, 0, 2,
    0, 0, 2,   0, 0, 0,   0, 0, 0,   0, 0, 0,   0, 0, 1,   1, 0, 0,
    0, 0, 0,   0, 0, 0,   2, 0, 2,   0, 0, 2,   0, 0, 0,   0, 0, 0,
    0, 1, 0,   0, 0, 0,   0, 0, 0,   0, 2, 0,   0, 0, 0,   0, 0, 0,
    0, 0, 1,   1, 0, 0,   0, 0, 1,   0, 1, 1,   2, 2, 2,   0, 0, 1,
    1, 0, 0,   2, 0, 4,   0, 0, 0,   0, 1, 0,   0, 2, 0,   0, 0, 1,
    1, 0, 0,   0, 0, 0,   0, 1, 0,   0, 2, 0,   0, 0, 2,   2, 0, 2,
    0, 0, 1,   1, 0, 0,   2, 0, 2,   0, 0, 2,   0, 0, 0,   0, 1, 0,
    0, 2, 0,   0, 0, 0,   0, 0, 1,   1, 0, 0,   0, 0, 0,   0, 1, 0,
    0, 2, 0,   0, 0, 2,   0, 0, 0,   0, 0, 0,   2, 0, 2,   0, 0, 1,
    1, 0, 0,   2, 0, 2,   0, 0, 2,   0, 0, 0,   0, 0, 0,   0, 1, 0,
    0, 2, 0,   0, 0, 0,   0, 0, 0,   0, 0, 1,   1, 0, 0,   0, 0, 0,
    0, 1, 0,   0, 2, 0,   0, 0, 2,   0, 0, 0,   0, 0, 0,   0, 0, 0,
    0, 0, 0,   2, 0, 2,   1, 0, 1,   0, 1, 2,   2, 2, 2,   1, 0, 1,
    2, 0, 4,   0, 1, 0,   0, 2, 0,   1, 0, 1,   0, 1, 0,   0, 2, 0,
    0, 0, 2,   2, 0, 2,   0, 0, 0,   0, 0, 2,   3, 3, 3,   0, 0, 0,
    0, 0, 2,   3, 0, 3,   0, 0, 0,   0, 0, 0,   0, 3, 0,   0, 0, 0,
    0, 0, 0,   0, 0, 2,   3, 0, 0,   0, 3, 3,   0, 0, 0,   0, 0, 0,
    0, 0, 2,   3, 3, 0,   0, 0, 3,   0, 0, 0,   0, 0, 0,   0, 0, 2,
    3, 0, 0,   0, 0, 3,   0, 0, 0,   0, 0, 0,   0, 0, 0,   0, 3, 0,
    0, 0, 0,   0, 0, 0,   3, 3, 0,   0, 0, 0,   0, 0, 0,   0, 0, 5,
    0, 0, 0,   3, 0, 0,   0, 0, 0,   0, 0, 0,   0, 0, 5,   0, 0, 0,
    0, 3, 0,   0, 0, 0,   0, 0, 0,   0, 0, 0
};

// Function to initialize nrotharm_ at compile time
auto initialize_nrotharm() -> std::array<std::array<int, 322>, 3> {
    return convert_fortran_to_cpp<3, 322>(nrotharm_1d);
}

// Initialize the static member
std::array<std::array<int, 322>, 3> nrotharm_ = initialize_nrotharm();


// Define the static member
const std::array<std::array<std::array<int, 55>, 4>, 14> nsymop_ = initialize_nsymop();


// Function to map point group and irrep indices to linear irsymb index
auto getIrrepSymbolIndex(int point_group_idx, int irrep_idx) -> int {
    // The mapping uses the nir_ array to find correct starting column for each point group

    if (point_group_idx < 0 || point_group_idx >= 55 || irrep_idx < 0) {
        return -1;
    }

    // Get the starting column and count for this point group from nir_ array
    int start_col = nir_[point_group_idx][0] - 1; // 0-based
    int irrep_count = nir_[point_group_idx][1];

    // Check if irrep_idx is valid for this point group
    if (irrep_idx >= irrep_count) {
        return -1;
    }

    // Calculate the actual index in the irsymb_ array
    int index = start_col + irrep_idx;

    // Bounds check
    if (index < 0 || index >= 322) {
        return -1;
    }

    return index;
}

// Accessor functions 
auto getIrrepSymbol(int point_group_idx, int irrep_idx) -> std::string {
    int index = getIrrepSymbolIndex(point_group_idx, irrep_idx);
    if (index < 0 || index >= 322) {
        return "";
    }

    return irsymb_[index];
}

auto getRotationalHarmonic(int point_group_idx, int irrep_idx) -> int {
    // Map to the correct irrep index (similar to irsymb mapping)
    int index = getIrrepSymbolIndex(point_group_idx, irrep_idx);
    if (index < 0 || index >= 322) {
        return 0;
    }

    /** @brief Return the first harmonic value (nrotharm(1,index) in Fortran) */
    return nrotharm_[0][index];
}

/**
 * @brief Get specific rotational harmonic value by type
 *
 * Returns a specific type of rotational harmonic (1, 2, or 3 in Fortran terms,
 * 0, 1, or 2 in C++) for advanced spectroscopic analysis.
 *
 * @param point_group_idx Point group index (0-56)
 * @param irrep_idx Irreducible representation index within the point group
 * @param harmonic_type Harmonic type (0, 1, or 2)
 * @return Specific harmonic value
 */
auto getRotationalHarmonic(int point_group_idx, int irrep_idx, int harmonic_type) -> int {
    int index = getIrrepSymbolIndex(point_group_idx, irrep_idx);
    if (index < 0 || index >= 322 || harmonic_type < 0 || harmonic_type >= 3) {
        return 0;
    }

    return nrotharm_[harmonic_type][index];
}




// Static point group data
const std::array<std::string, SymmetryData::max_pgs> SymmetryData::sg = {
    "C1 ", "Cs ", "Ci ", "C2 ", "C3 ", "C4 ", "C5 ", "C6 ", "C7 ", "C8 ",
    "D2 ", "D3 ", "D4 ", "D5 ", "D6 ", "D7 ", "D8 ", "C2v", "C3v", "C4v",
    "C5v", "C6v", "C7v", "C8v", "C2h", "C3h", "C4h", "C5h", "C6h", "C7h",
    "C8h", "D2h", "D3h", "D4h", "D5h", "D6h", "D7h", "D8h", "D2d", "D3d",
    "D4d", "D5d", "D6d", "D7d", "D8d", "S4 ", "S6 ", "S8 ", "T  ", "Th ",
    "Td ", "O  ", "Oh ", "I  ", "Ih ", "Civ", "Dih"
};

const std::array<std::array<int, 6>, SymmetryData::max_pgs> SymmetryData::ng = {{
    {  1, 0,  0,  0,  0,  0}, // C1
    {  2, 0,  1,  0,  0,  0}, // Cs
    {  2, 1,  0,  0,  0,  0}, // Ci
    {  2, 0,  0,  1,  0,  0}, // C2
    {  3, 0,  0,  2,  0,  0}, // C3
    {  4, 0,  0,  3,  0,  0}, // C4
    {  5, 0,  0,  4,  0,  0}, // C5
    {  6, 0,  0,  5,  0,  0}, // C6
    {  7, 0,  0,  6,  0,  0}, // C7
    {  8, 0,  0,  7,  0,  0}, // C8
    {  4, 0,  0,  3,  0,  0}, // D2
    {  6, 0,  0,  5,  0,  0}, // D3
    {  8, 0,  0,  7,  0,  0}, // D4
    { 10, 0,  0,  9,  0,  0}, // D5
    { 12, 0,  0, 11,  0,  0}, // D6
    { 14, 0,  0, 13,  0,  0}, // D7
    { 16, 0,  0, 15,  0,  0}, // D8
    {  4, 0,  2,  1,  0,  0}, // C2v
    {  6, 0,  3,  2,  0,  0}, // C3v
    {  8, 0,  4,  3,  0,  0}, // C4v
    { 10, 0,  5,  4,  0,  0}, // C5v
    { 12, 0,  6,  5,  0,  0}, // C6v
    { 14, 0,  7,  6,  0,  0}, // C7v
    { 16, 0,  8,  7,  0,  0}, // C8v
    {  4, 1,  1,  1,  0,  0}, // C2h
    {  6, 0,  1,  2,  2,  0}, // C3h
    {  8, 1,  1,  3,  2,  0}, // C4h
    { 10, 0,  1,  4,  4,  0}, // C5h
    { 12, 1,  1,  5,  4,  0}, // C6h
    { 14, 0,  1,  6,  6,  0}, // C7h
    { 16, 1,  1,  7,  6,  0}, // C8h
    {  8, 1,  3,  3,  0,  0}, // D2h
    { 12, 0,  4,  5,  2,  0}, // D3h
    { 16, 1,  5,  7,  2,  0}, // D4h
    { 20, 0,  6,  9,  4,  0}, // D5h
    { 24, 1,  7, 11,  4,  0}, // D6h
    { 28, 0,  8, 13,  6,  0}, // D7h
    { 32, 1,  9, 15,  6,  0}, // D8h
    {  8, 0,  2,  3,  2,  0}, // D2d
    { 12, 1,  3,  5,  2,  0}, // D3d
    { 16, 0,  4,  7,  4,  0}, // D4d
    { 20, 1,  5,  9,  4,  0}, // D5d
    { 24, 0,  6, 11,  6,  0}, // D6d
    { 28, 1,  7, 13,  6,  0}, // D7d
    { 32, 0,  8, 15,  8,  0}, // D8d
    {  4, 0,  0,  1,  2,  0}, // S4
    {  6, 1,  0,  2,  2,  0}, // S6
    {  8, 0,  0,  3,  4,  0}, // S8
    { 12, 0,  0, 11,  0,  0}, // T
    { 24, 1,  3, 11,  8,  0}, // Th
    { 24, 0,  6, 11,  6,  0}, // Td
    { 24, 0,  0, 23,  0,  0}, // O
    { 48, 1,  9, 23, 14,  0}, // Oh
    { 60, 0,  0, 59,  0,  0}, // I
    {120, 1, 15, 59, 44,  0}, // Ih
    { -1, 0,  1,  1,  0, -1}, // Civ
    { -1, 1,  1,  2,  1, -1}  // Dih
}};

const std::array<std::string, SymmetryData::max_pgs> SymmetryData::cg = {
    "{(E)}                                                  ",
    "{(E) (SG)}                                             ",
    "{(E) (i)}                                              ",
    "{(E) (C2)}                                             ",
    "{(E) (C3)}                                             ",
    "{(E) (C4) (C2)}                                        ",
    "{(E) (C5)}                                             ",
    "{(E) (C6) (C3) (C2)}                                   ",
    "{(E) (C7)}                                             ",
    "{(E) (C8) (C4) (C2)}                                   ",
    "{(E) 3*(C2)}                                          ",
    "{(E) (C3) 3*(C2)}                                      ",
    "{(E) (C4) 5*(C2)}                                      ",
    "{(E) (C5) 5*(C2)}                                      ",
    "{(E) (C6) (C3) 7*(C2)}                                 ",
    "{(E) (C7) 7*(C2)}                                      ",
    "{(E) (C8) (C4) 9*(C2)}                                 ",
    "{(E) (C2) 2*(SG)}                                      ",
    "{(E) (C3) 3*(SG)}                                      ",
    "{(E) (C4) (C2) 4*(SG)}                                 ",
    "{(E) (C5) 5*(SG)}                                      ",
    "{(E) (C6) (C3) (C2) 6*(SG)}                            ",
    "{(E) (C7) 7*(SG)}                                      ",
    "{(E) (C8) (C4) (C2) 8*(SG)}                            ",
    "{(E) (C2) (i) (SG)}                                    ",
    "{(E) (C3) (S3) (SG)}                                   ",
    "{(E) (i) (C4) (C2) (S4) (SG)}                          ",
    "{(E) (C5) (S5) (SG)}                                   ",
    "{(E) (i) (C6) (C3) (C2) (S6) (S3) (SG)}                ",
    "{(E) (C7) (S7) (SG)}                                   ",
    "{(E) (i) (C8) (C4) (C2) (S8) (S4) (SG)}                ",
    "{(E) (i) 3*(C2) 3*(SG)}                                ",
    "{(E) (C3) 3*(C2) (S3) 4*(SG)}                          ",
    "{(E) (i) (C4) 5*(C2) (S4) 5*(SG)}                      ",
    "{(E) (C5) 5*(C2) (S5) 6*(SG)}                          ",
    "{(E) (i) (C6) (C3) 7*(C2) (S6) (S3) 7*(SG)}            ",
    "{(E) (C7) 7*(C2) (S7) 8*(SG)}                          ",
    "{(E) (i) (C8) (C4) 9*(C2) (S8) (S4) 9*(SG)}            ",
    "{(E) 3*(C2) (S4) 2*(SG)}                               ",
    "{(E) (i) (C3) 3*(C2) (S6) 3*(SG)}                      ",
    "{(E) (C4) 5*(C2) (S8) 4*(SG)}                          ",
    "{(E) (i) (C5) 5*(C2) (S10) 5*(SG)}                     ",
    "{(E) (C6) (C3) 7*(C2) (S12) (S4) 6*(SG)}               ",
    "{(E) (i) (C7) 7*(C2) (S14) 7*(SG)}                     ",
    "{(E) (C8) (C4) 9*(C2) (S16) 8*(SG)}                    ",
    "{(E) (C2) (S4)}                                        ",
    "{(E) (i) (C3) (S6)}                                    ",
    "{(E) (C4) (C2) (S8)}                                   ",
    "{(E) 4*(C3) 3*(C2)}                                    ",
    "{(E) (i) 4*(C3) 3*(C2) 4*(S6) 3*(SG)}                  ",
    "{(E) 4*(C3) 3*(C2) 3*(S4) 6*(SG)}                      ",
    "{(E) 3*(C4) 4*(C3) 9*(C2)}                             ",
    "{(E) (i) 3*(C4) 4*(C3) 9*(C2) 4*(S6) 3*(S4) 9*(SG)}    ",
    "{(E) 6*(C5) 10*(C3) 15*(C2)}                           ",
    "{(E) (i) 6*(C5) 10*(C3) 15*(C2) 6*(S10) 10*(S6) 15*(SG)}",
    "{(E) 2*(Cinf) ... inf*(SG)}                             ",
    "{(E) (i) 2*(Cinf) ... inf*(C2) 2*(Sinf) ... inf*(SG)}   "
};

// Initialize static members 
std::array<std::string, SymmetryData::max_pgs> SymmetryData::pgsymb = {
    "C1 ", "Cs ", "Ci ", "C2 ", "C3 ", "C4 ", "C5 ", "C6 ", "C7 ", "C8 ",
    "D2 ", "D3 ", "D4 ", "D5 ", "D6 ", "D7 ", "D8 ", "C2v", "C3v", "C4v",
    "C5v", "C6v", "C7v", "C8v", "C2h", "C3h", "C4h", "C5h", "C6h", "C7h",
    "C8h", "D2h", "D3h", "D4h", "D5h", "D6h", "D7h", "D8h", "D2d", "D3d",
    "D4d", "D5d", "D6d", "D7d", "D8d", "S4 ", "S6 ", "S8 ", "T  ", "Th ",
    "Td ", "O  ", "Oh ", "I  ", "Ih ", "Civ", "Dih"
};

std::array<std::array<int, 2>, SymmetryData::max_pgs> SymmetryData::nsgb = initialize_nsgb();

std::array<int, SymmetryData::max_subgroups> SymmetryData::nsgr = {{
    1,  1,  2,  1,  3,  1,  4,  1,  5,  1,
    4,  6,  1,  7,  1,  4,  5,  8,  1,  9,
    1,  4,  6, 10,  1,  4, 11,  1,  4,  5,
    12, 1,  4,  6, 11, 13,  1,  4,  7, 14,
    1,  4,  5,  8, 11, 12, 15,  1,  4,  9,
    16,  1,  4,  6, 10, 11, 13, 17,  1,  2,
    4, 18,  1,  2,  5, 19,  1,  2,  4,  6,
    18, 20,  1,  2,  7, 21,  1,  2,  4,  5,
    8, 18, 19, 22,  1,  2,  9, 23,  1,  2,
    4,  6, 10, 18, 20, 24,  1,  2,  3,  4,
    25,  1,  2,  5, 26,  1,  2,  3,  4,  6,
    25, 46, 27,  1,  2,  7, 28,  1,  2,  3,
    4,  5,  8, 25, 26, 47, 29,  1,  2,  9,
    30,  1,  2,  3,  4,  6, 10, 25, 27, 48,
    31,  1,  2,  3,  4, 18, 25, 32,  1,  2,
    4,  5, 12, 18, 19, 26, 33,  1,  2,  3,
    4,  6, 11, 13, 18, 20, 25, 27, 32, 39,
    46, 34,  1,  2,  4,  7, 14, 18, 21, 28,
    35,  1,  2,  3,  4,  5,  8, 11, 12, 15,
    18, 19, 22, 25, 26, 29, 32, 33, 40, 47,
    36,  1,  2,  4,  9, 16, 18, 23, 30, 37,
    1,  2,  3,  4,  6, 10, 11, 13, 17, 18,
    20, 24, 25, 27, 31, 32, 34, 39, 41, 46,
    48, 38,  1,  2,  4, 11, 18, 46, 39,  1,
    2,  3,  4,  5, 12, 19, 47, 40,  1,  2,
    4,  6, 11, 13, 18, 20, 48, 41,  1,  2,
    3,  4,  7, 14, 21, 42,  1,  2,  4,  5,
    8, 11, 12, 15, 18, 19, 22, 46, 43,  1,
    2,  3,  4,  9, 16, 23, 44,  1,  2,  4,
    6, 10, 11, 13, 17, 18, 20, 24, 45,  1,
    4, 46,  1,  3,  5, 47,  1,  4,  6, 48,
    1,  4,  5, 11, 49,  1,  2,  3,  4,  5,
    11, 18, 25, 32, 47, 49, 50,  1,  2,  4,
    5, 11, 18, 19, 39, 46, 49, 51,  1,  4,
    5,  6, 11, 12, 13, 49, 52,  1,  2,  3,
    4,  5,  6, 11, 12, 13, 18, 19, 20, 25,
    27, 32, 34, 39, 40, 46, 47, 49, 50, 51,
    52, 53,  1,  4,  5,  7, 11, 12, 14, 49,
    54,  1,  2,  3,  4,  5,  7, 11, 12, 14,
    18, 19, 21, 25, 32, 40, 42, 47, 49, 50,
    54, 55,  1, 56,  1, 57
}};








/**
 * @brief Function to initialize nsgb_ from 1D Fortran data
 * @return 2D array of subgroup boundaries
 */
auto initialize_nsgb() -> std::array<std::array<int, 2>, 57> {
    return convert_fortran_to_cpp<57, 2>(nsgb_1d);
}









// 4. Symmetry operations (nsymop) 


// Helper function to convert 1D Fortran column-major array to 3D C++ row-major array
template<size_t Dim1, size_t Dim2, size_t Dim3>
auto convert_fortran_to_cpp_3d(const std::vector<int>& fortran_1d, size_t dim1, size_t dim2, size_t dim3) 
    -> std::array<std::array<std::array<int, Dim3>, Dim2>, Dim1> {
    // Validate input size
    if (fortran_1d.size() != dim1 * dim2 * dim3) {
        throw std::invalid_argument("Error: Input array size does not match dimensions.");
    }

    std::array<std::array<std::array<int, Dim3>, Dim2>, Dim1> result{};



    // Map column-major (Fortran) to row-major (C++)
    // Fortran index: nsymop(i,j,k) at (k-1)*dim1*dim2 + (j-1)*dim1 + (i-1)
    // C++ index: result[i-1][j-1][k-1]
    for (size_t k = 0; k < dim3; ++k) {
        for (size_t j = 0; j < dim2; ++j) {
            for (size_t i = 0; i < dim1; ++i) {
                size_t index = k * dim1 * dim2 + j * dim1 + i;
                result[i][j][k] = fortran_1d[index];
            }
        }
    }

    return result;
}

// Function to initialize nsymop_ from 1D Fortran data
auto initialize_nsymop() -> std::array<std::array<std::array<int, 55>, 4>, 14> {
    // Fortran 1D initializer (column-major, 14*4*55 = 3080 elements)
    // Expanded from the Fortran data statement, with n*0 replaced by explicit zeros
    const std::vector<int> fortran_1d = {
        // c        E
        4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  SGH
        4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E   i
        4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C2
        4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C3
        4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C4  C2
        4, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E   C5   C5^2
        4, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E   C6   C3  C2
        4, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 6, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E   C7  C7^2  C7^3
        4, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E   C8  C4  C8^3  C2
        4, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 8, 4, 8, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C2   C2'  C2"
        4, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C3  C2'
        4, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C4  C2  C2'  C2"
        4, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 4, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C5  C5^2  C2'
        4, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 5, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C6  C3  C2  C2'  C2"
        4, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 6, 3, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 1, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C7 C7^2  C7^3  C2'
        4, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 7, 7, 7, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C8  C4   C8^3  C2  C2'  C2"
        4, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0,
        0, 8, 4, 8, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 3, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 1, 4, 4, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C2  SGV  SHD
        4, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C3  SGV
        4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C4  C2  SGV  SGD
        4, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 4, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C5  C5^2  SGV
        4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 5, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C6  C3  C2  SGV  SGD
        4, 2, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 6, 3, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 1, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E   C7  C7^2  C7^3  SGV
        4, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 7, 7, 7, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E   C8  C4  C8^3  C2  SGV  SGD
        4, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 8, 4, 8, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 1, 4, 4, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C2  i  SGH
        4, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C3  SGH  S3
        4, 2, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C4  C2  i   S4 SGH
        4, 2, 2, 0, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 4, 2, 1, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C5  C5^2  SGH  S5  S5^3
        4, 2, 2, 1, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 5, 5, 1, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C6  C3  C2  i   S6  S3 SGH
        4, 2, 2, 2, 0, 3, 3, 1, 0, 0, 0, 0, 0, 0,
        0, 6, 3, 2, 1, 6, 3, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0,
        // c        E   C7  C7^2  C7^3  SGH   S7  S7^3  S7^5
        4, 2, 2, 2, 1, 3, 3, 3, 0, 0, 0, 0, 0, 0,
        0, 7, 7, 7, 1, 7, 7, 7, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 3, 0, 1, 3, 5, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0,
        // c        E   C8  C4  C8^3  C2  i    S8  S4  S8^3  SGH
        4, 2, 2, 2, 2, 0, 3, 3, 3, 1, 0, 0, 0, 0,
        0, 8, 4, 8, 2, 1, 8, 4, 8, 1, 0, 0, 0, 0,
        0, 1, 1, 3, 0, 0, 1, 1, 3, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 1, 1, 2, 2, 2, 1, 0, 0, 0, 0,
        // c        E  C2  C2'  C2"  i  SGH  SGV  SGD
        4, 2, 2, 2, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 2, 2, 2, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        // c        E  C3  C2'  SGH  S3  SGV
        4, 2, 2, 1, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 2, 1, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 3, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C4  C2  C2'  C2"   i   S4 SGH  SGV  SGD
        4, 2, 2, 2, 2, 0, 3, 1, 1, 1, 0, 0, 0, 0,
        0, 4, 2, 2, 2, 1, 4, 1, 2, 3, 0, 0, 0, 0,
        0, 1, 0, 2, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 1, 2, 2, 1, 2, 1, 2, 2, 0, 0, 0, 0,
        // c        E  C5  C5^2  C2'  SGH   S5  S5^3 SGV
        4, 2, 2, 2, 1, 3, 3, 1, 0, 0, 0, 0, 0, 0,
        0, 5, 5, 2, 1, 5, 5, 2, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 2, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 5, 1, 2, 2, 5, 0, 0, 0, 0, 0, 0,
        // c        E  C6  C3  C2  C2'  C2"   i   S6  S3 SGH  SGV  SGD
        4, 2, 2, 2, 2, 2, 0, 3, 3, 1, 1, 1, 0, 0,
        0, 6, 3, 2, 2, 2, 1, 6, 3, 1, 2, 3, 0, 0,
        0, 1, 1, 0, 2, 3, 0, 1, 1, 0, 0, 0, 0, 0,
        1, 2, 2, 1, 3, 3, 1, 2, 2, 1, 3, 3, 0, 0,
        // c        E  C7 C7^2  C7^3  C2'  SGH   S7  S7^3  S7^5  SGV
        4, 2, 2, 2, 2, 1, 3, 3, 3, 1, 0, 0, 0, 0,
        0, 7, 7, 7, 2, 1, 7, 7, 7, 2, 0, 0, 0, 0,
        0, 1, 2, 3, 2, 0, 1, 3, 5, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 7, 1, 2, 2, 2, 7, 0, 0, 0, 0,
        // c        E  C8  C4   C8^3  C2  C2'  C2"  i   S8  S4  S8^3  SGH  SGV  SGD
        4, 2, 2, 2, 2, 2, 2, 0, 3, 3, 3, 1, 1, 1,
        0, 8, 4, 8, 2, 2, 2, 1, 8, 4, 8, 1, 2, 3,
        0, 1, 1, 3, 0, 2, 3, 0, 1, 1, 3, 0, 0, 0,
        1, 2, 2, 2, 1, 4, 4, 1, 2, 2, 2, 1, 4, 4,
        // c        E  S4  C2  C2'  SGD
        4, 3, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 4, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C3  C2'  i  S6  SGD
        4, 2, 2, 0, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 2, 1, 6, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 3, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  S8  C4  S8^3  C2  C2'  SGD
        4, 3, 2, 3, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 8, 4, 8, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 1, 4, 4, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C5  C5^2  C2'  i   S10  S10^3 SGD
        4, 2, 2, 2, 0, 3, 3, 1, 0, 0, 0, 0, 0, 0,
        0, 5, 5, 2, 1, 10, 10, 3, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 2, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 5, 1, 2, 2, 5, 0, 0, 0, 0, 0, 0,
        // c        E  S12  C6  S4  C3  S12^5  C2  C2'  SGD
        4, 3, 2, 3, 2, 3, 2, 2, 1, 0, 0, 0, 0, 0,
        0, 12, 6, 4, 3, 12, 2, 2, 3, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 5, 0, 2, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 2, 2, 1, 6, 6, 0, 0, 0, 0, 0,
        // c        E  C7 C7^2  C7^3  C2'  i   S14  S14^3 S14^5  SGD
        4, 2, 2, 2, 2, 0, 3, 3, 3, 1, 0, 0, 0, 0,
        0, 7, 7, 7, 2, 1, 14, 14, 14, 3, 0, 0, 0, 0,
        0, 1, 2, 3, 2, 0, 1, 3, 5, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 7, 1, 2, 2, 2, 7, 0, 0, 0, 0,
        // c        E  S16  C8  S16^3  C4  S16^5 C8^3 S16^7  C2  C2' SGD
        4, 3, 2, 3, 2, 3, 2, 3, 2, 2, 1, 0, 0, 0,
        0, 16, 8, 16, 4, 16, 8, 16, 2, 2, 3, 0, 0, 0,
        0, 1, 1, 3, 1, 5, 3, 7, 0, 2, 0, 0, 0, 0,
        1, 2, 2, 2, 1, 4, 4, 4, 1, 8, 8, 0, 0, 0,
        // c        E  S4  C2
        4, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E   C3  i  S6
        4, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 1, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  S8  C4  S8^3  C2
        4, 3, 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 8, 4, 8, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C3  C2
        4, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 8, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C3  C2  i  S6  SGH
        4, 2, 2, 0, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 2, 1, 6, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 8, 3, 1, 8, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C3  C2 S4  SGD
        4, 2, 2, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 2, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 8, 3, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C3  C2  C4  C2'
        4, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 2, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 8, 3, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E  C3  C2  C4  C2'  i   S6  SGH  S4  SGD
        4, 2, 2, 2, 2, 0, 3, 1, 3, 1, 0, 0, 0, 0,
        0, 3, 2, 4, 2, 1, 6, 1, 4, 3, 0, 0, 0, 0,
        0, 1, 0, 1, 2, 0, 1, 0, 1, 0, 0, 0, 0, 0,
        1, 8, 3, 6, 6, 1, 8, 3, 6, 6, 0, 0, 0, 0,
        // c        E   C5  C5^2 C3  C2
        4, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 5, 5, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 12, 12, 20, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // c        E   C5  C5^2  C3  C2   i   S10  S10^3  S6   SGH
        4, 2, 2, 2, 2, 0, 3, 3, 3, 1, 0, 0, 0, 0,
        0, 5, 5, 3, 2, 1, 10, 10, 6, 1, 0, 0, 0, 0,
        0, 1, 2, 1, 0, 0, 1, 3, 1, 0, 0, 0, 0, 0,
        1, 12, 12, 20, 15, 1, 12, 12, 20, 15, 0, 0, 0, 0
    };



    return convert_fortran_to_cpp_3d<14, 4, 55>(fortran_1d, 14, 4, 55);
}




// Assuming element_symbols_ is a member or global array of strings for elements 1-118 or so.
// For simplicity, we'll use a static array here.

// Additional helpers needed
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
static void center3(const std::array<double, 3>& p1, const std::array<double, 3>& p2, const std::array<double, 3>& p3, std::array<double, 3>& p0) {
    p0[0] = (p1[0] + p2[0] + p3[0]) / 3.0;
    p0[1] = (p1[1] + p2[1] + p3[1]) / 3.0;
    p0[2] = (p1[2] + p2[2] + p3[2]) / 3.0;
}

static void normal(const std::array<double, 3>& p1, const std::array<double, 3>& p2, std::array<double, 3>& v1, const std::array<double, 3>& p0) {
    std::array<double, 3> vec1 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    std::array<double, 3> vec2 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    symm_crossp(vec1, vec2, v1); // Use existing crossp
}
#pragma GCC diagnostic pop



void symclass(int natoms, int nprm, const std::vector<std::vector<int>>& nper,
              int& nseq, std::vector<int>& nccl, std::vector<std::vector<int>>& nscl,
              const std::vector<int>& nat, const std::vector<std::string>& symb, int nout) {
    // Suppress unused parameter warnings
    (void)nat;
    (void)symb;
    (void)nout;

    nseq = 0;
    
    // Outer loop - equivalent to Fortran's labeled loop
    for (int i = 0; i < natoms; i++) {  // i from 1 to natoms in Fortran (0-based in C++)
        bool continue_outer = false;
        
        // Check if atom i is already in an existing class
        for (int j = 0; j < nseq; j++) {  // j from 1 to nseq in Fortran
            for (int k = 0; k < nccl[j]; k++) {  // k from 1 to nccl(j) in Fortran
                if (nscl[k][j] == i + 1) {  // Fortran uses 1-based indexing
                    continue_outer = true;
                    break;
                }
            }
            if (continue_outer) break;
        }
        
        if (continue_outer) continue;
        
        // Start a new equivalence class
        nseq++;
        nccl[nseq - 1] = 0;  // Adjust for 0-based indexing
        
        // Inner loop - build the equivalence class
        for (int j = 0; j < nprm; j++) {  // j from 1 to nprm in Fortran
            bool continue_inner = false;
            
            // Check if nper(i,j) is already in the current class
            for (int k = 0; k < nccl[nseq - 1]; k++) {  // k from 1 to nccl(nseq) in Fortran
                if (nper[i][j] == nscl[k][nseq - 1]) {
                    continue_inner = true;
                    break;
                }
            }
            
            if (continue_inner) continue;
            
            // Add nper(i,j) to the current class
            nscl[nccl[nseq - 1]][nseq - 1] = nper[i][j];
            nccl[nseq - 1]++;
        }
    }
    
    // Sort each equivalence class
    for (int i = 0; i < nseq; i++) {  // i from 1 to nseq in Fortran
        for (int j = 0; j < nccl[i] - 1; j++) {  // j from 1 to nccl(i)-1 in Fortran
            int ii = j;
            
            // Find minimum element from j+1 to nccl(i)
            for (int k = j + 1; k < nccl[i]; k++) {  // k from j+1 to nccl(i) in Fortran
                if (nscl[k][i] < nscl[ii][i]) {
                    ii = k;
                }
            }
            
            // Swap if necessary
            if (ii != j) {
                int itemp = nscl[j][i];
                nscl[j][i] = nscl[ii][i];
                nscl[ii][i] = itemp;
            }
        }
    }
}


void sym_elements(int natoms, const std::vector<int>& nat,
                 const std::vector<std::vector<double>>& coord,
                 const std::vector<std::string>& symb, double delta,
                 int& norder, int& ni, int& nsg, int& ncr, int& nsr, int& np,
                 std::vector<std::vector<double>>& symn,
                 std::vector<std::vector<int>>& nsym, int nout, int& nprm,
                 std::vector<std::vector<int>>& nper, int& nseq,
                 std::vector<int>& nccl, std::vector<std::vector<int>>& nscl) {
    // Suppress unused parameter warnings
    (void)nseq;
    (void)nccl;
    (void)nscl;

    const int nmax = 150;
    
    // Variable declarations matching Fortran
    double delta2 = 0.0;
    double delta3 = 0.0;
    bool symcen = false;
    bool linear = false;
    bool planar = false;
    int neq;
    int nrot = 0;
    int icent = 0;
    
    // Arrays
    std::vector<int> ntrans(natoms);
    std::vector<int> meq(natoms);
    std::vector<std::vector<int>> ieq(10, std::vector<int>(natoms));
    std::vector<std::array<double, 3>> sigman(nmax);
    std::vector<std::array<double, 3>> rotn(nmax);
    std::vector<double> rota(nmax);
    
    // Points and vectors
    std::array<double, 3> p0, p1, p2, p3;
    std::array<double, 3> v1, v2, v3, v0;
    std::array<double, 3> a, b, c;
    
    // Generate PI
    const double pi = 4.0 * std::atan(1.0);
    
    // Initialize permutation arrays
    nprm = 1;
    for (int i = 0; i < natoms; ++i) {
        nper[i][0] = i + 1; // 1-based indexing
    }
    
    // Initialize symn and nsym
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < nmax; ++j) {
            symn[i][j] = 0.0;
        }
    }
    for (int i = 0; i < nmax; ++i) {
        for (int j = 0; j < 5; ++j) {
            nsym[i][j] = 0;
        }
    }
    
    // Partitioning atoms into equivalence classes
    neq = 1;
    meq[0] = 1;
    ieq[0][0] = 1; // 1-based
    
    for (int i = 1; i < natoms; ++i) {
        int nati = nat[i];
        bool found = false;
        
        for (int j = 0; j < neq; ++j) {
            int natj = nat[ieq[j][0] - 1]; // Convert to 0-based for nat array
            if (nati == natj) {
                meq[j]++;
                ieq[j][meq[j] - 1] = i + 1; // 1-based
                found = true;
                break;
            }
        }
        
        if (!found) {
            meq[neq] = 1;
            ieq[neq][0] = i + 1; // 1-based
            neq++;
        }
    }
    
    if (nout == 2) {
        std::cout << "\n-- Equivalence classes of atoms: " << neq << "\n";
        for (int i = 0; i < neq; ++i) {
            std::cout << "\n     #" << (i+1) << " (atom " 
                      << symb[nat[ieq[i][0] - 1] - 1] << ")" << "\n";
            std::cout << "     ";
            for (int j = 0; j < meq[i]; ++j) {
                std::cout << std::setw(4) << ieq[i][j];
            }
            std::cout << "\n";
        }
    }
    
    nsg = 0;
    
    // Centre of symmetry
    int nc;
    double del = 0.0;
    symm_inversion(natoms, nat, coord, delta, nc, ntrans, del);
    if (nc == natoms) {
        symcen = true;
        if (del > delta3) delta3 = del;
        nprm = 2;
        for (int i = 0; i < natoms; ++i) {
            nper[i][1] = ntrans[i] + 1; // Convert to 1-based
        }
    }
    
    // Initialize nsym[0]
    for (int j = 0; j < 5; ++j) {
        nsym[0][j] = 0;
    }
    
    if (symcen) {
        if (nout >= 1) std::cout << "\n-- CENTRE OF SYMMETRY: {i}" << "\n";
        nsym[0][1] = 1;
    }
    
    // Find atom at center of mass
    for (int i = 0; i < natoms; ++i) {
        p0[0] = coord[0][i];
        p0[1] = coord[1][i];
        p0[2] = coord[2][i];
        double sp = std::sqrt(symm_dot(p0.data(), p0.data(), 3));
        if (sp <= delta) {
            icent = i + 1; // 1-based
            if (sp > delta3) delta3 = sp;
        }
    }
    
    if (icent != 0 && nout == 2) {
        std::cout << "\n-- Atom " << symb[nat[icent - 1] - 1] 
                  << " (" << icent << ") in the COM" << "\n";
    }
    nsym[0][2] = icent;
    if (icent > 0) nsym[0][4] = 1;
    
    // Check for linear molecule
    linear = true;
    for (int i = 0; i < natoms - 1 && linear; ++i) {
        if (i + 1 == icent) continue; // Skip if this is the center atom

        p1[0] = coord[0][i];
        p1[1] = coord[1][i];
        p1[2] = coord[2][i];

        for (int j = i + 1; j < natoms && linear; ++j) {
            if (j + 1 == icent) continue; // Skip if this is the center atom

            p2[0] = coord[0][j];
            p2[1] = coord[1][j];
            p2[2] = coord[2][j];

            symm_crossp(p1, p2, p0);
            double vn = std::sqrt(symm_dot(p0.data(), p0.data(), 3));
            if (vn > delta) {
                linear = false;
                v1[0] = p0[0] / vn;
                v1[1] = p0[1] / vn;
                v1[2] = p0[2] / vn;
                break;
            }
            if (vn > delta3) delta3 = vn;
        }
    }

    bool is_linear = linear;
    if (is_linear) {
        if (nout >= 1) std::cout << "\n-- LINEAR MOLECULE" << "\n";

        if (symcen) {
            if (nout >= 1) {
                std::cout << "\n-- The structure should belong to the Dinf_h point group." << "\n";
                std::cout << "\n-- PLANES OF SYMMETRY --" << "\n";
            }

            nsg = 1;
            nsym[1][0] = 1; // SG
            nsym[1][1] = 0;
            nsym[1][2] = 0;
            nsym[1][3] = 2;
            nsym[1][4] = natoms;
            for (int k = 0; k < 3; ++k) {
                symn[k][1] = 0.0;
            }

            if (nout >= 1) std::cout << "\n-- Infinite planes" << "\n";
            if (nout == 2) std::cout << "     All atoms included." << "\n";

            if (nout >= 1) std::cout << "\n-- Distinct PROPER ROTATIONAL AXES --" << "\n";

            ncr = 2;
            nsym[2][0] = 2; // C2
            nsym[2][1] = 0;
            nsym[2][2] = 1;
            nsym[2][3] = 1;
            nsym[2][4] = natoms;

            // Set axis direction
            int ref_atom = (icent != 1) ? 0 : 1;
            double norm = std::sqrt(coord[0][ref_atom]*coord[0][ref_atom] +
                                   coord[1][ref_atom]*coord[1][ref_atom] +
                                   coord[2][ref_atom]*coord[2][ref_atom]);
            for (int k = 0; k < 3; ++k) {
                symn[k][2] = coord[k][ref_atom] / norm;
            }

            nsym[3][0] = 2;
            nsym[3][1] = 2;
            nsym[3][2] = 1;
            nsym[3][3] = 2;
            nsym[3][4] = 0;
            for (int k = 0; k < 3; ++k) {
                symn[k][3] = 0.0;
            }

            if (nout >= 1) {
                std::cout << "\n-- Axis #1: C(" << std::fixed << std::setprecision(2) << 0.0 << ")" << "\n";
            }
            if (nout == 2) {
                std::cout << "  d: " << std::fixed << std::setprecision(5);
                for (int k = 0; k < 3; ++k) {
                    std::cout << std::setw(12) << symn[k][2];
                }
                std::cout << "\n     All atoms included." << "\n";
            }

            if (nout >= 1) {
                std::cout << "\n-- Axis #2: C(" << std::fixed << std::setprecision(2) << 180.0 << ")" << "\n";
            }
            if (nout == 2) {
                std::cout << "  d: " << std::fixed << std::setprecision(5);
                for (int k = 0; k < 3; ++k) {
                    std::cout << std::setw(12) << symn[k][3];
                }
                std::cout << "\n     Atoms included:" << "\n";
                if (icent != 0) {
                    std::cout << "          " << symb[nat[nsym[0][2] - 1] - 1]
                              << " (" << nsym[0][2] << ")" << "\n";
                }
            }

            nsr = 1;
            nsym[4][0] = 3; // S
            nsym[4][1] = 0;
            nsym[4][2] = 1;
            nsym[4][3] = 0;
            nsym[4][4] = 1;

            for (int k = 0; k < 3; ++k) {
                symn[k][4] = coord[k][ref_atom] / norm;
            }

            ni = 1;
            norder = -1;
            np = -1;

            if (nout >= 1) std::cout << "\n-- Number of symmetry operations = infinite" << "\n";

        } else {
            if (nout >= 1) {
                std::cout << "\n-- The structure should belong to the Cinf_v point group." << "\n";
                std::cout << "\n-- PLANES OF SYMMETRY --" << "\n";
            }

            nsg = 1;
            nsym[1][0] = 1;
            nsym[1][1] = 0;
            nsym[1][2] = 0;
            nsym[1][3] = 2;
            nsym[1][4] = natoms;
            for (int k = 0; k < 3; ++k) {
                symn[k][1] = 0.0;
            }

            if (nout >= 1) std::cout << "\n-- Infinite planes" << "\n";
            if (nout == 2) std::cout << "     All atoms included." << "\n";

            if (nout >= 1) std::cout << "\n-- Distinct PROPER ROTATIONAL AXES --" << "\n";

            ncr = 1;
            nsym[2][0] = 2;
            nsym[2][1] = 0;
            nsym[2][2] = 1;
            nsym[2][3] = 1;
            nsym[2][4] = natoms;

            int ref_atom = (icent != 1) ? 0 : 1;
            double norm = std::sqrt(coord[0][ref_atom]*coord[0][ref_atom] +
                                   coord[1][ref_atom]*coord[1][ref_atom] +
                                   coord[2][ref_atom]*coord[2][ref_atom]);
            for (int k = 0; k < 3; ++k) {
                symn[k][2] = coord[k][ref_atom] / norm;
            }

            if (nout >= 1) {
                std::cout << "\n-- Axis #1: C(" << std::fixed << std::setprecision(2) << 0.0 << ")" << "\n";
            }
            if (nout == 2) {
                std::cout << "  d: " << std::fixed << std::setprecision(5);
                for (int k = 0; k < 3; ++k) {
                    std::cout << std::setw(12) << symn[k][2];
                }
                std::cout << "\n     All atoms included." << "\n";
            }

            nsr = 0;
            norder = -1;
            ni = 0;
            np = -1;

            if (nout >= 1) std::cout << "\n-- Number of symmetry operations = infinite" << "\n";
        }

        // Skip the rest of the function for linear molecules
    }

    if (!is_linear) {
        // Check for planar molecule
    planar = true;
    for (int i = 0; i < natoms; ++i) {
        p3[0] = coord[0][i];
        p3[1] = coord[1][i];
        p3[2] = coord[2][i];
        double sp = symm_dot(v1.data(), p3.data(), 3);
        if (std::abs(sp) > delta) {
            planar = false;
            break;
        }
        if (std::abs(sp) > delta3) delta3 = std::abs(sp);
    }
    
    if (planar) {
        nsg = nsg + 1;
        sigman[nsg - 1][0] = v1[0];
        sigman[nsg - 1][1] = v1[1];
        sigman[nsg - 1][2] = v1[2];
        
        if (nout >= 1) std::cout << "\n-- PLANAR MOLECULE" << "\n";
        if (nout == 2) {
            std::cout << "  n: " << std::fixed << std::setprecision(5);
            for (int k = 0; k < 3; ++k) {
                std::cout << std::setw(12) << v1[k];
            }
            std::cout << "\n";
        }
        
        if (symcen && planar) {
            for (int i = 12; i >= 2; --i) {
                double alpha = 2.0 * pi / static_cast<double>(i);
                double sp = alpha * 180.0 / pi;
                double sina = std::sin(alpha);
                double cosa = std::cos(alpha);
                
                symm_rotate(natoms, nat, coord, v1, sina, cosa, delta, nc, ntrans, del);
                if (nc == natoms) {
                    add_Cn(nrot, rotn, rota, v1, p3, sp, delta);
                    add_perm(natoms, ntrans, nprm, nper);
                    if (del > delta3) delta3 = del;
                }
            }
        }
    }
    
    // Planes of symmetry from equivalent atom pairs
    for (int i = 0; i < neq; ++i) {
        int meqi = meq[i];
        if (meqi == 1) continue;
        
        for (int j1 = 0; j1 < meqi - 1; ++j1) {
            int i1 = ieq[i][j1] - 1; // Convert to 0-based
            p1[0] = coord[0][i1];
            p1[1] = coord[1][i1];
            p1[2] = coord[2][i1];
            
            for (int j2 = j1 + 1; j2 < meqi; ++j2) {
                int i2 = ieq[i][j2] - 1; // Convert to 0-based
                p2[0] = coord[0][i2];
                p2[1] = coord[1][i2];
                p2[2] = coord[2][i2];
                
                // Midpoint
                p0[0] = (p1[0] + p2[0]) / 2.0;
                p0[1] = (p1[1] + p2[1]) / 2.0;
                p0[2] = (p1[2] + p2[2]) / 2.0;
                
                // Direction vector
                v1[0] = p2[0] - p0[0];
                v1[1] = p2[1] - p0[1];
                v1[2] = p2[2] - p0[2];
                
                double vn = std::sqrt(symm_dot(v1.data(), v1.data(), 3));
                if (vn > delta) {
                    v1[0] = v1[0] / vn;
                    v1[1] = v1[1] / vn;
                    v1[2] = v1[2] / vn;
                    
                    std::array<double, 3> neg_p0 = {-p0[0], -p0[1], -p0[2]};
                    double sp = symm_dot(v1.data(), neg_p0.data(), 3);
                    if (std::abs(sp) < delta) {
                        symm_reflect(natoms, nat, coord, v1, p0, delta, nc, ntrans, del);
                        if (nc == natoms) {
                            if (del > delta3) delta3 = del;
                            add_SG(nsg, sigman, v1, delta);
                            add_perm(natoms, ntrans, nprm, nper);
                        }
                    }
                }
                
                // Cross product plane
                symm_crossp(p1, p2, v2);
                vn = std::sqrt(symm_dot(v2.data(), v2.data(), 3));
                if (vn > delta) {
                    v2[0] = v2[0] / vn;
                    v2[1] = v2[1] / vn;
                    v2[2] = v2[2] / vn;
                    
                    std::array<double, 3> neg_p0 = {-p0[0], -p0[1], -p0[2]};
                    double sp = symm_dot(v2.data(), neg_p0.data(), 3);
                    if (std::abs(sp) < delta) {
                        symm_reflect(natoms, nat, coord, v2, p0, delta, nc, ntrans, del);
                        if (nc == natoms) {
                            if (del > delta3) delta3 = del;
                            add_SG(nsg, sigman, v2, delta);
                            add_perm(natoms, ntrans, nprm, nper);
                        }
                    }
                }
                
                // Third plane
                symm_crossp(v1, v2, v3);
                vn = std::sqrt(symm_dot(v3.data(), v3.data(), 3));
                if (vn > delta) {
                    v3[0] = v3[0] / vn;
                    v3[1] = v3[1] / vn;
                    v3[2] = v3[2] / vn;
                    
                    std::array<double, 3> neg_p0 = {-p0[0], -p0[1], -p0[2]};
                    double sp = symm_dot(v3.data(), neg_p0.data(), 3);
                    if (std::abs(sp) < delta) {
                        symm_reflect(natoms, nat, coord, v3, p0, delta, nc, ntrans, del);
                        if (nc == natoms) {
                            if (del > delta3) delta3 = del;
                            add_SG(nsg, sigman, v3, delta);
                            add_perm(natoms, ntrans, nprm, nper);
                        }
                    }
                }
            }
        }
    }
    
    // Output planes of symmetry
    if (nout >= 1) std::cout << "\n-- PLANES OF SYMMETRY --" << "\n";
    
    for (int i = 0; i < nsg; ++i) {
        if (nout >= 1) std::cout << "\n-- Plane #" << (i + 1) << "\n";
        
        v1[0] = sigman[i][0];
        v1[1] = sigman[i][1];
        v1[2] = sigman[i][2];
        
        if (nout == 2) {
            std::cout << "  n: " << std::fixed << std::setprecision(5);
            for (int k = 0; k < 3; ++k) {
                std::cout << std::setw(12) << v1[k];
            }
            std::cout << "\n     Atoms included:" << "\n";
        }
        
        int m = 0;
        for (int j = 0; j < natoms; ++j) {
            p3[0] = coord[0][j];
            p3[1] = coord[1][j];
            p3[2] = coord[2][j];
            double sp = symm_dot(v1.data(), p3.data(), 3);
            if (std::abs(sp) <= delta) {
                if (nout == 2) {
                    std::cout << "          " << symb[nat[j] - 1] 
                              << " (" << (j + 1) << ")" << "\n";
                }
                m++;
                if (std::abs(sp) > delta3) delta3 = std::abs(sp);
            }
        }
        
        // Store in symn and nsym arrays
        for (int k = 0; k < 3; ++k) {
            symn[k][i + 1] = v1[k];
        }
        nsym[i + 1][0] = 1; // SG
        nsym[i + 1][1] = 0;
        nsym[i + 1][2] = 0;
        nsym[i + 1][3] = 0; // Will be classified later
        nsym[i + 1][4] = m;
    }
    
    // Proper rotations from plane intersections
    if (nout >= 1) std::cout << "\n-- Proper rotations due to the centre and planes of symmetry --" << "\n";
    
    for (int i = 0; i < nsg - 1; ++i) {
        v1[0] = sigman[i][0];
        v1[1] = sigman[i][1];
        v1[2] = sigman[i][2];
        
        for (int j = i + 1; j < nsg; ++j) {
            v2[0] = sigman[j][0];
            v2[1] = sigman[j][1];
            v2[2] = sigman[j][2];
            
            symm_crossp(v1, v2, v3);
            double vn = std::sqrt(symm_dot(v3.data(), v3.data(), 3));
            if (vn > delta) {
                v3[0] = v3[0] / vn;
                v3[1] = v3[1] / vn;
                v3[2] = v3[2] / vn;
                
                double sp = symm_dot(v1.data(), v2.data(), 3);
                sp = std::acos(sp) * 180.0 / pi;
                if (sp > 90.0) sp = 180.0 - sp;
                sp = 2.0 * sp;
                
                add_Cn(nrot, rotn, rota, v3, p3, sp, delta);
            }
        }
    }
    
    // Output rotations from plane intersections
    for (int i = 0; i < nrot; ++i) {
        int m = 0; // Suppress unused variable warning
        (void)m;
        double sp = rota[i];
        if (nout >= 1) {
            std::cout << "\n-- Rotation #" << (i + 1) << ": C(" 
                      << std::fixed << std::setprecision(2) << sp << ")" << "\n";
        }
        
        v1[0] = rotn[i][0];
        v1[1] = rotn[i][1];
        v1[2] = rotn[i][2];
        
        if (nout == 2) {
            std::cout << "  d: " << std::fixed << std::setprecision(5);
            for (int k = 0; k < 3; ++k) {
                std::cout << std::setw(12) << v1[k];
            }
            std::cout << "\n     Atoms included:" << "\n";
        }
        
        for (int j = 0; j < natoms; ++j) {
            p3[0] = coord[0][j];
            p3[1] = coord[1][j];
            p3[2] = coord[2][j];
            
            v2[0] = p0[0] + v1[0];
            v2[1] = p0[1] + v1[1];
            v2[2] = p0[2] + v1[2];
            
            symm_crossp(v2, p3, v0);
            double vn = std::sqrt(symm_dot(v0.data(), v0.data(), 3));
            if (std::abs(vn) <= delta) {
                if (nout == 2) {
                    std::cout << "          " << symb[nat[j] - 1] 
                              << " (" << (j + 1) << ")" << "\n";
                }
                m++;
                if (std::abs(vn) > delta3) delta3 = std::abs(vn);
            }
        }
    }
    
    // Proper rotational axes - single atoms
    if (nout >= 1) std::cout << "\n-- Distinct PROPER ROTATIONAL AXES --" << "\n";
    
    // Cn (for each atom)
    for (int i = 0; i < neq; ++i) {
        int meqi = meq[i];
        for (int j = 0; j < meqi; ++j) {
            int i1 = ieq[i][j] - 1; // Convert to 0-based
            p0[0] = coord[0][i1];
            p0[1] = coord[1][i1];
            p0[2] = coord[2][i1];
            
            double vn = std::sqrt(symm_dot(p0.data(), p0.data(), 3));
            if (vn < delta) continue;
            
            v1[0] = p0[0] / vn;
            v1[1] = p0[1] / vn;
            v1[2] = p0[2] / vn;
            
            for (int k = 12; k >= 2; --k) {
                double alpha = 2.0 * pi / static_cast<double>(k);
                double sp = alpha * 180.0 / pi;
                double sina = std::sin(alpha);
                double cosa = std::cos(alpha);
                
                symm_rotate(natoms, nat, coord, v1, sina, cosa, delta, nc, ntrans, del);
                if (nc == natoms) {
                    add_Cn(nrot, rotn, rota, v1, p3, sp, delta);
                    add_perm(natoms, ntrans, nprm, nper);
                    if (del > delta3) delta3 = del;
                }
            }
        }
    }
    
    // Cn (for each pair of atoms)
    for (int i = 0; i < neq; ++i) {
        int meqi = meq[i];
        if (meqi < 2) continue;
        
        for (int j1 = 0; j1 < meqi - 1; ++j1) {
            int i1 = ieq[i][j1] - 1; // Convert to 0-based
            p1[0] = coord[0][i1];
            p1[1] = coord[1][i1];
            p1[2] = coord[2][i1];
            
            for (int j2 = j1 + 1; j2 < meqi; ++j2) {
                int i2 = ieq[i][j2] - 1; // Convert to 0-based
                p2[0] = coord[0][i2];
                p2[1] = coord[1][i2];
                p2[2] = coord[2][i2];
                
                p0[0] = (p1[0] + p2[0]) / 2.0;
                p0[1] = (p1[1] + p2[1]) / 2.0;
                p0[2] = (p1[2] + p2[2]) / 2.0;
                
                double vn = std::sqrt(symm_dot(p0.data(), p0.data(), 3));
                if (vn < delta) continue;
                
                v1[0] = p0[0] / vn;
                v1[1] = p0[1] / vn;
                v1[2] = p0[2] / vn;
                
                double alpha = pi;
                double sp = 180.0;
                double sina = std::sin(alpha);
                double cosa = std::cos(alpha);
                
                symm_rotate(natoms, nat, coord, v1, sina, cosa, delta, nc, ntrans, del);
                if (nc == natoms) {
                    add_Cn(nrot, rotn, rota, v1, p3, sp, delta);
                    add_perm(natoms, ntrans, nprm, nper);
                    if (del > delta3) delta3 = del;
                }
            }
        }
    }
    
    // Cn (n > 2) from triplets
    for (int i = 0; i < neq; ++i) {
        int meqi = meq[i];
        if (meqi < 3) continue;
        
        for (int j1 = 0; j1 < meqi - 2; ++j1) {
            int i1 = ieq[i][j1] - 1; // Convert to 0-based
            a[0] = coord[0][i1];
            a[1] = coord[1][i1];
            a[2] = coord[2][i1];
            
            for (int j2 = j1 + 1; j2 < meqi - 1; ++j2) {
                int i2 = ieq[i][j2] - 1; // Convert to 0-based
                b[0] = coord[0][i2];
                b[1] = coord[1][i2];
                b[2] = coord[2][i2];
                
                p1[0] = (a[0] + b[0]) / 2.0;
                p1[1] = (a[1] + b[1]) / 2.0;
                p1[2] = (a[2] + b[2]) / 2.0;
                
                p3[0] = b[0] - p1[0];
                p3[1] = b[1] - p1[1];
                p3[2] = b[2] - p1[2];
                
                double vn = std::sqrt(symm_dot(p3.data(), p3.data(), 3));
                v1[0] = p3[0] / vn;
                v1[1] = p3[1] / vn;
                v1[2] = p3[2] / vn;
                
                for (int j3 = j2 + 1; j3 < meqi; ++j3) {
                    int i3 = ieq[i][j3] - 1; // Convert to 0-based
                    c[0] = coord[0][i3];
                    c[1] = coord[1][i3];
                    c[2] = coord[2][i3];
                    
                    p2[0] = (b[0] + c[0]) / 2.0;
                    p2[1] = (b[1] + c[1]) / 2.0;
                    p2[2] = (b[2] + c[2]) / 2.0;
                    
                    p3[0] = c[0] - p2[0];
                    p3[1] = c[1] - p2[1];
                    p3[2] = c[2] - p2[2];
                    
                    vn = std::sqrt(symm_dot(p3.data(), p3.data(), 3));
                    v2[0] = p3[0] / vn;
                    v2[1] = p3[1] / vn;
                    v2[2] = p3[2] / vn;
                    
                    symm_crossp(v1, v2, v3);
                    vn = std::sqrt(symm_dot(v3.data(), v3.data(), 3));
                    if (vn < delta) continue;
                    
                    v3[0] = v3[0] / vn;
                    v3[1] = v3[1] / vn;
                    v3[2] = v3[2] / vn;
                    
                    double sp = std::acos(symm_dot(v1.data(), v2.data(), 3));
                    if (std::abs(sp) < delta) continue;
                    
                    int m = static_cast<int>(2.0 * pi / sp + delta);
                    if ((m * sp) < (2.0 * pi - delta)) continue;
                    if ((m < 3) || (m > meqi)) continue;
                    
                    double alpha = sp;
                    sp = alpha * 180.0 / pi;
                    double sina = std::sin(alpha);
                    double cosa = std::cos(alpha);
                    
                    symm_rotate(natoms, nat, coord, v3, sina, cosa, delta, nc, ntrans, del);
                    if (nc == natoms) {
                        add_Cn(nrot, rotn, rota, v3, p3, sp, delta);
                        add_perm(natoms, ntrans, nprm, nper);
                        if (del > delta3) delta3 = del;
                    }
                }
            }
        }
    }
    
    // Output rotation axes found so far
    for (int i = 0; i < nrot; ++i) {
        int m = 0; // Suppress unused variable warning
        (void)m;
        double sp = rota[i];
        if (nout >= 1) {
            std::cout << "\n-- Axis #" << (i + 1) << ": C(" 
                      << std::fixed << std::setprecision(2) << sp << ")" << "\n";
        }
        
        v1[0] = rotn[i][0];
        v1[1] = rotn[i][1];
        v1[2] = rotn[i][2];
        
        if (nout == 2) {
            std::cout << "  d: " << std::fixed << std::setprecision(5);
            for (int k = 0; k < 3; ++k) {
                std::cout << std::setw(12) << v1[k];
            }
            std::cout << "\n     Atoms included:" << "\n";
        }
        
        for (int j = 0; j < natoms; ++j) {
            p3[0] = coord[0][j];
            p3[1] = coord[1][j];
            p3[2] = coord[2][j];
            
            symm_crossp(v1, p3, v0);
            double vn = std::sqrt(symm_dot(v0.data(), v0.data(), 3));
            if (std::abs(vn) <= delta) {
                if (nout == 2) {
                    std::cout << "          " << symb[nat[j] - 1] 
                              << " (" << (j + 1) << ")" << "\n";
                }
                m++;
                if (std::abs(vn) > delta3) delta3 = std::abs(vn);
            }
        }
    }
    
    // Generate all proper rotations
    if (nout >= 1) std::cout << "\n-- PROPER ROTATIONAL AXES & ROTATIONS --" << "\n";
    
    int nsgi = nsg + 1;
    int ii = 0;
    
    for (int i = 0; i < nrot; ++i) {
        v1[0] = rotn[i][0];
        v1[1] = rotn[i][1];
        v1[2] = rotn[i][2];
        
        for (int k = 12; k >= 2; --k) {
            double alpha = 2.0 * pi / static_cast<double>(k);
            double sina = std::sin(alpha);
            double cosa = std::cos(alpha);
            
            symm_rotate(natoms, nat, coord, v1, sina, cosa, delta, nc, ntrans, del);
            if (nc == natoms) {
                if (del > delta3) delta3 = del;
                add_perm(natoms, ntrans, nprm, nper);
                ii++;
                
                int m = 0;
                for (int j = 0; j < natoms; ++j) {
                    p3[0] = coord[0][j];
                    p3[1] = coord[1][j];
                    p3[2] = coord[2][j];
                    
                    symm_crossp(v1, p3, v0);
                    double vn = std::sqrt(symm_dot(v0.data(), v0.data(), 3));
                    if (std::abs(vn) <= delta) {
                        m++;
                        if (std::abs(vn) > delta3) delta3 = std::abs(vn);
                    }
                }
                
                // Store vector n
                for (int idx = 0; idx < 3; ++idx) {
                    symn[idx][nsgi + ii - 1] = v1[idx];
                }
                
                // Cn
                nsym[nsgi + ii - 1][0] = 2;
                // n-fold (n ^ 1) 
                nsym[nsgi + ii - 1][1] = k;
                nsym[nsgi + ii - 1][2] = 1;
                // Principal axis
                nsym[nsgi + ii - 1][3] = 0;
                // Number of atoms included (unmoved atoms) 
                nsym[nsgi + ii - 1][4] = m;
                
                if (nout >= 1) {
                    std::cout << "-- #" << (i + 1) << "-" << ii << ": C(" << k << ")" << "\n";
                }
                
                if (k > 2) {
                    for (int kk = 2; kk < k; ++kk) {
                        double alpha2 = static_cast<double>(kk) * alpha;
                        double sina2 = std::sin(alpha2);
                        double cosa2 = std::cos(alpha2);
                        
                        int nc2;
                        symm_rotate(natoms, nat, coord, v1, sina2, cosa2, delta, nc2, ntrans, del);
                        if (nc2 == natoms) {
                            if (del > delta3) delta3 = del;
                            add_perm(natoms, ntrans, nprm, nper);
                            ii++;
                            
                            // Store vector n
                            for (int idx = 0; idx < 3; ++idx) {
                                symn[idx][nsgi + ii - 1] = v1[idx];
                            }
                            
                            // Cn
                            nsym[nsgi + ii - 1][0] = 2;
                            // Principal axis
                            nsym[nsgi + ii - 1][3] = 0;
                            // Number of atoms included 
                            nsym[nsgi + ii - 1][4] = m;
                            
                            int ngcd = symm_igcd(k, kk);
                            nsym[nsgi + ii - 1][1] = k / ngcd;
                            nsym[nsgi + ii - 1][2] = kk / ngcd;
                            
                            if (nout >= 1) {
                                std::cout << "-- #" << (i + 1) << "-" << ii << ": C(" << k 
                                          << " ^" << kk << ")" << "\n";
                            }
                        }
                    }
                    break;
                }
            }
        }
    }
    ncr = ii;
    
    // Improper rotational axes
    if (nout >= 1) std::cout << "\n-- IMPROPER ROTATIONAL AXES & ROTATIONS --" << "\n";
    
    int nsgicn = nsgi + ncr;
    ii = 0;
    
    for (int i = 0; i < nrot; ++i) {
        v1[0] = rotn[i][0];
        v1[1] = rotn[i][1];
        v1[2] = rotn[i][2];
        
        for (int k = 24; k >= 2; --k) {
            double alpha = 2.0 * pi / static_cast<double>(k);
            double sina = std::sin(alpha);
            double cosa = std::cos(alpha);
            
            symm_srotate(natoms, nat, coord, v1, sina, cosa, delta, nc, ntrans, del);
            if (nc == natoms && k > 2) {
                if (del > delta3) delta3 = del;
                add_perm(natoms, ntrans, nprm, nper);
                ii++;
                
                int m = 0;
                if (icent > 0) m = 1;
                
                // Store vector n
                for (int idx = 0; idx < 3; ++idx) {
                    symn[idx][nsgicn + ii - 1] = v1[idx];
                }
                
                // Sn
                nsym[nsgicn + ii - 1][0] = 3;
                // n-fold (n ^ 1) 
                nsym[nsgicn + ii - 1][1] = k;
                nsym[nsgicn + ii - 1][2] = 1;
                nsym[nsgicn + ii - 1][3] = 0;
                // Number of atoms included 
                nsym[nsgicn + ii - 1][4] = m;
                
                if (nout >= 1) {
                    std::cout << "-- #" << (i + 1) << "-" << ii << ": S(" << k << ")" << "\n";
                }
                
                int kv = k - 1;
                if (k % 2 != 0) kv = 2 * k - 1;
                
                for (int kk = 2; kk <= kv; ++kk) {
                    double alpha2 = static_cast<double>(kk) * alpha;
                    double sina2 = std::sin(alpha2);
                    double cosa2 = std::cos(alpha2);
                    
                    int nc2;
                    symm_srotate(natoms, nat, coord, v1, sina2, cosa2, delta, nc2, ntrans, del);
                    if ((nc2 == natoms) && (kk % 2 != 0) && (kk != k) && (symm_igcd(k, kk) == 1)) {
                        if (del > delta3) delta3 = del;
                        add_perm(natoms, ntrans, nprm, nper);
                        ii++;
                        
                        for (int idx = 0; idx < 3; ++idx) {
                            symn[idx][nsgicn + ii - 1] = v1[idx];
                        }
                        
                        nsym[nsgicn + ii - 1][0] = 3;
                        nsym[nsgicn + ii - 1][4] = m;
                        nsym[nsgicn + ii - 1][1] = k;
                        nsym[nsgicn + ii - 1][2] = kk;
                        
                        if (nout >= 1) {
                            std::cout << "-- #" << (i + 1) << "-" << ii << ": S(" << k 
                                      << "^" << kk << ")" << "\n";
                        }
                    }
                }
            }
        }
    }
    nsr = ii;
    ni = 0;
    if (symcen) ni = 1;
    norder = 1 + ni + nsg + ncr + nsr;
    
    if (nout >= 1) {
        std::cout << "\n-- Number of symmetry operations (including E) = " << norder << "\n";
    }
    
    // Determination of the principal axis
    np = 0;
    for (int i = nsg + 1; i < nsg + ncr + 1; ++i) {
        if ((nsym[i][1] > np) && (nsym[i][2] == 1)) {
            np = nsym[i][1];
        }
    }
    for (int i = nsg + 1; i < nsg + ncr + 1; ++i) {
        if ((nsym[i][1] == np) && (nsym[i][2] == 1)) {
            nsym[i][3] = 1;
        }
    }
    
    // Rotation axes: principal axes & orthogonal C2 axes
    int nnp = 0;
    for (int i = nsg + 1; i < nsg + ncr + 1; ++i) {
        if (nsym[i][3] == 1) nnp++;
    }
    
    if (nnp == 3 && np == 2) {
        for (int i = nsg + 1; i < nsg + ncr + 1; ++i) {
            nsym[i][3] = 2;
        }
        
        // D2, D2d, D2h
        if (nsg == 3) {
            // D2h - Fixed bug from original Fortran (nm was 0, changed to 2)
            int nm = 2;
            for (int i = 1; i < nsg + 1; ++i) {
                if (nsym[i][4] > nsym[nm][4]) nm = i;
            }
            nsym[nm][3] = 1;
            
            std::array<double, 3> v2;
            for (int k = 0; k < 3; ++k) {
                v2[k] = symn[k][nm];
            }
            
            for (int i = nsg + 1; i < nsg + ncr + 1; ++i) {
                std::array<double, 3> v1_temp;
                for (int k = 0; k < 3; ++k) {
                    v1_temp[k] = symn[k][i];
                }
                
                double vk = symm_dot(v1_temp.data(), v2.data(), 3);
                if (std::abs(std::abs(vk) - 1.0) <= delta) {
                    nsym[i][3] = 1;
                    if (std::abs(std::abs(vk) - 1.0) > delta2) {
                        delta2 = std::abs(std::abs(vk) - 1.0);
                    }
                }
            }
        } else if (nsg == 2) {
            // D2d
            for (int i = nsg + 1; i < nsg + ncr + 1; ++i) {
                std::array<double, 3> v1_temp;
                for (int k = 0; k < 3; ++k) {
                    v1_temp[k] = symn[k][i];
                }
                
                for (int j = nsg + ncr + 1; j < nsg + ncr + nsr + 1; ++j) {
                    std::array<double, 3> v2_temp;
                    for (int k = 0; k < 3; ++k) {
                        v2_temp[k] = symn[k][j];
                    }
                    
                    double vk = symm_dot(v1_temp.data(), v2_temp.data(), 3);
                    if (std::abs(std::abs(vk) - 1.0) <= delta) {
                        nsym[i][3] = 1;
                        if (std::abs(std::abs(vk) - 1.0) > delta2) {
                            delta2 = std::abs(std::abs(vk) - 1.0);
                        }
                    }
                }
            }
        }
    }
    
    // Orthogonal C2 axes
    for (int i = nsg + 1; i < nsg + ncr + 1; ++i) {
        if (nsym[i][3] != 1) continue;
        
        std::array<double, 3> v1_temp;
        for (int k = 0; k < 3; ++k) {
            v1_temp[k] = symn[k][i];
        }
        
        for (int j = nsg + 1; j < nsg + ncr + 1; ++j) {
            if ((nsym[j][1] == 2) && (nsym[j][2] == 1)) {
                std::array<double, 3> v2_temp;
                for (int k = 0; k < 3; ++k) {
                    v2_temp[k] = symn[k][j];
                }
                
                double vk = symm_dot(v1_temp.data(), v2_temp.data(), 3);
                if (std::abs(vk) < delta) {
                    nsym[j][3] = 2;
                    if (std::abs(vk) > delta2) delta2 = std::abs(vk);
                }
            }
        }
    }
    
    // Perpendicular C2 axes
    int maxcn = 0;
    if (ncr > 0) {
        for (int i = nsg + 1; i < nsg + ncr + 1; ++i) {
            if (nsym[i][3] != 2) continue;
            if (nsym[i][4] > maxcn) maxcn = nsym[i][4];
        }
        for (int i = nsg + 1; i < nsg + ncr + 1; ++i) {
            if (nsym[i][3] != 2) continue;
            if (nsym[i][4] < maxcn) nsym[i][3] = 3;
        }
    }
    
    // Planes of symmetry: SGH, SGV, SGD
    if (nsg > 0) {
        // SGH
        for (int i = nsg + 1; i < nsg + ncr + 1; ++i) {
            if (nsym[i][3] != 1) continue;
            
            std::array<double, 3> v1_temp;
            for (int k = 0; k < 3; ++k) {
                v1_temp[k] = symn[k][i];
            }
            
            for (int j = 1; j < nsg + 1; ++j) {
                std::array<double, 3> v2_temp;
                for (int k = 0; k < 3; ++k) {
                    v2_temp[k] = symn[k][j];
                }
                
                double vk = symm_dot(v1_temp.data(), v2_temp.data(), 3);
                if (std::abs(std::abs(vk) - 1.0) <= delta) {
                    nsym[j][3] = 1;
                    if (std::abs(std::abs(vk) - 1.0) > delta2) {
                        delta2 = std::abs(std::abs(vk) - 1.0);
                    }
                }
            }
        }
        
        // SGV
        for (int i = 1; i < nsg + 1; ++i) {
            if (nsym[i][3] == 1) continue;
            
            std::array<double, 3> v1_temp;
            for (int k = 0; k < 3; ++k) {
                v1_temp[k] = symn[k][i];
            }
            
            for (int j = nsg + 1; j < nsg + ncr + 1; ++j) {
                if (nsym[j][3] != 1) continue;
                
                std::array<double, 3> v2_temp;
                for (int k = 0; k < 3; ++k) {
                    v2_temp[k] = symn[k][j];
                }
                
                double vk = symm_dot(v1_temp.data(), v2_temp.data(), 3);
                if (std::abs(vk) < delta) {
                    nsym[i][3] = 2;
                    if (std::abs(vk) > delta2) delta2 = std::abs(vk);
                }
            }
        }
        
        int maxsg = 0;
        for (int i = 1; i < nsg + 1; ++i) {
            if (nsym[i][3] != 2) continue;
            if (nsym[i][4] > maxsg) maxsg = nsym[i][4];
        }
        for (int i = 1; i < nsg + 1; ++i) {
            if (nsym[i][3] != 2) continue;
            if (nsym[i][4] < maxsg) nsym[i][3] = 3;
        }
    } // end if (!is_linear)
    
    if (nout >= 1) std::cout << "\n-- SYMMETRY OPERATIONS --" << "\n";
    
    // COM and Inversion Center
    if (nout == 2) {
        if (nsym[0][1] == 0) {
            if (nsym[0][2] > 0) {
                std::cout << "               #1: COM    -- with atom " 
                          << symb[nat[nsym[0][2] - 1] - 1] << " (#" << nsym[0][2] << ")" << "\n";
            } else {
                std::cout << "               #1: COM" << "\n";
            }
        } else if (nsym[0][1] == 1) {
            if (nsym[0][2] > 0) {
                std::cout << "               #1: INVERSION CENTER  -- with atom " 
                          << symb[nat[nsym[0][2] - 1] - 1] << " (#" << nsym[0][2] << ")" << "\n";
            } else {
                std::cout << "               #1: INVERSION CENTER " << "\n";
            }
        }
        
        // SG
        for (int k = 1; k < nsg + 1; ++k) {
            std::string symel;
            if (nsym[k][0] == 1 && nsym[k][3] == 0) {
                symel = "SG ";
            } else if (nsym[k][0] == 1 && nsym[k][3] == 1) {
                symel = "SGH";
            } else if (nsym[k][0] == 1 && nsym[k][3] == 2) {
                symel = "SGV";
            } else if (nsym[k][0] == 1 && nsym[k][3] == 3) {
                symel = "SGD";
            }
            
            if (nsym[k][4] == 0) {
                std::cout << "               #" << (k + 1) << ": " << symel << "\n";
            } else {
                std::cout << "               #" << (k + 1) << ": " << symel 
                          << "     -- with " << nsym[k][4] << " unmoved atoms" << "\n";
            }
        }
        
        // C(n^k)
        for (int k = nsg + 1; k < nsg + ncr + 1; ++k) {
            if (nsym[k][2] > 1) {
                std::cout << "               #" << (k + 1) << ": C(" << nsym[k][1] 
                          << "^" << nsym[k][2] << ")" << "\n";
            } else {
                if (nsym[k][4] > 0) {
                    if (nsym[k][3] == 2) {
                        std::cout << "               #" << (k + 1) << ": C'(" << nsym[k][1] 
                                  << ")  -- with " << nsym[k][4] << " unmoved atoms" << "\n";
                    } else if (nsym[k][3] == 3) {
                        std::cout << "               #" << (k + 1) << ": C\"(" << nsym[k][1] 
                                  << ")  -- with " << nsym[k][4] << " unmoved atoms" << "\n";
                    } else {
                        std::cout << "               #" << (k + 1) << ": C(" << nsym[k][1] 
                                  << ")   -- with " << nsym[k][4] << " unmoved atoms" << "\n";
                    }
                } else {
                    if (nsym[k][3] == 2) {
                        std::cout << "               #" << (k + 1) << ": C'(" << nsym[k][1] << ")" << "\n";
                    } else if (nsym[k][3] == 3) {
                        std::cout << "               #" << (k + 1) << ": C\"(" << nsym[k][1] << ")" << "\n";
                    } else {
                        std::cout << "               #" << (k + 1) << ": C(" << nsym[k][1] << ")" << "\n";
                    }
                }
            }
        }
        
        // S(n^k)
        for (int k = nsg + ncr + 1; k < nsg + ncr + nsr + 1; ++k) {
            if (nsym[k][2] > 1) {
                std::cout << "               #" << (k + 1) << ": S(" << nsym[k][1] 
                          << "^" << nsym[k][2] << ")" << "\n";
            } else {
                if (nsym[k][4] > 0) {
                    std::cout << "               #" << (k + 1) << ": S(" << nsym[k][1] 
                              << ")   -- with " << nsym[k][4] << " unmoved atoms" << "\n";
                } else {
                    std::cout << "               #" << (k + 1) << ": S(" << nsym[k][1] << ")" << "\n";
                }
            }
        } 
}
    }
    
    // Note: The original Fortran code ends here, but typically would call 
    // additional functions for symmetry classification if needed
};




/**
 * @brief Complete character table data 
 *
 * This is the 14×322 character table 
 * table contains the character values (typically 1, -1, 2, etc.) for all
 * irreducible representations across all 57 supported point groups.
 *
 * @return Complete character table as 14×322 array
 */
// Converted from Fortran 2D array 'chtab'
// Dimensions: [322][14]
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-const-variable"
const std::array<std::array<double, 14>, 322> chtab = {{
//
// Character tables
//
//data chtab(:,1:96)/
//
// C1
//  A
    { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// Cs
//  A'
    { 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A"
    { 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// Ci
//  Ag
    { 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Au
    { 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C2
//  A
    { 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B
    { 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C3
//  A
    { 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C4
//  A
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B
    { 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C5
//  A
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 0.618, -1.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, -1.618, 0.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C6
//  A
    { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B
    { 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.0, -1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, -1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C7
//  A
    { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.247, -0.445, -1.802, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, -0.445, -1.802, 1.247, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3
    { 2.0, -1.802, 1.247, -0.445, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C8
//  A
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B
    { 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.414, 0.0, -1.414, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3
    { 2.0, -1.414, 0.0, 1.414, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D2
//  A
    { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B3
    { 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D3
//  A1
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D4
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D5
//  A1
    { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 0.618, -1.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, -1.618, 0.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D6
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.0, -1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, -1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D7
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.247, -0.445, -1.802, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, -0.445, -1.802, 1.247, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3
    { 2.0, -1.802, 1.247, -0.445, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D8
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.414, 0.0, -1.414, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3
    { 2.0, -1.414, 0.0, 1.414, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C2v
//  A1
    { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C3v
//  A1
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C4v
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C5v
//  A1
    { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 0.618, -1.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, -1.618, 0.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C6v
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.0, -1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, -1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C7v
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.247, -0.445, -1.802, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, -0.445, -1.802, 1.247, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3
    { 2.0, -1.802, 1.247, -0.445, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C8v
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.414, 0.0, -1.414, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3
    { 2.0, -1.414, 0.0, 1.414, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//data chtab(:,97:142) /
//
// C2h
//  Ag
    { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Bg
    { 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Au
    { 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Bu
    { 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C3h
//  A'
    { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A"
    { 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E'
    { 2.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E"
    { 2.0, -1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C4h
//  Ag
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Bg
    { 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Eg
    { 2.0, 0.0, -2.0, 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Au
    { 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Bu
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Eu
    { 2.0, 0.0, -2.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C5h
//  A'
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A"
    { 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1'
    { 2.0, 0.618, -1.618, 2.0, 0.618, -1.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1"
    { 2.0, 0.618, -1.618, -2.0, -0.618, 1.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2'
    { 2.0, -1.618, 0.618, 2.0, -1.618, 0.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2"
    { 2.0, -1.618, 0.618, -2.0, 1.618, -0.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C6h
//  Ag
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Bg
    { 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1g
    { 2.0, 1.0, -1.0, -2.0, 2.0, -1.0, 1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2g
    { 2.0, -1.0, -1.0, 2.0, 2.0, -1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Au
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Bu
    { 1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1u
    { 2.0, 1.0, -1.0, -2.0, -2.0, 1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2u
    { 2.0, -1.0, -1.0, 2.0, -2.0, 1.0, 1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C7h
//  A'
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A"
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1'
    { 2.0, 1.247, -0.445, -1.802, 2.0, 1.247, -1.802, -0.445, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1"
    { 2.0, 1.247, -0.445, -1.802, -2.0, -1.247, 1.802, 0.445, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2'
    { 2.0, -0.445, -1.802, 1.247, 2.0, -0.445, 1.247, -1.802, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2"
    { 2.0, -0.445, -1.802, 1.247, -2.0, 0.445, -1.247, 1.802, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3'
    { 2.0, -1.802, 1.247, -0.445, 2.0, -1.802, -0.445, 1.247, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3"
    { 2.0, -1.802, 1.247, -0.445, -2.0, 1.802, 0.445, -1.247, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// C8h
//  Ag
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  Bg
    { 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  E1g
    { 2.0, 1.414, 0.0, -1.414, -2.0, 2.0, -1.414, 0.0, 1.414, -2.0, 0.0, 0.0, 0.0, 0.0 },
//  E2g
    { 2.0, 0.0, -2.0, 0.0, 2.0, 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0 },
//  E3g
    { 2.0, -1.414, 0.0, 1.414, -2.0, 2.0, 1.414, 0.0, -1.414, -2.0, 0.0, 0.0, 0.0, 0.0 },
//  Au
    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  Bu
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  E1u
    { 2.0, 1.414, 0.0, -1.414, -2.0, -2.0, 1.414, 0.0, -1.414, 2.0, 0.0, 0.0, 0.0, 0.0 },
//  E2u
    { 2.0, 0.0, -2.0, 0.0, 2.0, -2.0, 0.0, 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0 },
//  E3u
    { 2.0, -1.414, 0.0, 1.414, -2.0, -2.0, -1.414, 0.0, 1.414, 2.0, 0.0, 0.0, 0.0, 0.0 },
//data chtab(:,143:210) /
//
// D2h
//  Ag
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1g
    { 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2g
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B3g
    { 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Au
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1u
    { 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2u
    { 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B3u
    { 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D3h
//  A1'
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A1"
    { 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2'
    { 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2"
    { 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E'
    { 2.0, -1.0, 0.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E"
    { 2.0, -1.0, 0.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D4h
//  A1g
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  A2g
    { 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  B1g
    { 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  B2g
    { 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  Eg
    { 2.0, 0.0, -2.0, 0.0, 0.0, 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A1u
    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  A2u
    { 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  B1u
    { 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  B2u
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  Eu
    { 2.0, 0.0, -2.0, 0.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D5h
//  A1'
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A1"
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2'
    { 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2"
    { 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1'
    { 2.0, 0.618, -1.618, 0.0, 2.0, 0.618, -1.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1"
    { 2.0, 0.618, -1.618, 0.0, -2.0, -0.618, 1.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2'
    { 2.0, -1.618, 0.618, 0.0, 2.0, -1.618, 0.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2"
    { 2.0, -1.618, 0.618, 0.0, -2.0, 1.618, -0.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D6h
//  A1g
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0 },
//  A2g
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0 },
//  B1g
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 0.0, 0.0 },
//  B2g
    { 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 0.0, 0.0 },
//  E1g
    { 2.0, 1.0, -1.0, -2.0, 0.0, 0.0, 2.0, -1.0, 1.0, -2.0, 0.0, 0.0, 0.0, 0.0 },
//  E2g
    { 2.0, -1.0, -1.0, 2.0, 0.0, 0.0, 2.0, -1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0 },
//  A1u
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0 },
//  A2u
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0 },
//  B1u
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, 0.0, 0.0 },
//  B2u
    { 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 0.0, 0.0 },
//  E1u
    { 2.0, 1.0, -1.0, -2.0, 0.0, 0.0, -2.0, 1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0 },
//  E2u
    { 2.0, -1.0, -1.0, 2.0, 0.0, 0.0, -2.0, 1.0, 1.0, -2.0, 0.0, 0.0, 0.0, 0.0 },
//
// D7h
//  A1'
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  A1"
    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  A2'
    { 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  A2"
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  E1'
    { 2.0, 1.247, -0.445, -1.802, 0.0, 2.0, 1.247, -1.802, -0.445, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1"
    { 2.0, 1.247, -0.445, -1.802, 0.0, -2.0, -1.247, 1.802, 0.445, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2'
    { 2.0, -0.445, -1.802, 1.247, 0.0, 2.0, -0.445, 1.247, -1.802, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2"
    { 2.0, -0.445, -1.802, 1.247, 0.0, -2.0, 0.445, -1.247, 1.802, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3'
    { 2.0, -1.802, 1.247, -0.445, 0.0, 2.0, -1.802, -0.445, 1.247, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3"
    { 2.0, -1.802, 1.247, -0.445, 0.0, -2.0, 1.802, 0.445, -1.247, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D8h
//  A1g

    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },
//  A2g

    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0 },
//  B1g

    { 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0 },
//  B2g

    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0 },
//  E1g

    { 2.0, 1.414, 0.0, -1.414, -2.0, 0.0, 0.0, 2.0, -1.414, 0.0, 1.414, -2.0, 0.0, 0.0 },
//  E2g

    { 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0 },
//  E3g

    { 2.0, -1.414, 0.0, 1.414, -2.0, 0.0, 0.0, 2.0, 1.414, 0.0, -1.414, -2.0, 0.0, 0.0 },
//  A1u

    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
//  A2u

    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0 },
//  B1u

    { 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0 },
//  B2u

    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 },
//  E1u

    { 2.0, 1.414, 0.0, -1.414, -2.0, 0.0, 0.0, -2.0, 1.414, 0.0, -1.414, 2.0, 0.0, 0.0 },
//  E2u

    { 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, -2.0, 0.0, 2.0, 0.0, -2.0, 0.0, 0.0 },
//  E3u
    { 2.0, -1.414, 0.0, 1.414, -2.0, 0.0, 0.0, -2.0, -1.414, 0.0, 1.414, 2.0, 0.0, 0.0 },
//
//data chtab(:,211:307) /
//
// D2d
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D3d
//  A1g
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2g
    { 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Eg
    { 2.0, -1.0, 0.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A1u
    { 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2u
    { 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Eu
    { 2.0, -1.0, 0.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D4d
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.414, 0.0, -1.414, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3
    { 2.0, -1.414, 0.0, 1.414, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D5d
//  A1g
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2g
    { 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1g
    { 2.0, 0.618, -1.618, 0.0, 2.0, -1.618, 0.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2g
    { 2.0, -1.618, 0.618, 0.0, 2.0, 0.618, -1.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A1u
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2u
    { 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1u
    { 2.0, 0.618, -1.618, 0.0, -2.0, 1.618, -0.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2u
    { 2.0, -1.618, 0.618, 0.0, -2.0, -0.618, 1.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D6d
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.732, 1.0, 0.0, -1.0, -1.732, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, 1.0, -1.0, -2.0, -1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3
    { 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E4
    { 2.0, -1.0, -1.0, 2.0, -1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E5
    { 2.0, -1.732, 1.0, 0.0, -1.0, 1.732, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D7d
//  A1g
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  A2g
    { 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  E1g
    { 2.0, 1.247, -0.445, -1.802, 0.0, 2.0, -1.802, -0.445, 1.247, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2g
    { 2.0, -0.445, -1.802, 1.247, 0.0, 2.0, 1.247, -1.802, -0.445, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3g
    { 2.0, -1.802, 1.247, -0.445, 0.0, 2.0, -0.445, 1.247, -1.802, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A1u
    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  A2u
    { 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  E1u
    { 2.0, 1.247, -0.445, -1.802, 0.0, -2.0, 1.802, 0.445, -1.247, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2u
    { 2.0, -0.445, -1.802, 1.247, 0.0, -2.0, -1.247, 1.802, 0.445, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3u
    { 2.0, -1.802, 1.247, -0.445, 0.0, -2.0, 0.445, -1.247, 1.802, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// D8d
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0 },
//  B1
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0 },
//  B2
    { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.848, 1.414, 0.765, 0.0, -0.765, -1.414, -1.848, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, 1.414, 0.0, -1.414, -2.0, -1.414, 0.0, 1.414, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3
    { 2.0, 0.765, -1.414, -1.848, 0.0, 1.848, 1.414, -0.765, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E4
    { 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E5
    { 2.0, -0.765, -1.414, 1.848, 0.0, -1.848, 1.414, 0.765, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E6
    { 2.0, -1.414, 0.0, 1.414, -2.0, 1.414, 0.0, -1.414, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E7
    { 2.0, -1.848, 1.414, -0.765, 0.0, 0.765, -1.414, 1.848, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// S4
//  A
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B
    { 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// S6
//  Ag
    { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Eg
    { 2.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Au
    { 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Eu
    { 2.0, -1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// S8
//  A
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  B
    { 1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E1
    { 2.0, 1.414, 0.0, -1.414, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E2
    { 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E3
    { 2.0, -1.414, 0.0, 1.414, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// T
//  A
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  T
    { 3.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// Th
//  Ag
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Eg
    { 2.0, -1.0, 2.0, 2.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Tg
    { 3.0, 0.0, -1.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Au
    { 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Eu
    { 2.0, -1.0, 2.0, -2.0, 1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Tu
    { 3.0, 0.0, -1.0, -3.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// Td
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  T1
    { 3.0, 0.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  T2
    { 3.0, 0.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// O
//  A1
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  A2
    { 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  E
    { 2.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  T1
    { 3.0, 0.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  T2
    { 3.0, 0.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// Oh
//  A1g
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  A2g
    { 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  Eg
    { 2.0, -1.0, 2.0, 0.0, 0.0, 2.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  T1g
    { 3.0, 0.0, -1.0, 1.0, -1.0, 3.0, 0.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  T2g
    { 3.0, 0.0, -1.0, -1.0, 1.0, 3.0, 0.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  A1u
    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  A2u
    { 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  Eu
    { 2.0, -1.0, 2.0, 0.0, 0.0, -2.0, 1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  T1u
    { 3.0, 0.0, -1.0, 1.0, -1.0, -3.0, 0.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  T2u
    { 3.0, 0.0, -1.0, -1.0, 1.0, -3.0, 0.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//
//data chtab(:,308:) /
//
// I
//  A
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  T1
    { 3.0, 1.618, -0.618, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  T2
    { 3.0, -0.618, 1.618, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  G
    { 4.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  H
    { 5.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//
// Ih
//  Ag
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  T1g
    { 3.0, 1.618, -0.618, 0.0, -1.0, 3.0, -0.618, 1.618, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  T2g
    { 3.0, -0.618, 1.618, 0.0, -1.0, 3.0, 1.618, -0.618, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  Gg
    { 4.0, -1.0, -1.0, 1.0, 0.0, 4.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
//  Hg
    { 5.0, 0.0, 0.0, -1.0, 1.0, 5.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  Au
    { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
//  T1u
    { 3.0, 1.618, -0.618, 0.0, -1.0, -3.0, 0.618, -1.618, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  T2u
    { 3.0, -0.618, 1.618, 0.0, -1.0, -3.0, -1.618, 0.618, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
//  Gu
    { 4.0, -1.0, -1.0, 1.0, 0.0, -4.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
  //  Hu
    { 5.0, 0.0, 0.0, -1.0, 1.0, -5.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0 },
}};
#pragma clang diagnostic pop
//
///**
// * @brief Get the character table for a specific point group
// * @param point_group_idx Point group index (0-56)
// * @return Character table as vector of vectors (irreps x operations)
// */
//std::vector<std::vector<double>> getCharacterTable(int point_group_idx) {
//    if (point_group_idx < 0 || point_group_idx >= 57) {
//        return {};
//    }
//    int start = nir_[point_group_idx][0] - 1; // 0-based
//    int count = nir_[point_group_idx][1];
//    std::vector<std::vector<double>> table(count, std::vector<double>(14));
//    for (int i = 0; i < count; ++i) {
//        for (int j = 0; j < 14; ++j) {
//            table[i][j] = chtab[start + i][j];
//        }
//    }
//    return table;
//}
//
///**
// * @brief Get a specific character value
// * @param point_group_idx Point group index (0-56)
// * @param irrep_idx Irrep index within the point group (0-based)
// * @param operation_idx Operation index (0-13)
// * @return Character value
// */
//double getCharacter(int point_group_idx, int irrep_idx, int operation_idx) {
//    if (point_group_idx < 0 || point_group_idx >= 57 || operation_idx < 0 || operation_idx >= 14) {
//        return 0.0;
//    }
//    int start = nir_[point_group_idx][0] - 1;
//    int count = nir_[point_group_idx][1];
//    if (irrep_idx < 0 || irrep_idx >= count) {
//        return 0.0;
//    }
//    return chtab[start + irrep_idx][operation_idx];
//}




} // namespace symmetry