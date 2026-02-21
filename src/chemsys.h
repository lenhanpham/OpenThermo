/**
 * @file chemsys.h
 * @brief Core definitions and data structures for OpenThermo molecular thermochemistry calculations
 * @author Le Nhan Pham
 * @date 2025
 *
 * This header file contains fundamental constants, data structures, and type definitions
 * used throughout the OpenThermo program for molecular thermochemistry calculations.
 * It provides the basic building blocks for representing molecular systems and their properties.
 */

#ifndef CHEMSYS_H
#define CHEMSYS_H

#include <array>
#include <string>
#include <vector>
#include <cstdint>

/**
 * @brief Structure representing an atom in a molecular system
 *
 * This structure encapsulates the essential properties of an atom including its
 * position in 3D space, atomic number, and mass. It serves as the fundamental
 * building block for molecular representations in OpenThermo.
 */
struct Atom
{
    int    index;   /**< Atomic number (index in periodic table) */
    double x, y, z; /**< Cartesian coordinates in Angstrom units */
    double mass;    /**< Relative atomic mass in atomic mass units (amu) */
};

/**
 * @defgroup PhysicalConstants Physical Constants
 * @brief Fundamental physical constants used in thermochemistry
 * calculations
 *
 * These constants are used throughout OpenThermo for unit conversions and
 * thermodynamic
 * calculations. All values are in SI units unless otherwise specified.
 * @{
 */

/** Ideal gas constant in J/(mol·K) */
const double R = 8.3144648;

/** Boltzmann constant in J/K */
const double kb = 1.3806503e-23;

/** Avogadro constant (particles/mol) */
const double NA = 6.02214179e23;

/** Conversion factor: Hartree to electron volts */
const double au2eV = 27.2113838;

/** Conversion factor: Hartree to kcal/mol */
const double au2kcal_mol = 627.51;

/** Conversion factor: Hartree to kJ/mol */
const double au2kJ_mol = 2625.5;

/** Conversion factor: Hartree to Joules */
const double au2J = 4.359744575e-18;

/** Conversion factor: Hartree to cm⁻¹ (wavenumbers) */
const double au2cm_1 = 219474.6363;

/** Conversion factor: calories to Joules */
const double cal2J = 4.184;

/** Conversion factor: cm⁻¹ to Hz (frequency) */
const double wave2freq = 2.99792458e10;

/** Planck constant in J·s */
const double h = 6.62606896e-34;

/** Conversion factor: atomic mass units to kg */
const double amu2kg = 1.66053878e-27;

/** Mathematical constant π */
const double pi = 3.141592653589793;

/** Conversion factor: Bohr radius to Angstrom */
const double b2a = 0.52917720859;

/** Conversion factor: atmospheres to Pascals */
const double atm2Pa = 101325;

/** Maximum number of supported chemical elements */
const int nelesupp = 150;

/** Maximum number of isotopes per element */
const int maxiso = 300;

/** @} */

// Element names (index 0 is Bq, 1-150 are elements, 121-150 are "??")
const std::array<std::string, nelesupp + 1> ind2name{
    "Bq", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",  // 0-10
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",                    // 11-18
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr",  // 19-36
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
    "Xe",                                                                                                  // 37-54
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",  // 55-71
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",              // 72-86
    "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",  // 87-103
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "Un", "Ux",  // 104-120
    "??", "??", "??", "??", "??", "??", "??", "??", "??", "??",                                            // 121-130
    "??", "??", "??", "??", "??", "??", "??", "??", "??", "??",                                            // 131-140
    "??", "??", "??", "??", "??", "??", "??", "??", "??", "??"                                             // 141-150
};

/**
 * @brief Enumeration for low vibrational frequency treatment methods
 */
enum class LowVibTreatment : std::uint8_t
{
    Harmonic    = 0, /**< Standard rigid rotor harmonic oscillator */
    Truhlar     = 1, /**< Quasi-rigid rotor harmonic oscillator (raises low frequencies) */
    Grimme      = 2, /**< Grimme's interpolation between harmonic and free rotor */
    Minenkov    = 3, /**< Minenkov's interpolation scheme */
    HeadGordon  = 4  /**< Head-Gordon's interpolation for energy (+ optional entropy) */
};

/**
 * @brief Enumeration for the average moment of inertia (Bav) preset
 *        used in the free-rotor entropy term of quasi-RRHO methods.
 */
enum class BavPreset : std::uint8_t
{
    QChem  = 0, /**< B_av = 1 cm^-1 => I_av = 2.79928e-46 kg m^2 (Q-Chem manual) */
    Grimme = 1  /**< mu'_av = 1e-44 kg m^2 (Grimme 2012, ORCA/xtb/GoodVibes/Shermo) */
};

/// Return the numeric Bav value (kg m^2) for a given preset.
inline double bavPresetValue(BavPreset p)
{
    switch (p)
    {
        case BavPreset::QChem:  return 2.79928e-46;
        case BavPreset::Grimme: return 1e-44;
        default:                return 2.79928e-46;
    }
}

/// Return a human-readable label for a BavPreset.
inline const char* bavPresetName(BavPreset p)
{
    switch (p)
    {
        case BavPreset::QChem:  return "qchem";
        case BavPreset::Grimme: return "grimme";
        default:                return "qchem";
    }
}

/**
 * @brief Hardware and OpenMP parallelisation configuration
 *
 * Encapsulates all fields related to thread detection, validation, and
 * strategy selection. Separated from SystemData because execution
 * configuration is orthogonal to molecular chemistry: thermochemical
 * functions should not need to know about thread counts or HPC schedulers.
 */
struct ExecutionConfig
{
    int  omp_threads_requested   = 0;      ///< 0 = auto (default), >0 = user-specified
    int  omp_threads_actual      = 0;      ///< Final thread count after validation
    int  physical_cores_detected = 0;      ///< Hardware info for notifications
    int  scheduler_cpus_detected = 0;      ///< HPC scheduler allocated CPUs (0 = no scheduler)
    bool omp_user_override       = false;  ///< Whether user explicitly set thread count via -omp-threads
    int  omp_strategy            = 0;      ///< 0=outer (T/P scan), 1=inner (vibrational loop)
};

// Structure to hold system data
struct SystemData
{
    // Molecular information
    int                                  ncenter = 0;                 // Number of atoms
    std::vector<Atom>                    a;                           // Array of atoms
    int                                  ilinear = 0;                 // 0/1: Not linear / linear molecule
    int                                  nfreq   = 0;                 // Number of frequencies
    std::vector<double>                  freq;                        // Frequencies (s^-1, Hz)
    std::vector<double>                  wavenum;                     // Wavenumbers (cm^-1)
    std::array<std::array<double, 3>, 3> inertmat = {{{0.0}}};        // Moments of inertia matrix (amu*Bohr^2)
    std::array<double, 3>                inert    = {0.0, 0.0, 0.0};  // Inertia of principal axes (amu*Bohr^2)
    double                               E        = 0.0;              // Electronic energy (a.u.)
    double                               thermG;  // Thermal correction to Gibbs energy (kJ/mol, set by calcthermo)
    double                               totmass  = 0.0;  // Total mass (amu)
    int                                  spinmult = 0;    // Spin multiplicity (0 for .otm input)
    int                                  rotsym   = 0;    // Rotational symmetry number
    int                                  nelevel  = 0;    // Number of electronic excitation levels
    std::vector<double>                  elevel;          // Electronic exication energy of every considered level (eV)
    std::vector<int>                     edegen;          // Degeneracy of electronic energy levels

    // Parameters loaded from settings.ini or arguments
    std::string     PGnameinit     = "?";                      // Initial point group label
    std::string     PGname         = "?";                      // Point group label used in OpenThermo
    std::string     concstr         = "0";                      // Concentration string
    int             prtlevel        = 1;                        // Output verbosity: 0=minimal, 1=default, 2=verbose, 3=full
    int             prtvib          = 0;                        // Print vibration contributions
    LowVibTreatment lowVibTreatment = LowVibTreatment::Grimme;  // Low frequency treatment
    int             massmod         = 3;                        // Mass assignment mode
    int             outotm          = 0;                        // Output .otm file flag
    int             ipmode          = 0;       // 0 for isolated systems, and 1 for periodic solid state systems
    double          T               = 298.15;  // Temperature (K)
    double          Tlow = 0.0, Thigh = 0.0, Tstep = 0.0;  // Temperature range
    double          P    = 1.0;                            // Pressure (atm)
    double          Plow = 0.0, Phigh = 0.0, Pstep = 0.0;  // Pressure range
    double          sclZPE   = 1.0;                        // ZPE scaling factor
    double          sclheat  = 1.0;                        // Heat capacity scaling factor
    double          sclS     = 1.0;                        // Entropy scaling factor
    double          sclCV    = 1.0;                        // CV scaling factor
    double          ravib    = 100.0;                      // Vibrational averaging parameter
    double          intpvib  = 100.0;                      // Vibrational interpolation parameter
    bool            hgEntropy = true;                      // Enable entropy interpolation for Head-Gordon method
    BavPreset       bavPreset = BavPreset::Grimme;          // Bav preset for free-rotor entropy
    double          Bav       = 1e-44;                     // Average moment of inertia (kg m^2)
    bool            bavUserOverride = false;                // Whether user explicitly set -bav
    double          imagreal = 0.0;                        // Imaginary frequency threshold
    double          Eexter   = 0.0;                        // External electronic energy
    int vasp_energy_select   = 0;  // VASP energy selection: 0=energy  without entropy (default), 1=energy(sigma->0)

    // OpenMP / execution configuration
    ExecutionConfig exec;  ///< Hardware and OpenMP parallelisation settings

    // Special
    int inoset = 0;  // Skip settings.ini if 1

    // Others
    bool        alive = false;  // File existence flag
    std::string inputfile;      // Input file path

#ifdef _WIN32
    int isys = 1;  // Windows
#else
    int isys = 2;  // Linux/MacOS
#endif
};

// Arrays (to be initialized elsewhere, e.g., initmass)
// Note: Arrays are sized +1 to accommodate 1-based indexing used in the code
extern std::array<std::array<double, maxiso + 1>, nelesupp + 1> isomass;  // Isotope masses
extern std::array<std::array<double, maxiso + 1>, nelesupp + 1> isowei;   // Isotope composition
extern std::array<double, nelesupp + 1>                         elemass;  // Element masses

#endif  // CHEMSYS_H