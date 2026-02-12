/**
 * @file calc.h
 * @brief Header for thermochemistry calculation functions
 * @author Le Nhan Pham
 * @date 2025
 *
 * This header file declares functions for performing various thermochemistry
 * calculations including ensemble averaging, moment of inertia calculations,
 * thermodynamic property calculations, and vibrational contributions.
 */

#ifndef CALC_H
#define CALC_H

#include <vector>
#include <string>
#include <iostream>
#include "chemsys.h" // For SystemData

/**
 * @brief Namespace containing thermochemistry calculation functions
 *
 * This namespace encapsulates all functions related to molecular thermochemistry
 * calculations, including statistical mechanics, thermodynamic properties,
 * and ensemble averaging.
 */
namespace calc {

/**
 * @brief Per-contribution breakdown of thermodynamic properties
 *
 * Stores individual contributions (translational, rotational, vibrational,
 * electronic) as well as totals, avoiding the need for separate computation
 * in display functions.
 */
struct ThermoResult {
    // Translation
    double q_trans  = 1.0;
    double U_trans  = 0.0;
    double H_trans  = 0.0;
    double CV_trans = 0.0;
    double CP_trans = 0.0;
    double S_trans  = 0.0;

    // Rotation
    double q_rot  = 1.0;
    double U_rot  = 0.0;
    double CV_rot = 0.0;
    double S_rot  = 0.0;

    // Vibration
    double qvib_v0  = 1.0;
    double qvib_bot = 1.0;
    double U_vib      = 0.0;
    double U_vib_heat = 0.0;
    double ZPE        = 0.0;
    double CV_vib     = 0.0;
    double S_vib      = 0.0;

    // Electronic
    double q_ele  = 1.0;
    double U_ele  = 0.0;
    double CV_ele = 0.0;
    double S_ele  = 0.0;

    // Totals
    double corrU  = 0.0;
    double corrH  = 0.0;
    double corrG  = 0.0;
    double S_tot  = 0.0;
    double CV_tot = 0.0;
    double CP_tot = 0.0;
    double QV     = 0.0;   // total q(V=0) product
    double Qbot   = 0.0;   // total q(bot) product
};

/**
 * @brief Check if a file exists on the filesystem
 *
 * @param filename Path to the file to check
 * @return true if file exists and is accessible, false otherwise
 */
auto file_exists(const std::string& filename) -> bool;

/**
 * @brief Perform ensemble averaging calculations across multiple input files
 *
 * This function processes multiple molecular geometry files and computes
 * averaged thermodynamic properties using Boltzmann weighting.
 *
 * @param sys SystemData structure containing calculation parameters
 * @param filelist Vector of input file paths
 * @param Elist Vector to store electronic energies for each file
 * @param Ulist Vector to store internal energies for each file
 * @param Hlist Vector to store enthalpies for each file
 * @param Glist Vector to store Gibbs energies for each file
 * @param Slist Vector to store entropies for each file
 * @param CVlist Vector to store constant volume heat capacities for each file
 * @param CPlist Vector to store constant pressure heat capacities for each file
 * @param QVlist Vector to store vibrational partition functions for each file
 * @param Qbotlist Vector to store bottom partition functions for each file
 */
void ensemble(SystemData& sys, const std::vector<std::string>& filelist, std::vector<double>& Elist,
              std::vector<double>& Ulist, std::vector<double>& Hlist, std::vector<double>& Glist,
              std::vector<double>& Slist, std::vector<double>& CVlist, std::vector<double>& CPlist,
              std::vector<double>& QVlist, std::vector<double>& Qbotlist);

/**
 * @brief Calculate moments of inertia for the molecular system
 *
 * Computes the moment of inertia tensor and principal moments of inertia
 * based on the molecular geometry and atomic masses.
 *
 * @param sys SystemData structure containing molecular geometry and masses
 * @note Updates sys.inertmat and sys.inert arrays
 */
void calcinertia(SystemData& sys);

/**
 * @brief Calculate thermodynamic properties at given temperature and pressure
 *
 * This is the single source of truth for all thermodynamic computation.
 * Returns a ThermoResult with per-contribution breakdown and totals.
 *
 * @param sys SystemData structure with molecular data and parameters
 * @param T Temperature in Kelvin
 * @param P Pressure in atmospheres
 * @return ThermoResult with per-contribution breakdown and totals
 */
ThermoResult calcthermo(const SystemData& sys, double T, double P);

/**
 * @brief Calculate thermodynamic properties (output-parameter overload)
 *
 * Convenience overload that extracts totals from ThermoResult into
 * individual output parameters. Used by the T/P scan loop.
 *
 * @param sys SystemData structure with molecular data and parameters
 * @param T Temperature in Kelvin
 * @param P Pressure in atmospheres
 * @param corrU [out] Thermal correction to internal energy (kJ/mol)
 * @param corrH [out] Thermal correction to enthalpy (kJ/mol)
 * @param corrG [out] Thermal correction to Gibbs energy (kJ/mol)
 * @param S [out] Total entropy (J/mol-K)
 * @param CV [out] Constant volume heat capacity (J/mol-K)
 * @param CP [out] Constant pressure heat capacity (J/mol-K)
 * @param QV [out] Vibrational partition function
 * @param Qbot [out] Bottom partition function (rotational + electronic)
 */
void calcthermo(const SystemData& sys, double T, double P, double& corrU, double& corrH, double& corrG,
                double& S, double& CV, double& CP, double& QV, double& Qbot);

/**
 * @brief Display detailed thermodynamic properties for a single (T, P) point
 *
 * Calls calcthermo() internally and formats the per-contribution breakdown
 * to the console. Sets sys.thermG as a side effect.
 *
 * @param sys SystemData structure containing the results
 */
void showthermo(SystemData& sys);

/**
 * @brief Display per-mode vibrational detail table
 *
 * Prints per-mode partition functions and thermodynamic contributions
 * (ZPE, U, CV, S) for each vibrational mode. Separated from showthermo()
 * to support future verbosity-level control.
 *
 * @param sys SystemData with vibrational data
 * @param T Temperature in Kelvin
 * @param out Output stream (cout or file stream for .vibcon export)
 */
void showvibdetail(const SystemData& sys, double T, std::ostream& out);

/**
 * @brief Calculate concentration-dependent Gibbs energy correction
 *
 * Computes the Gibbs energy correction due to concentration effects
 * in solution or gas-phase calculations.
 *
 * @param sys SystemData structure with system parameters
 * @param concnow [out] Current concentration
 * @param concspec [out] Species concentration
 * @param Gconc [out] Concentration correction to Gibbs energy
 */
void getGconc(const SystemData &sys, const double &concnow, const double &concspec, double &Gconc);

/**
 * @brief Calculate vibrational contributions for a specific mode
 *
 * Computes the thermodynamic contributions from a single vibrational mode
 * including zero-point energy, thermal energy, heat capacity, and entropy.
 *
 * @param sys SystemData structure with vibrational data
 * @param i Index of the vibrational mode (0-based)
 * @param T Temperature in Kelvin
 * @param tmpZPE [out] Zero-point energy contribution
 * @param tmpheat [out] Thermal energy contribution
 * @param tmpCV [out] Constant volume heat capacity contribution
 * @param tmpS [out] Entropy contribution
 */
void getvibcontri(const SystemData& sys, int i, double T, double& tmpZPE, double& tmpheat, double& tmpCV, double& tmpS);

/**
 * @brief Calculate total vibrational contributions at temperature
 *
 * Computes the total thermodynamic contributions from all vibrational modes
 * at a given temperature using the harmonic oscillator approximation.
 *
 * @param sys SystemData structure with vibrational data
 * @param T Temperature in Kelvin
 * @param U_vib [out] Total vibrational internal energy
 * @param CV_vib [out] Total vibrational heat capacity at constant volume
 * @param S_vib [out] Total vibrational entropy
 * @param QV [out] Vibrational partition function
 */
void getvibcontri(const SystemData& sys, double T, double& U_vib, double& CV_vib, double& S_vib, double& QV);

} // namespace calc

#endif // CALC_H