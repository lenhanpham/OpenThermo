
/**
 * @file calc.cpp
 * @brief Implementation of thermochemistry calculation functions
 * @author Le Nhan Pham
 * @date 2025
 *
 * This file contains the core implementation of molecular thermochemistry
 * calculations including statistical mechanics, thermodynamic property
 * calculations, moment of inertia computations, and ensemble averaging.
 */

#include "calc.h"
#include "chemsys.h"
#include "loadfile.h"
#include "symmetry.h"
#include "util.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace util;


#ifndef M_PI
    #define M_PI 3.141592653589793 /**< Define π if not already defined */
#endif

namespace calc
{

    // Forward declaration
    void elecontri(const SystemData& sys, double& tmpq, double& tmpheat, double& tmpCV, double& tmpS);


    /**
     * @brief Check if a file exists and is accessible
     *
     * @param filename Path to the file to check

     * * @return true if file exists and can be opened, false otherwise
     */
    bool file_exists(const std::string& filename)
    {
        std::ifstream f(filename);
        return f.good();
    }

    /**
     * @brief Extract basename without extension from a file path
     * @param filepath Full file path
     *
     * @return Basename without extension
     */
    std::string get_basename_without_extension(const std::string& filepath)
    {
        // Find last directory separator
        size_t      last_slash = filepath.find_last_of("/\\");
        std::string filename   = (last_slash != std::string::npos) ? filepath.substr(last_slash + 1) : filepath;

        // Find last dot for extension
        size_t last_dot = filename.rfind('.');
        return (last_dot != std::string::npos) ? filename.substr(0, last_dot) : filename;
    }


    /**
     * @brief Perform ensemble averaging calculations across multiple input files
     *
     * This function processes multiple molecular geometry files and computes
     * averaged thermodynamic properties using Boltzmann weighting. It loads
     * each file, performs thermochemistry calculations, and accumulates results
     * for ensemble averaging.
     *
     * @param sys SystemData structure containing calculation parameters
     * @param filelist Vector of input file paths to process
     * @param Elist [out] Electronic energies for each file
     * @param Ulist [out] Internal energies for each file
     * @param Hlist [out] Enthalpies for each file
     * @param Glist [out] Gibbs energies for each file
     * @param Slist [out] Entropies for each file
     * @param CVlist [out] Constant volume heat capacities for each file
     * @param CPlist [out] Constant pressure heat capacities for each file
     * @param QVlist [out] Vibrational partition functions for each file
     * @param Qbotlist [out] Bottom partition functions for each file
     *
     * @note Supports multiple file formats: .otm, Gaussian, ORCA, GAMESS, NWChem, CP2K, xTB
     * @note Automatically detects file format and loads appropriate parser
     */
    void ensemble(SystemData&                     sys,
                  const std::vector<std::string>& filelist,
                  std::vector<double>&            Elist,
                  std::vector<double>&            Ulist,
                  std::vector<double>&            Hlist,
                  std::vector<double>&            Glist,
                  std::vector<double>&            Slist,
                  std::vector<double>&            CVlist,
                  std::vector<double>&            CPlist,
                  std::vector<double>&            QVlist,
                  std::vector<double>&            Qbotlist)
    {
        int                 nfile = filelist.size();
        std::vector<double> wei(nfile);
        for (int ifile = 0; ifile < nfile; ++ifile)
        {
            sys.inputfile = filelist[ifile];
            if (!file_exists(sys.inputfile))
            {
                std::cerr << "Error: Unable to find " << sys.inputfile << "\n";
                std::cerr << "Press ENTER button to exit program" << "\n";
                std::cin.get();
                std::exit(1);
            }

            std::cout << "Processing " << sys.inputfile << "... (" << (ifile + 1) << " of " << nfile << " )" << "\n";

            size_t otm_pos = sys.inputfile.find(".otm");
            if (otm_pos != std::string::npos)
            {
                LoadFile::loadotm(sys);
            }
            else
            {
                auto pcprog = deterprog(sys);
                if (pcprog == QuantumChemistryProgram::Unknown)
                {
                    throw std::runtime_error("Invalid program type for file " + sys.inputfile);
                }
                if (pcprog == QuantumChemistryProgram::Gaussian)
                    LoadFile::loadgau(sys);
                else if (pcprog == QuantumChemistryProgram::Orca)
                    LoadFile::loadorca(sys);
                else if (pcprog == QuantumChemistryProgram::Gamess)
                    LoadFile::loadgms(sys);
                else if (pcprog == QuantumChemistryProgram::Nwchem)
                    LoadFile::loadnw(sys);
                else if (pcprog == QuantumChemistryProgram::Cp2k)
                    LoadFile::loadCP2K(sys);
                else if (pcprog == QuantumChemistryProgram::Xtb)
                    LoadFile::loadxtb(sys);
                else if (pcprog == QuantumChemistryProgram::Vasp)
                    LoadFile::loadvasp(sys);

                modmass(sys);
                sys.nelevel = 1;
                sys.elevel  = {0.0};
                int deg     = std::max(sys.spinmult, 1);
                sys.edegen  = {deg};
                if (sys.outotm == 1)
                    outotmfile(sys);
            }

            if (sys.Eexter != 0.0)
            {
                sys.E = sys.Eexter;
                std::cout << "Note: The electronic energy specified by Eexter will be used" << "\n";
            }
            else if (Elist[ifile] != 0.0)
            {
                sys.E = Elist[ifile];
            }
            else
            {
                Elist[ifile] = sys.E;
                if (sys.E != 0.0)
                {
                    std::cout << "Note: The electronic energy extracted from file will be used" << "\n";
                }
            }

            sys.totmass = 0.0;
            for (const auto& atom : sys.a)
            {
                sys.totmass += atom.mass;
            }

            calcinertia(sys);

            sys.ilinear = 0;
            for (double in : sys.inert)
            {
                if (in < 0.001)
                {
                    sys.ilinear = 1;
                    break;
                }
            }

            // Symmetry detection
            symmetry::SymmetryDetector symDetector;
            symDetector.ncenter = sys.a.size();
            symDetector.a       = sys.a;
            symDetector.a_index.resize(sys.a.size());
            for (size_t i = 0; i < sys.a.size(); ++i)
            {
                symDetector.a_index[i] = i;
            }
            symDetector.detectPG(sys.prtvib ? 1 : 0);
            sys.rotsym  = symDetector.rotsym;
            sys.PGname = symDetector.PGname;

            // Handle imaginary frequencies
            if (sys.imagreal != 0.0)
            {
                for (int j = 0; j < sys.nfreq; ++j)
                {
                    if (sys.wavenum[j] < 0 && std::abs(sys.wavenum[j]) < sys.imagreal)
                    {
                        sys.wavenum[j] = std::abs(sys.wavenum[j]);
                        std::cout << "Note: Imaginary frequency " << sys.wavenum[j] << " cm^-1 set to real frequency!"
                                  << "\n";
                    }
                }
            }

            std::vector<double> freq(sys.nfreq);
            for (int j = 0; j < sys.nfreq; ++j)
            {
                freq[j] = sys.wavenum[j] * wave2freq;
            }
            sys.freq = freq;

            double thermU, thermH, thermG, CP_tot, QV, Qbot;
            calcthermo(sys, sys.T, sys.P, thermU, thermH, thermG, Slist[ifile], CVlist[ifile], CP_tot, QV, Qbot);
            sys.thermG = thermG;

            Glist[ifile] = thermG / au2kJ_mol + sys.E;

            Ulist[ifile]    = thermU / au2kJ_mol + sys.E;
            Hlist[ifile]    = thermH / au2kJ_mol + sys.E;
            CPlist[ifile]   = CP_tot / au2kJ_mol;
            QVlist[ifile]   = QV / NA;
            Qbotlist[ifile] = Qbot / NA;

            // Deallocate
            sys.a.clear();
            sys.elevel.clear();
            sys.edegen.clear();
            sys.freq.clear();
            sys.wavenum.clear();
        }

        std::cout << "\n";
        double qall = 0.0;
        double Gmin = *std::min_element(Glist.begin(), Glist.end());
        std::cout << "#System       U               H               G             S          CV" << "\n";
        std::cout << "             a.u.            a.u.            a.u.        J/mol/K     J/mol/K" << "\n";
        for (int ifile = 0; ifile < nfile; ++ifile)
        {
            double dG = (Glist[ifile] - Gmin) * au2kJ_mol * 1000.0;  // Relative free energy in J/mol
            qall += std::exp(-dG / (R * sys.T));
            std::cout << std::fixed << std::setprecision(6) << std::setw(5) << (ifile + 1) << std::setw(16)
                      << Ulist[ifile] << std::setw(16) << Hlist[ifile] << std::setw(16) << Glist[ifile]
                      << std::setprecision(3) << std::setw(12) << Slist[ifile] << std::setw(12) << CVlist[ifile]
                      << "\n";
        }

        std::cout << "\n";
        for (int ifile = 0; ifile < nfile; ++ifile)
        {
            double dG  = (Glist[ifile] - Gmin) * au2kJ_mol * 1000.0;
            wei[ifile] = std::exp(-dG / (R * sys.T)) / qall;
            std::cout << " System" << std::setw(5) << (ifile + 1) << "     Relative G=" << std::fixed
                      << std::setprecision(3) << std::setw(9) << (Glist[ifile] - Gmin) * au2kJ_mol
                      << " kJ/mol     Boltzmann weight=" << std::setprecision(3) << std::setw(8) << wei[ifile] * 100.0
                      << " %" << "\n";
        }

        double weiE = 0.0, weiU = 0.0, weiH = 0.0, weiS = 0.0, weiCV = 0.0;
        for (int ifile = 0; ifile < nfile; ++ifile)
        {
            weiE += wei[ifile] * Elist[ifile];
            weiU += wei[ifile] * Ulist[ifile];
            weiH += wei[ifile] * Hlist[ifile];
            weiS += wei[ifile] * Slist[ifile];
            weiCV += wei[ifile] * CVlist[ifile];
        }
        double confS = 0.0;
        for (int ifile = 0; ifile < nfile; ++ifile)
        {
            if (wei[ifile] > 0.0)
                confS -= R * wei[ifile] * std::log(wei[ifile]);
        }
        weiS += confS;
        double weiG = weiH - sys.T * weiS / 1000.0 / au2kJ_mol;
        std::cout << "\n";
        std::cout << "Conformation weighted data:" << "\n";
        std::cout << " Electronic energy: " << std::fixed << std::setprecision(6) << std::setw(16) << weiE << " a.u."
                  << "\n";
        std::cout << " U: " << std::setw(16) << weiU << " a.u." << "\n";
        std::cout << " H: " << std::setw(16) << weiH << " a.u." << "\n";
        std::cout << " G: " << std::setw(16) << weiG << " a.u." << "\n";
        std::cout << " S: " << std::fixed << std::setprecision(3) << std::setw(13) << weiS
                  << " J/mol/K    Conformation entropy:" << std::setw(10) << confS << " J/mol/K" << "\n";
        std::cout << " CV:" << std::setw(13) << weiCV << " J/mol/K" << "\n";
        std::cout << " CP:" << std::setw(13) << weiCV + R << " J/mol/K" << "\n";

        if (sys.concstr != "0")
        {
            double concnow = 0.0, concspec = std::stod(sys.concstr), Gconc;
            getGconc(sys, concnow, concspec, Gconc);
            std::cout << "\n";
            std::cout << " Present concentration (estimated by ideal gas model):" << std::fixed << std::setprecision(6)
                      << std::setw(10) << concnow << " mol/L" << "\n";
            std::cout << " Concentration specified by \"conc\" parameter:" << std::setw(12) << concspec << " mol/L"
                      << "\n";
            std::cout << " delta-G of conc. change:" << std::fixed << std::setprecision(3) << std::setw(11) << Gconc
                      << " kJ/mol" << std::setw(11) << Gconc / cal2J << " kcal/mol" << std::setprecision(6)
                      << std::setw(11) << Gconc / au2kJ_mol << " a.u." << "\n";
            std::cout << " Weighted Gibbs free energy at specified concentration: " << std::fixed
                      << std::setprecision(7) << std::setw(19) << (weiG + Gconc / au2kJ_mol) << " a.u." << "\n";
        }
    }


    /**
     * @brief Calculate thermodynamic properties at given temperature and pressure
     *
     * This is the main thermodynamic calculation function that computes all
     * thermochemical properties using statistical mechanics. It calculates
     * translational, rotational, vibrational, and electronic contributions
     * to thermodynamic functions.
     *
     * @param sys SystemData structure with molecular data and parameters
     * @param T Temperature in Kelvin
     * @param P Pressure in atmospheres (currently unused)
     * @param corrU [out] Thermal correction to internal energy (kJ/mol)
     * @param corrH [out] Thermal correction to enthalpy (kJ/mol)
     * @param corrG [out] Thermal correction to Gibbs energy (kJ/mol)
     * @param S [out] Total entropy (J/mol·K)
     * @param CV [out] Constant volume heat capacity (J/mol·K)
     * @param CP [out] Constant pressure heat capacity (J/mol·K)
     * @param QV [out] Vibrational partition function
     * @param Qbot [out] Bottom partition function (rotational + electronic)
     *
     * @note Uses statistical mechanics formulas for ideal gas
     * @note Includes corrections for low-frequency vibrational modes
     * @note Accounts for molecular symmetry and electronic degeneracy
     */
    void calcthermo(SystemData& sys,
                    double      T,
                    double /*P*/,
                    double& corrU,
                    double& corrH,
                    double& corrG,
                    double& S,
                    double& CV,
                    double& CP,
                    double& QV,
                    double& Qbot)
    {
        double q_trans = 1.0, U_trans = 0.0, CV_trans = 0.0, CP_trans = 0.0, H_trans = 0.0, S_trans = 0.0;
        double q_rot = 1.0, U_rot = 0.0, CV_rot = 0.0, S_rot = 0.0;
        double qvib_v0 = 1.0, qvib_bot = 1.0, U_vib_heat = 0.0, CV_vib = 0.0, S_vib = 0.0, ZPE = 0.0;
        double q_ele = 1.0, U_ele = 0.0, CV_ele = 0.0, S_ele = 0.0;

        // Translation contribution
        if (sys.ipmode == 0)
        {
            double P_Pa = sys.P * atm2Pa;
            q_trans =
                std::pow(2.0 * M_PI * (sys.totmass * amu2kg) * kb * sys.T / (h * h), 3.0 / 2.0) * R * sys.T / P_Pa;
            CV_trans = 3.0 / 2.0 * R;
            CP_trans = 5.0 / 2.0 * R;
            U_trans  = 3.0 / 2.0 * R * sys.T / 1000.0;
            H_trans  = 5.0 / 2.0 * R * sys.T / 1000.0;
            S_trans  = R * (std::log(q_trans / NA) + 5.0 / 2.0);
        }
        else if (sys.ipmode == 1)
        {
            q_trans  = 1.0;
            CV_trans = 0.0;
            CP_trans = 0.0;
            U_trans  = 0.0;
            H_trans  = 0.0;
            S_trans  = 0.0;
        }

        // Rotation contribution
        if (sys.ipmode == 0)
        {
            double sum_inert = sys.inert[0] + sys.inert[1] + sys.inert[2];
            if (sum_inert < 1e-10)
            {  // Single atom
                q_rot  = 1.0;
                U_rot  = 0.0;
                CV_rot = 0.0;
                S_rot  = 0.0;
            }
            else
            {
                std::vector<double> inertkg(3);
                for (int i = 0; i < 3; ++i)
                {
                    inertkg[i] = sys.inert[i] * amu2kg * std::pow(b2a * 1e-10, 2);  // Convert to kg*m^2
                }
                if (sys.ilinear == 1)
                {  // Linear molecule
                    q_rot  = 8.0 * M_PI * M_PI * inertkg[2] * kb * sys.T / sys.rotsym / (h * h);
                    U_rot  = R * sys.T / 1000.0;
                    CV_rot = R;
                    S_rot  = R * (std::log(q_rot) + 1.0);
                }
                else
                {  // Non-linear molecule
                    q_rot = 8.0 * M_PI * M_PI / sys.rotsym / std::pow(h, 3) *
                            std::pow(2.0 * M_PI * kb * sys.T, 3.0 / 2.0) *
                            std::sqrt(inertkg[0] * inertkg[1] * inertkg[2]);
                    U_rot  = 3.0 * R * sys.T / 2.0 / 1000.0;
                    CV_rot = 3.0 * R / 2.0;
                    S_rot  = R * (std::log(q_rot) + 3.0 / 2.0);
                }
            }
        }
        else if (sys.ipmode == 1)
        {
            q_rot  = 1.0;
            U_rot  = 0.0;
            CV_rot = 0.0;
            S_rot  = 0.0;
        }

        // Vibration contribution
        for (int i = 0; i < sys.nfreq; ++i)
        {
            if (sys.freq[i] <= 0.0)
                continue;
            double freqtmp = sys.freq[i];
            if (sys.lowVibTreatment == LowVibTreatment::Truhlar && sys.wavenum[i] < sys.ravib)
            {
                freqtmp = sys.ravib * wave2freq;
            }
            double tmpv0  = 1.0 / (1.0 - std::exp(-h * freqtmp / (kb * sys.T)));
            double tmpbot = std::exp(-h * freqtmp / (kb * 2.0 * sys.T)) / (1.0 - std::exp(-h * freqtmp / (kb * sys.T)));
            qvib_v0 *= tmpv0;
            qvib_bot *= tmpbot;
        }

        for (int i = 0; i < sys.nfreq; ++i)
        {
            double tmpZPE, tmpheat, tmpCV, tmpS;
            getvibcontri(sys, i, tmpZPE, tmpheat, tmpCV, tmpS);
            ZPE += tmpZPE;
            U_vib_heat += tmpheat;
            CV_vib += tmpCV;
            S_vib += tmpS;
        }
        double U_vib = U_vib_heat + ZPE;

        // Electron contribution
        elecontri(sys, q_ele, U_ele, CV_ele, S_ele);

        // Total values
        double CV_tot = CV_trans + CV_rot + CV_vib + CV_ele;
        double CP_tot = CP_trans + CV_rot + CV_vib + CV_ele;
        double S_tot  = S_trans + S_rot + S_vib + S_ele;
        double thermU = U_trans + U_rot + U_vib + U_ele;
        double thermH = H_trans + U_rot + U_vib + U_ele;
        double thermG;
        if (T == 0.0)
        {
            thermG = thermH;
        }
        else
        {
            thermG = thermH - T * S_tot / 1000.0;
        }

        // Assign output parameters
        corrU = thermU;
        corrH = thermH;
        corrG = thermG;
        S     = S_tot;
        CV    = CV_tot;
        CP    = CP_tot;
        QV    = q_trans * q_rot * qvib_v0 * q_ele;
        Qbot  = q_trans * q_rot * qvib_bot * q_ele;
    }


    /**
     * @brief Display calculated thermodynamic properties in formatted output
     *
     * This function
     * computes and displays all thermochemical properties including
     * translational, rotational, vibrational, and
     * electronic contributions to
     * thermodynamic functions. The output includes partition functions, energies,

     * * entropies, heat capacities, and final corrected values.
     *
     * @param sys SystemData structure
     * containing molecular data and calculation parameters
     *
     * @note Outputs results to stdout with detailed
     * breakdown by contribution type
     * @note Includes vibrational mode analysis if prtvib is enabled
     * @note
     * Handles different low-frequency treatment methods (harmonic, Truhlar, Grimme)
     * @note Accounts for
     * concentration-dependent corrections if specified
     */
    void showthermo(SystemData& sys)
    {
        double q_trans = 1.0, U_trans = 0.0, CV_trans = 0.0, CP_trans = 0.0, H_trans = 0.0, S_trans = 0.0;
        double q_rot = 1.0, U_rot = 0.0, CV_rot = 0.0, S_rot = 0.0;
        double qvib_v0 = 1.0, qvib_bot = 1.0, U_vib_heat = 0.0, CV_vib = 0.0, S_vib = 0.0, ZPE = 0.0;
        double q_ele = 1.0, U_ele = 0.0, CV_ele = 0.0, S_ele = 0.0;
        double thermU, thermH, thermG, CV_tot, CP_tot, S_tot;

        // Translation contribution
        if (sys.ipmode == 0)
        {
            std::cout << "\nNote: Only for translation, U is different to H, and CV is different to CP\n"
                      << "\n";
            std::cout << "                        ------- Translation -------\n"
                      << "                        ---------------------------\n";
            double P_Pa = sys.P * atm2Pa;
            q_trans =
                std::pow(2.0 * M_PI * (sys.totmass * amu2kg) * kb * sys.T / (h * h), 3.0 / 2.0) * R * sys.T / P_Pa;
            CV_trans = 3.0 / 2.0 * R;
            CP_trans = 5.0 / 2.0 * R;
            U_trans  = 3.0 / 2.0 * R * sys.T / 1000.0;
            H_trans  = 5.0 / 2.0 * R * sys.T / 1000.0;
            S_trans  = R * (std::log(q_trans / NA) + 5.0 / 2.0);
            std::cout << std::scientific << std::setprecision(6) << " Translational q: " << std::setw(16) << q_trans
                      << "     q/NA: " << std::setw(16) << q_trans / NA << "\n";
            std::cout << std::fixed << std::setprecision(3) << " Translational U: " << std::setw(10) << U_trans
                      << " kJ/mol " << std::setw(10) << U_trans / cal2J << " kcal/mol\n";
            std::cout << " Translational H: " << std::setw(10) << H_trans << " kJ/mol " << std::setw(10)
                      << H_trans / cal2J << " kcal/mol\n";
            std::cout << " Translational S: " << std::setw(10) << S_trans << " J/mol/K" << std::setw(10)
                      << S_trans / cal2J << " cal/mol/K  -TS:" << std::setw(8) << -S_trans / cal2J / 1000.0 * sys.T
                      << " kcal/mol\n";
            std::cout << " Translational CV:" << std::setw(10) << CV_trans << " J/mol/K" << std::setw(10)
                      << CV_trans / cal2J << " cal/mol/K\n";
            std::cout << " Translational CP:" << std::setw(10) << CP_trans << " J/mol/K" << std::setw(10)
                      << CP_trans / cal2J << " cal/mol/K\n";
        }
        else if (sys.ipmode == 1)
        {
            std::cout << "\nTranslation contribution is ignored since ipmode=1\n";
            q_trans  = 1.0;
            CV_trans = 0.0;
            CP_trans = 0.0;
            U_trans  = 0.0;
            H_trans  = 0.0;
            S_trans  = 0.0;
        }

        // Rotation contribution
        if (sys.ipmode == 0)
        {
            std::cout << "\n                        -------- Rotation --------\n"
                      << "                        --------------------------\n";
            double sum_inert = sys.inert[0] + sys.inert[1] + sys.inert[2];
            if (sum_inert < 1e-10)
            {  // Single atom
                q_rot  = 1.0;
                U_rot  = 0.0;
                CV_rot = 0.0;
                S_rot  = 0.0;
            }
            else
            {
                std::vector<double> inertkg(3);
                for (int i = 0; i < 3; ++i)
                {
                    inertkg[i] = sys.inert[i] * amu2kg * std::pow(b2a * 1e-10, 2);  // Convert to kg*m^2
                }
                if (sys.ilinear == 1)
                {  // Linear molecule
                   // Use the LARGEST moment of inertia (not inertkg[2]!)
                    double largest_inertkg = *std::max_element(inertkg.begin(), inertkg.end());
                    q_rot                  = 8.0 * M_PI * M_PI * largest_inertkg * kb * sys.T / sys.rotsym / (h * h);
                    U_rot                  = R * sys.T / 1000.0;
                    CV_rot                 = R;
                    S_rot                  = R * (std::log(q_rot) + 1.0);
                }
                else
                {  // Non-linear molecule
                    q_rot = 8.0 * M_PI * M_PI / sys.rotsym / std::pow(h, 3) *
                            std::pow(2.0 * M_PI * kb * sys.T, 3.0 / 2.0) *
                            std::sqrt(inertkg[0] * inertkg[1] * inertkg[2]);
                    U_rot  = 3.0 * R * sys.T / 2.0 / 1000.0;
                    CV_rot = 3.0 * R / 2.0;
                    S_rot  = R * (std::log(q_rot) + 3.0 / 2.0);
                }
            }
            std::cout << std::scientific << std::setprecision(6) << " Rotational q: " << std::setw(16) << q_rot << "\n";
            std::cout << std::fixed << std::setprecision(3) << " Rotational U: " << std::setw(10) << U_rot << " kJ/mol "
                      << std::setw(10) << U_rot / cal2J << " kcal/mol    =H\n";
            std::cout << " Rotational S: " << std::setw(10) << S_rot << " J/mol/K" << std::setw(10) << S_rot / cal2J
                      << " cal/mol/K   -TS:" << std::setw(8) << -S_rot / cal2J / 1000.0 * sys.T << " kcal/mol\n";
            std::cout << " Rotational CV:" << std::setw(10) << CV_rot << " J/mol/K" << std::setw(10) << CV_rot / cal2J
                      << " cal/mol/K   =CP\n";
        }
        else if (sys.ipmode == 1)
        {
            std::cout << "\nRotation contribution is ignored since ipmode=1\n";
            q_rot  = 1.0;
            U_rot  = 0.0;
            CV_rot = 0.0;
            S_rot  = 0.0;
        }

        // Vibration contribution
        std::ofstream vibfile;
        std::ostream* ivibout = &std::cout;
        std::string   vibcon_filename;
        if (sys.prtvib == -1)
        {
            vibcon_filename = get_basename_without_extension(sys.inputfile) + ".vibcon";
            vibfile.open(vibcon_filename);
            if (!vibfile.is_open())
                throw std::runtime_error("showthermo: Unable to open " + vibcon_filename);
            ivibout = &vibfile;
        }

        std::cout << "\n                        -------- Vibration --------\n"
                  << "                        ---------------------------\n";
        if (sys.lowVibTreatment == LowVibTreatment::Truhlar)
        {
            int nlow = 0;
            for (double wn : sys.wavenum)
            {
                if (std::abs(wn) < sys.ravib)
                    ++nlow;
            }
            if (nlow > 0)
            {
                std::cout << "Note: " << nlow << " low frequencies are raised to " << std::fixed << std::setprecision(1)
                          << sys.ravib << " cm^-1 during calculating S, U(T)-U(0), CV and q\n\n";
            }
        }
        else if (sys.lowVibTreatment == LowVibTreatment::Grimme)
        {
            std::cout << "Note: Interpolation between harmonic oscillator model and free rotor model is \n"
                         "      used to evaluate S, other terms are identical to harmonic oscillator model\n\n";
        }
        else if (sys.lowVibTreatment == LowVibTreatment::Minenkov)
        {
            std::cout << "Note: Interpolation between harmonic oscillator model and free rotor model is \n"
                         "      used to evaluate S and U(T). "
                      << "In this case ZPE and U(T)-U(0) cannot be separated and thus not shown. \n"
                         "Other terms are identical to harmonic oscillator model\n\n";
        }

        // Calculate partition function
        if (std::abs(sys.prtvib) == 1)
        {
            if (sys.sclZPE != 1.0 || sys.sclheat != 1.0 || sys.sclS != 1.0 || sys.sclCV != 1.0)
            {
                *ivibout << "Note: The wavenumbers shown below are unscaled ones\n\n";
            }
            *ivibout << " Mode  Wavenumber    Freq        Vib. Temp.    q(V=0)        q(bot)\n";
            *ivibout << "         cm^-1        GHz            K\n";
        }

        for (int i = 0; i < sys.nfreq; ++i)
        {
            if (sys.freq[i] <= 0.0)
                continue;
            double freqtmp = sys.freq[i];
            if (sys.lowVibTreatment == LowVibTreatment::Truhlar && sys.wavenum[i] < sys.ravib)
            {
                freqtmp = sys.ravib * wave2freq;
            }
            double tmpv0  = 1.0 / (1.0 - std::exp(-h * freqtmp / (kb * sys.T)));
            double tmpbot = std::exp(-h * freqtmp / (kb * 2.0 * sys.T)) / (1.0 - std::exp(-h * freqtmp / (kb * sys.T)));
            qvib_v0 *= tmpv0;
            qvib_bot *= tmpbot;
            if (std::abs(sys.prtvib) == 1)
            {
                *ivibout << std::fixed << std::setprecision(2) << std::setw(5) << (i + 1) << std::setw(11)
                         << sys.wavenum[i] << std::scientific << std::setprecision(5) << std::setw(14)
                         << sys.freq[i] / 1e9 << std::fixed << std::setprecision(2) << std::setw(12)
                         << sys.freq[i] * h / kb << std::setprecision(8) << std::setw(14) << tmpv0 << std::setw(14)
                         << tmpbot << "\n";
            }
        }

        // Calculate contribution to thermochemistry quantities
        if (std::abs(sys.prtvib) == 1)
        {
            *ivibout << "\n Mode  Wavenumber     ZPE      U(T)-U(0)    U(T)      CV(T)       S(T)\n";
            *ivibout << "         cm^-1      kcal/mol   kcal/mol   kcal/mol  cal/mol/K  cal/mol/K\n";
        }

        for (int i = 0; i < sys.nfreq; ++i)
        {
            double tmpZPE, tmpheat, tmpCV, tmpS;
            getvibcontri(sys, i, tmpZPE, tmpheat, tmpCV, tmpS);
            ZPE += tmpZPE;
            U_vib_heat += tmpheat;
            CV_vib += tmpCV;
            S_vib += tmpS;
            if (std::abs(sys.prtvib) == 1)
            {
                *ivibout << std::fixed << std::setprecision(2) << std::setw(5) << (i + 1) << std::setw(11)
                         << sys.wavenum[i] << std::setprecision(5) << std::setw(11) << tmpZPE / cal2J << std::setw(11)
                         << tmpheat / cal2J << std::setw(11) << (tmpheat + tmpZPE) / cal2J << std::setw(11)
                         << tmpCV / cal2J << std::setw(11) << tmpS / cal2J << "\n";
            }
        }
        double U_vib = U_vib_heat + ZPE;

        if (sys.prtvib == -1)
        {
            vibfile.close();
            std::cout << "Contributions to thermochemistry quantities from every frequency mode have been exported to "
                      << vibcon_filename << " in current folder\n\n";
        }

        std::cout << std::scientific << std::setprecision(6) << " Vibrational q(V=0): " << std::setw(16) << qvib_v0
                  << "\n";
        std::cout << " Vibrational q(bot): " << std::setw(16) << qvib_bot << "\n";
        if (sys.lowVibTreatment != LowVibTreatment::Minenkov)
        {
            std::cout << std::fixed << std::setprecision(3) << " Vibrational U(T)-U(0):" << std::setw(10) << U_vib_heat
                      << " kJ/mol" << std::setw(10) << U_vib_heat / cal2J << " kcal/mol   =H(T)-H(0)\n";
        }
        std::cout << " Vibrational U: " << std::setw(10) << U_vib << " kJ/mol " << std::setw(10) << U_vib / cal2J
                  << " kcal/mol    =H\n";
        std::cout << " Vibrational S: " << std::setw(10) << S_vib << " J/mol/K" << std::setw(10) << S_vib / cal2J
                  << " cal/mol/K   -TS:" << std::setw(8) << -S_vib / cal2J / 1000.0 * sys.T << " kcal/mol\n";
        std::cout << " Vibrational CV:" << std::setw(10) << CV_vib << " J/mol/K" << std::setw(10) << CV_vib / cal2J
                  << " cal/mol/K   =CP\n";
        if (sys.lowVibTreatment != LowVibTreatment::Minenkov)
        {
            std::cout << std::setprecision(2) << " Zero-point energy (ZPE):" << std::setw(10) << ZPE << " kJ/mol,"
                      << std::setw(10) << ZPE / cal2J << " kcal/mol" << std::setprecision(6) << std::setw(12)
                      << ZPE / au2kJ_mol << " a.u.\n";
        }

        // Electron contribution
        std::cout << "\n                 -------- Electronic excitation --------\n"
                  << "                 ---------------------------------------\n";
        // Calculate electronic contributions using elecontri function
        elecontri(sys, q_ele, U_ele, CV_ele, S_ele);
        std::cout << std::scientific << std::setprecision(6) << " Electronic q: " << std::setw(16) << q_ele << "\n";
        std::cout << std::fixed << std::setprecision(3) << " Electronic U: " << std::setw(10) << U_ele << " kJ/mol "
                  << std::setw(10) << U_ele / cal2J << " kcal/mol    =H\n";
        std::cout << " Electronic S: " << std::setw(10) << S_ele << " J/mol/K" << std::setw(10) << S_ele / cal2J
                  << " cal/mol/K   -TS:" << std::setw(8) << -S_ele / cal2J / 1000.0 * sys.T << " kcal/mol\n";
        std::cout << " Electronic CV:" << std::setw(10) << CV_ele << " J/mol/K" << std::setw(10) << CV_ele / cal2J
                  << " cal/mol/K   =CP\n";

        // Total result
        std::cout << "\n\n"
                  << "                       ----------------------------\n"
                  << "                       -------- Final data --------\n"
                  << "                       ----------------------------\n";
        std::cout << std::scientific << std::setprecision(6) << " Total q(V=0):    " << std::setw(16)
                  << q_trans * q_rot * qvib_v0 * q_ele << "\n";
        std::cout << " Total q(bot):    " << std::setw(16) << q_trans * q_rot * qvib_bot * q_ele << "\n";
        std::cout << " Total q(V=0)/NA: " << std::setw(16) << q_trans * q_rot * qvib_v0 * q_ele / NA << "\n";
        std::cout << " Total q(bot)/NA: " << std::setw(16) << q_trans * q_rot * qvib_bot * q_ele / NA << "\n";

        CV_tot = CV_trans + CV_rot + CV_vib + CV_ele;
        CP_tot = CP_trans + CV_rot + CV_vib + CV_ele;
        S_tot  = S_trans + S_rot + S_vib + S_ele;
        std::cout << std::fixed << std::setprecision(3) << " Total CV:" << std::setw(12) << CV_tot << " J/mol/K"
                  << std::setw(12) << CV_tot / cal2J << " cal/mol/K\n";
        std::cout << " Total CP:" << std::setw(12) << CP_tot << " J/mol/K" << std::setw(12) << CP_tot / cal2J
                  << " cal/mol/K\n";
        std::cout << " Total S: " << std::setw(12) << S_tot << " J/mol/K" << std::setw(12) << S_tot / cal2J
                  << " cal/mol/K    -TS:" << std::setw(10) << -S_tot / cal2J / 1000.0 * sys.T << " kcal/mol\n";

        if (sys.T == 0.0)
        {
            thermU = ZPE;
            thermH = ZPE;
            thermG = ZPE;
        }
        else
        {
            thermU = U_trans + U_rot + U_vib + U_ele;
            thermH = H_trans + U_rot + U_vib + U_ele;
            thermG = thermH - sys.T * S_tot / 1000.0;
        }
        sys.thermG = thermG;

        if (sys.lowVibTreatment != LowVibTreatment::Minenkov)
        {
            std::cout << std::setprecision(3) << " Zero point energy (ZPE):" << std::setw(11) << ZPE << " kJ/mol"
                      << std::setw(11) << ZPE / cal2J << " kcal/mol" << std::setprecision(6) << std::setw(11)
                      << ZPE / au2kJ_mol << " a.u.\n";
        }
        std::cout << std::setprecision(3) << " Thermal correction to U:" << std::setw(11) << thermU << " kJ/mol"
                  << std::setw(11) << thermU / cal2J << " kcal/mol" << std::setprecision(6) << std::setw(11)
                  << thermU / au2kJ_mol << " a.u.\n";
        std::cout << " Thermal correction to H:" << std::setw(11) << thermH << " kJ/mol" << std::setw(11)
                  << thermH / cal2J << " kcal/mol" << std::setprecision(6) << std::setw(11) << thermH / au2kJ_mol
                  << " a.u.\n";
        std::cout << " Thermal correction to G:" << std::setw(11) << thermG << " kJ/mol" << std::setw(11)
                  << thermG / cal2J << " kcal/mol" << std::setprecision(6) << std::setw(11) << thermG / au2kJ_mol
                  << " a.u.\n";

        double U0      = sys.E + ZPE / au2kJ_mol;
        double U_final = sys.E + thermU / au2kJ_mol;
        double H_final = sys.E + thermH / au2kJ_mol;
        double G_final = sys.E + thermG / au2kJ_mol;

        std::cout << std::fixed << std::setprecision(7) << " Electronic energy:" << std::setw(19) << sys.E << " a.u.\n";
        if (sys.lowVibTreatment != LowVibTreatment::Minenkov)
        {
            std::cout << " Sum of electronic energy and ZPE, namely U/H/G at 0 K:" << std::setw(19) << U0 << " a.u.\n";
        }
        std::cout << " Sum of electronic energy and thermal correction to U: " << std::setw(19) << U_final << " a.u.\n";
        std::cout << " Sum of electronic energy and thermal correction to H: " << std::setw(19) << H_final << " a.u.\n";
        std::cout << " Sum of electronic energy and thermal correction to G: " << std::setw(19) << G_final << " a.u.\n";

        if (sys.concstr != "0")
        {
            double concnow = 0.0, concspec = std::stod(sys.concstr), Gconc;
            getGconc(sys, concnow, concspec, Gconc);
            std::cout << "\n"
                      << " Present concentration (estimated by ideal gas model):" << std::fixed << std::setprecision(6)
                      << std::setw(10) << concnow << " mol/L\n";
            std::cout << " Concentration specified by \"conc\" parameter:" << std::setw(12) << concspec << " mol/L\n";
            std::cout << " delta-G of conc. change:" << std::fixed << std::setprecision(3) << std::setw(11) << Gconc
                      << " kJ/mol" << std::setw(11) << Gconc / cal2J << " kcal/mol" << std::setprecision(6)
                      << std::setw(11) << Gconc / au2kJ_mol << " a.u.\n";
            std::cout << " Gibbs free energy at specified concentration: " << std::fixed << std::setprecision(7)
                      << std::setw(19) << (G_final + Gconc / au2kJ_mol) << " a.u.\n";
        }
    }


    /**
     * @brief Calculate Gibbs energy correction due to concentration changes
     *
     * Computes the
     * correction to Gibbs free energy when changing from current
     * concentration to a specified concentration
     * using the formula:
     * ΔG_conc = G + RT * ln(conc_current / conc_specified)
     *
     * @param sys
     * SystemData structure with system parameters
     * @param concnow Current concentration in mol/L
     * @param
     * concspec Specified concentration in mol/L
     * @param Gconc [out] Concentration correction to Gibbs energy in
     * kJ/mol
     *
     * @note Used for solution-phase calculations where concentration affects free energy
     *
     * @note If concentrations are invalid (<=0), returns the base Gibbs energy without correction
     */
    void getGconc(const SystemData& sys, const double& concnow, const double& concspec, double& Gconc)
    {
        // Compute concentration correction to Gibbs free energy
        // delta-G_conc = RT * ln(concspec / concnow)
        double calculated_concnow      = sys.P * atm2Pa / (R * sys.T) / 1000.0;  // mol/L
        *const_cast<double*>(&concnow) = calculated_concnow;
        if (calculated_concnow > 0.0 && concspec > 0.0)
        {
            Gconc = R * sys.T * std::log(concspec / calculated_concnow) / 1000.0;  // kJ/mol
        }
        else
        {
            Gconc = 0.0;  // No concentration adjustment if invalid
        }
    }


    /**
     * @brief Calculate electronic contributions to thermodynamic properties
     *
     * Computes the
     * partition function and thermodynamic contributions from
     * electronic energy levels using statistical
     * mechanics. Handles both
     * single electronic states and multiple electronic levels with degeneracies.
     *

     * * @param sys SystemData structure containing electronic level information
     * @param tmpq [out] Electronic
     * partition function
     * @param tmpheat [out] Electronic contribution to internal energy (kJ/mol)
     * @param
     * tmpCV [out] Electronic contribution to heat capacity at constant volume (J/mol·K)
     * @param tmpS [out]
     * Electronic contribution to entropy (J/mol·K)
     *
     * @note Uses Boltzmann statistics for electronic level
     * populations
     * @note Accounts for electronic degeneracy and excitation energies
     * @note For spin
     * multiplicity only, treats as single degenerate ground state
     */
    void elecontri(const SystemData& sys, double& tmpq, double& tmpheat, double& tmpCV, double& tmpS)
    {
        tmpq      = 0.0;
        double t1 = 0.0, t2 = 0.0;

        if (sys.nelevel <= 0 || sys.elevel.size() != static_cast<size_t>(sys.nelevel) ||
            sys.edegen.size() != static_cast<size_t>(sys.nelevel))
        {
            throw std::runtime_error("elecontri: Invalid electron level data");
        }

        for (int ie = 0; ie < sys.nelevel; ++ie)
        {
            double exc = sys.elevel[ie] / au2eV * au2J;
            double ekt = (sys.T > 0.0) ? exc / (kb * sys.T) : 0.0;
            double qi  = sys.edegen[ie] * std::exp(-ekt);
            tmpq += qi;
            t1 += ekt * qi;
            t2 += ekt * ekt * qi;
        }

        if (tmpq <= 0.0)
        {
            throw std::runtime_error("elecontri: Partition function (tmpq) is zero or negative");
        }

        tmpS    = R * std::log(tmpq) + R * t1 / tmpq;
        tmpheat = (sys.T == 0.0) ? 0.0 : R * sys.T * t1 / tmpq / 1000.0;
        tmpCV   = R * t2 / tmpq - R * std::pow(t1 / tmpq, 2);
    }


    // Get contribution of vibration mode i to ZPE, U(T)-U(0), CV, S. In kJ/mol or kJ/mol/K
    //! Imaginary frequecies are ignored


    void getvibcontri(const SystemData& sys, int i, double& tmpZPE, double& tmpheat, double& tmpCV, double& tmpS)
    {
        tmpZPE  = 0.0;
        tmpheat = 0.0;
        tmpCV   = 0.0;
        tmpS    = 0.0;

        if (i < 0 || i >= sys.nfreq || sys.freq[i] <= 0.0)
        {
            return;  // Ignore imaginary or invalid frequencies
        }

        double prefac_trunc = 0.0, term_trunc = 0.0;
        if (sys.lowVibTreatment == LowVibTreatment::Truhlar)
        {
            double freqtrunc = 0.0;
            // Truhlar's QRRHO
            freqtrunc    = sys.ravib * wave2freq;
            prefac_trunc = h * freqtrunc / (kb * sys.T);
            term_trunc   = std::exp(-h * freqtrunc / (kb * sys.T));
        }

        // ZPE
        tmpZPE = sys.wavenum[i] * sys.sclZPE / 2.0 / au2cm_1 * au2kJ_mol;

        // Heating contribution to U, namely U(T)-U(0)
        if (sys.T > 0.0)
        {
            double prefac = h * sys.freq[i] * sys.sclheat / (kb * sys.T);
            double term   = std::exp(-h * sys.freq[i] * sys.sclheat / (kb * sys.T));
            if (sys.lowVibTreatment == LowVibTreatment::Truhlar && sys.wavenum[i] < sys.ravib)
            {
                prefac = prefac_trunc;
                term   = term_trunc;
            }
            if (sys.lowVibTreatment == LowVibTreatment::Minenkov)
            {  // Minenkov: Interpolation between RRHO and free rotor
                double UvRRHO =
                    tmpZPE + R * sys.T * prefac * term / (1.0 - term) / 1000.0;          // Harmonic-oscillator result
                tmpZPE        = 0.0;                                                     // Free rotor has no ZPE
                double Ufree  = R * sys.T / 2.0 / 1000.0;                                // Free-rotor result
                double tmpval = 1.0 + std::pow(sys.intpvib / sys.wavenum[i], 4);         // Denominator part of Eq. 6
                tmpheat       = (1.0 / tmpval) * UvRRHO + (1.0 - 1.0 / tmpval) * Ufree;  // Interpolation
            }
            else
            {  // RRHO
                tmpheat = R * sys.T * prefac * term / (1.0 - term) / 1000.0;
            }
        }

        // CV
        double prefac = h * sys.freq[i] * sys.sclCV / (kb * sys.T);
        double term   = std::exp(-h * sys.freq[i] * sys.sclCV / (kb * sys.T));
        if (sys.lowVibTreatment == LowVibTreatment::Truhlar && sys.wavenum[i] < sys.ravib)
        {
            prefac = prefac_trunc;
            term   = term_trunc;
        }
        tmpCV = R * std::pow(prefac, 2) * term / std::pow(1.0 - term, 2);

        // S
        prefac = h * sys.freq[i] * sys.sclS / (kb * sys.T);
        term   = std::exp(-h * sys.freq[i] * sys.sclS / (kb * sys.T));
        if (sys.lowVibTreatment == LowVibTreatment::Truhlar && sys.wavenum[i] < sys.ravib)
        {
            prefac = prefac_trunc;
            term   = term_trunc;
        }
        tmpS = R * (prefac * term / (1.0 - term) - std::log(1.0 - term));  // RRHO

        if (sys.lowVibTreatment == LowVibTreatment::Grimme || sys.lowVibTreatment == LowVibTreatment::Minenkov)
        {  // Grimme's entropy interpolation
            double miu   = h / (8.0 * M_PI * M_PI * sys.freq[i]);
            double Bav   = 1e-44;  // kg*m^2
            double miup  = miu * Bav / (miu + Bav);
            double Sfree = R * (0.5 + std::log(std::sqrt(8.0 * M_PI * M_PI * M_PI * miup * kb * sys.T / (h * h))));
            double wei   = 1.0 / (1.0 + std::pow(sys.intpvib / sys.wavenum[i], 4));
            tmpS         = wei * tmpS + (1.0 - wei) * Sfree;
        }
    }


    /**
    * @brief Calculate thermodynamic properties using current system temperature and pressure
    *
    *
    * Convenience overload that calls the main calcthermo function using the
    * temperature and pressure values
    * stored in the SystemData structure.
    *
    * @param sys SystemData structure with molecular data and T/P
    * parameters
    * @param corrU [out] Thermal correction to internal energy (kJ/mol)
    * @param corrH [out]
    * Thermal correction to enthalpy (kJ/mol)
    * @param corrG [out] Thermal correction to Gibbs energy (kJ/mol)

    * * @param S [out] Total entropy (J/mol·K)
    * @param CV [out] Constant volume heat capacity (J/mol·K)
    *
    * @param CP [out] Constant pressure heat capacity (J/mol·K)
    * @param QV [out] Vibrational partition function

    * * @param Qbot [out] Bottom partition function (rotational + electronic)
    *
    * @note Delegates to the full
    * calcthermo function with sys.T and sys.P
    */
    void calcthermo(SystemData& sys,
                    double&     corrU,
                    double&     corrH,
                    double&     corrG,
                    double&     S,
                    double&     CV,
                    double&     CP,
                    double&     QV,
                    double&     Qbot)
    {
        // Call the overloaded version with current T and P from sys
        calcthermo(sys, sys.T, sys.P, corrU, corrH, corrG, S, CV, CP, QV, Qbot);
    }

    /**
     * @brief Calculate total vibrational contributions to thermodynamic properties
     *
     * Computes the
     * total thermodynamic contributions from all vibrational modes
     * at a given temperature, including zero-point
     * energy, thermal energy,
     * heat capacity, entropy, and partition function.
     *
     * @param sys
     * SystemData structure containing vibrational frequency data
     * @param T Temperature in Kelvin
     * @param
     * U_vib [out] Total vibrational internal energy (including ZPE)
     * @param CV_vib [out] Total vibrational heat
     * capacity at constant volume
     * @param S_vib [out] Total vibrational entropy
     * @param QV [out] Total
     * vibrational partition function
     *
     * @note Applies scaling factors for ZPE, thermal energy, entropy, and
     * heat capacity
     * @note Handles low-frequency modes according to user settings
     * @note Ignores imaginary
     * frequencies (negative values)
     */
    void getvibcontri(const SystemData& sys, double T, double& U_vib, double& CV_vib, double& S_vib, double& QV)
    {
        U_vib  = 0.0;
        CV_vib = 0.0;
        S_vib  = 0.0;
        QV     = 1.0;

        if (sys.nfreq == 0)
            return;

        double ZPE = 0.0;
        for (int i = 0; i < sys.nfreq; ++i)
        {
            if (sys.freq[i] > 0)
            {
                double freq_scaled = sys.freq[i] * sys.sclZPE;
                ZPE += 0.5 * h * freq_scaled;
            }
        }

        // Calculate vibrational contributions using scaled frequencies
        double U_vib_heat = 0.0;
        for (int i = 0; i < sys.nfreq; ++i)
        {
            if (sys.freq[i] > 0)
            {
                double freq_heat = sys.freq[i] * sys.sclheat;
                double x         = h * freq_heat / (kb * T);
                if (x < 700)
                {  // Avoid overflow
                    double exp_x = std::exp(x);
                    U_vib_heat += h * freq_heat / (exp_x - 1.0);
                    CV_vib += h * freq_heat * x * exp_x / (kb * T * (exp_x - 1.0) * (exp_x - 1.0));
                }
            }
        }

        // Calculate entropy
        for (int i = 0; i < sys.nfreq; ++i)
        {
            if (sys.freq[i] > 0)
            {
                double freq_S = sys.freq[i] * sys.sclS;
                double x      = h * freq_S / (kb * T);
                if (x < 700)
                {  // Avoid overflow
                    double exp_x = std::exp(x);
                    S_vib += kb * (x / (exp_x - 1.0) - std::log(1.0 - 1.0 / exp_x));
                }
            }
        }

        // Calculate partition function
        for (int i = 0; i < sys.nfreq; ++i)
        {
            if (sys.freq[i] > 0)
            {
                double freq_CV = sys.freq[i] * sys.sclCV;
                double x       = h * freq_CV / (kb * T);
                if (x < 700)
                {  // Avoid overflow
                    QV *= 1.0 / (1.0 - std::exp(-x));
                }
            }
        }

        U_vib = U_vib_heat + ZPE;
    }

    /**
     * @brief Calculate moments of inertia for the molecular system
     *
     * Computes the moment of inertia tensor and principal moments of inertia
     * based on the molecular geometry and atomic masses. The calculation
     * involves finding the center of mass and constructing the inertia tensor.
     *
     * @param sys SystemData structure containing molecular geometry and masses
     * @note Updates sys.inertmat (3x3 inertia tensor) and sys.inert (principal moments)
     * @note Units: amu·Bohr² for inertia tensor, converted to standard units internally
     */
    void calcinertia(SystemData& sys)
    {
        // Calculate center of mass coordinates
        double cenmassx = 0.0, cenmassy = 0.0, cenmassz = 0.0;
        for (const auto& atom : sys.a)
        {
            cenmassx += atom.x * atom.mass;
            cenmassy += atom.y * atom.mass;
            cenmassz += atom.z * atom.mass;
        }
        cenmassx /= sys.totmass;
        cenmassy /= sys.totmass;
        cenmassz /= sys.totmass;

        // Build inertia matrix
        double sumxx = 0.0, sumyy = 0.0, sumzz = 0.0;
        double sumxy = 0.0, sumxz = 0.0, sumyz = 0.0;
        for (const auto& atom : sys.a)
        {
            double dx = atom.x - cenmassx;
            double dy = atom.y - cenmassy;
            double dz = atom.z - cenmassz;
            sumxx += atom.mass * (dy * dy + dz * dz);
            sumyy += atom.mass * (dx * dx + dz * dz);
            sumzz += atom.mass * (dx * dx + dy * dy);
            sumxy -= atom.mass * dx * dy;
            sumxz -= atom.mass * dx * dz;
            sumyz -= atom.mass * dy * dz;
        }

        sys.inertmat[0][0] = sumxx;
        sys.inertmat[1][1] = sumyy;
        sys.inertmat[2][2] = sumzz;
        sys.inertmat[0][1] = sys.inertmat[1][0] = sumxy;
        sys.inertmat[0][2] = sys.inertmat[2][0] = sumxz;
        sys.inertmat[1][2] = sys.inertmat[2][1] = sumyz;

        // Convert from amu*Ang^2 to amu*Bohr^2
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                sys.inertmat[i][j] /= (b2a * b2a);
            }
        }

        // Diagonalize
        std::vector<std::vector<double>> mat(3, std::vector<double>(3));
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                mat[i][j] = sys.inertmat[i][j];
            }
        }
        std::vector<std::vector<double>> eigvec(3, std::vector<double>(3));
        util::diagmat(mat, eigvec, sys.inert, 300, 1e-12);

        // Update inertmat with diagonalized matrix
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                sys.inertmat[i][j] = mat[i][j];
            }
        }
    }

}  // namespace calc
