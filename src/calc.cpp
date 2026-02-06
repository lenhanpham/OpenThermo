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
#include "omp_config.h"
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

    // Forward declarations
    void elecontri(const SystemData& sys, double& tmpq, double& tmpheat, double& tmpCV, double& tmpS);
    void elecontri(const SystemData& sys, double T, double& tmpq, double& tmpheat, double& tmpCV, double& tmpS);
    void getvibcontri(const SystemData& sys, int i, double T, double& tmpZPE, double& tmpheat, double& tmpCV, double& tmpS);


    /**
     * @brief Check if a file exists and is accessible
     *
     * @param filename Path to the file to check
     * @return true if file exists and can be opened, false otherwise
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
     * Computes all thermochemical properties (internal energy, enthalpy, Gibbs energy,
     * entropy, heat capacities, partition functions) using statistical mechanics.
     * Includes translational, rotational, vibrational, and electronic contributions.
     *
     * @param sys SystemData structure with molecular data (const, not modified)
     * @param T Temperature in Kelvin
     * @param P Pressure in atmospheres
     * @param corrU [out] Thermal correction to internal energy (kJ/mol)
     * @param corrH [out] Thermal correction to enthalpy (kJ/mol)
     * @param corrG [out] Thermal correction to Gibbs energy (kJ/mol)
     * @param S [out] Total entropy (J/mol/K)
     * @param CV [out] Heat capacity at constant volume (J/mol/K)
     * @param CP [out] Heat capacity at constant pressure (J/mol/K)
     * @param QV [out] Vibrational partition function (q(V=0)/NA)
     * @param Qbot [out] Bottom-of-well partition function (q(bot)/NA)
     *
     * @note Uses fused vibrational loop with precomputed invariants for performance
     * @note Supports Harmonic, Truhlar, Grimme, and Minenkov low-frequency treatments
     */
    void calcthermo(const SystemData& sys,
                    double      T,
                    double      P,
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

        const double b2m_sq = (b2a * 1e-10) * (b2a * 1e-10);

        if (sys.ipmode == 0)
        {
            double P_Pa = P * atm2Pa;
            double trans_base = 2.0 * M_PI * (sys.totmass * amu2kg) * kb * T / (h * h);
            q_trans = trans_base * std::sqrt(trans_base) * R * T / P_Pa;
            CV_trans = 1.5 * R;
            CP_trans = 2.5 * R;
            U_trans  = 1.5 * R * T / 1000.0;
            H_trans  = 2.5 * R * T / 1000.0;
            S_trans  = R * (std::log(q_trans / NA) + 2.5);
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

        if (sys.ipmode == 0)
        {
            double sum_inert = sys.inert[0] + sys.inert[1] + sys.inert[2];
            if (sum_inert < 1e-10)
            {
                q_rot  = 1.0;
                U_rot  = 0.0;
                CV_rot = 0.0;
                S_rot  = 0.0;
            }
            else
            {
                std::array<double, 3> inertkg;
                for (int i = 0; i < 3; ++i)
                {
                    inertkg[i] = sys.inert[i] * amu2kg * b2m_sq;
                }
                if (sys.ilinear == 1)
                {  // Linear molecule
                    q_rot  = 8.0 * M_PI * M_PI * inertkg[2] * kb * T / sys.rotsym / (h * h);
                    U_rot  = R * T / 1000.0;
                    CV_rot = R;
                    S_rot  = R * (std::log(q_rot) + 1.0);
                }
                else
                {  // Non-linear molecule
                    double rot_base = 2.0 * M_PI * kb * T;
                    q_rot = 8.0 * M_PI * M_PI / sys.rotsym / (h * h * h) *
                            rot_base * std::sqrt(rot_base) *
                            std::sqrt(inertkg[0] * inertkg[1] * inertkg[2]);
                    U_rot  = 1.5 * R * T / 1000.0;
                    CV_rot = 1.5 * R;
                    S_rot  = R * (std::log(q_rot) + 1.5);
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

        const double h_over_kbT = h / (kb * T);
        const double RT_1000 = R * T / 1000.0;
        const double zpe_factor = sys.sclZPE / 2.0 / au2cm_1 * au2kJ_mol;
        const int nfreq = sys.nfreq;
        const auto lowVib = sys.lowVibTreatment;
        const double ravib_freq = sys.ravib * wave2freq;
        const double sclheat = sys.sclheat;
        const double sclCV = sys.sclCV;
        const double sclS = sys.sclS;
        const bool uniform_scaling = (sclheat == 1.0 && sclCV == 1.0 && sclS == 1.0);
        const double prefac_trunc = (lowVib == LowVibTreatment::Truhlar) ? h_over_kbT * ravib_freq : 0.0;
        const double term_trunc = (lowVib == LowVibTreatment::Truhlar) ? std::exp(-prefac_trunc) : 0.0;
        const bool do_grimme_interp = (lowVib == LowVibTreatment::Grimme || lowVib == LowVibTreatment::Minenkov);
        constexpr double eight_pi2 = 8.0 * M_PI * M_PI;
        const double grimme_log_base = 8.0 * M_PI * M_PI * M_PI * 1e-44 * kb * T / (h * h);

#ifdef _OPENMP
#pragma omp parallel for reduction(*:qvib_v0,qvib_bot) reduction(+:ZPE,U_vib_heat,CV_vib,S_vib) if(nfreq > 50)
#endif
        for (int i = 0; i < nfreq; ++i)
        {
            double fi = sys.freq[i];
            if (fi <= 0.0)
                continue;
            double wi = sys.wavenum[i];
            bool truhlar_active = (lowVib == LowVibTreatment::Truhlar && wi < sys.ravib);

            double freqtmp = truhlar_active ? ravib_freq : fi;
            double x = h_over_kbT * freqtmp;
            double exp_neg_x = std::exp(-x);
            double one_minus_exp = 1.0 - exp_neg_x;
            qvib_v0 *= 1.0 / one_minus_exp;
            qvib_bot *= std::sqrt(exp_neg_x) / one_minus_exp;

            double local_ZPE = wi * zpe_factor;

            double pf_base = 0.0, tm_base = 0.0;
            if (uniform_scaling)
            {
                pf_base = h_over_kbT * fi;
                tm_base = (truhlar_active) ? term_trunc : std::exp(-pf_base);
                if (truhlar_active) pf_base = prefac_trunc;
            }

            double local_heat = 0.0;
            if (T > 0.0)
            {
                double pf_h, tm_h;
                if (uniform_scaling)
                {
                    pf_h = pf_base; tm_h = tm_base;
                }
                else
                {
                    pf_h = h_over_kbT * fi * sclheat;
                    tm_h = std::exp(-pf_h);
                    if (truhlar_active) { pf_h = prefac_trunc; tm_h = term_trunc; }
                }

                if (lowVib == LowVibTreatment::Minenkov)
                {
                    double UvRRHO = local_ZPE + RT_1000 * pf_h * tm_h / (1.0 - tm_h);
                    local_ZPE = 0.0;
                    double Ufree = RT_1000 * 0.5;
                    double r = sys.intpvib / wi;
                    double r2 = r * r;
                    double tmpval = 1.0 + r2 * r2;
                    local_heat = (1.0 / tmpval) * UvRRHO + (1.0 - 1.0 / tmpval) * Ufree;
                }
                else
                {
                    local_heat = RT_1000 * pf_h * tm_h / (1.0 - tm_h);
                }
            }

            double pf_cv, tm_cv;
            if (uniform_scaling)
            {
                pf_cv = pf_base; tm_cv = tm_base;
            }
            else
            {
                pf_cv = h_over_kbT * fi * sclCV;
                tm_cv = std::exp(-pf_cv);
                if (truhlar_active) { pf_cv = prefac_trunc; tm_cv = term_trunc; }
            }
            double omt_cv = 1.0 - tm_cv;
            double local_CV = R * pf_cv * pf_cv * tm_cv / (omt_cv * omt_cv);

            double pf_s, tm_s;
            if (uniform_scaling)
            {
                pf_s = pf_base; tm_s = tm_base;
            }
            else
            {
                pf_s = h_over_kbT * fi * sclS;
                tm_s = std::exp(-pf_s);
                if (truhlar_active) { pf_s = prefac_trunc; tm_s = term_trunc; }
            }
            double local_S = R * (pf_s * tm_s / (1.0 - tm_s) - std::log(1.0 - tm_s));

            if (do_grimme_interp)
            {
                double miu = h / (eight_pi2 * fi);
                constexpr double Bav = 1e-44;
                double miup = miu * Bav / (miu + Bav);
                double Sfree = R * (0.5 + 0.5 * std::log(grimme_log_base * miup / Bav));
                double gr = sys.intpvib / wi;
                double gr2 = gr * gr;
                double wei = 1.0 / (1.0 + gr2 * gr2);
                local_S = wei * local_S + (1.0 - wei) * Sfree;
            }

            ZPE += local_ZPE;
            U_vib_heat += local_heat;
            CV_vib += local_CV;
            S_vib += local_S;
        }
        double U_vib = U_vib_heat + ZPE;

        elecontri(sys, T, q_ele, U_ele, CV_ele, S_ele);

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
     * @brief Display detailed thermodynamic properties for a single (T, P) point
     *
     * Computes and prints all thermochemical contributions (translational, rotational,
     * vibrational, electronic) and final corrected energies to stdout. Also handles
     * concentration corrections and optional per-mode vibrational output.
     *
     * @param sys SystemData structure with molecular data and parameters
     *
     * @note Writes thermG back to sys for downstream use
     * @note When sys.prtvib != 0, outputs individual vibrational mode contributions
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
            double trans_base = 2.0 * M_PI * (sys.totmass * amu2kg) * kb * sys.T / (h * h);
            q_trans = trans_base * std::sqrt(trans_base) * R * sys.T / P_Pa;
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
                const double b2m_sq_show = (b2a * 1e-10) * (b2a * 1e-10);
                std::array<double, 3> inertkg;
                for (int i = 0; i < 3; ++i)
                {
                    inertkg[i] = sys.inert[i] * amu2kg * b2m_sq_show;
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
                    double rot_base = 2.0 * M_PI * kb * sys.T;
                    q_rot = 8.0 * M_PI * M_PI / sys.rotsym / (h * h * h) *
                            rot_base * std::sqrt(rot_base) *
                            std::sqrt(inertkg[0] * inertkg[1] * inertkg[2]);
                    U_rot  = 1.5 * R * sys.T / 1000.0;
                    CV_rot = 1.5 * R;
                    S_rot  = R * (std::log(q_rot) + 1.5);
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
     * @brief Calculate concentration-dependent Gibbs energy correction
     *
     * Computes the correction to Gibbs energy when converting from ideal gas
     * standard state to a specified solution concentration using RT*ln(c_spec/c_gas).
     *
     * @param sys SystemData with temperature and pressure
     * @param concnow [out] Calculated gas-phase concentration (mol/L)
     * @param concspec Specified target concentration (mol/L)
     * @param Gconc [out] Concentration correction to Gibbs energy (kJ/mol)
     */
    void getGconc(const SystemData& sys, const double& concnow, const double& concspec, double& Gconc)
    {
        double calculated_concnow      = sys.P * atm2Pa / (R * sys.T) / 1000.0;  // mol/L
        *const_cast<double*>(&concnow) = calculated_concnow;
        if (calculated_concnow > 0.0 && concspec > 0.0)
        {
            Gconc = R * sys.T * std::log(concspec / calculated_concnow) / 1000.0;  // kJ/mol
        }
        else
        {
            Gconc = 0.0;
        }
    }


    /**
     * @brief Calculate electronic contributions to thermodynamic properties
     *
     * Computes the electronic partition function and its contributions to
     * entropy, thermal energy, and heat capacity from multi-level electronic states.
     *
     * @param sys SystemData with electronic energy levels and degeneracies
     * @param T Temperature in Kelvin
     * @param tmpq [out] Electronic partition function
     * @param tmpheat [out] Electronic thermal energy contribution (kJ/mol)
     * @param tmpCV [out] Electronic heat capacity contribution (J/mol/K)
     * @param tmpS [out] Electronic entropy contribution (J/mol/K)
     */
    void elecontri(const SystemData& sys, double T, double& tmpq, double& tmpheat, double& tmpCV, double& tmpS)
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
            double ekt = (T > 0.0) ? exc / (kb * T) : 0.0;
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
        tmpheat = (T == 0.0) ? 0.0 : R * T * t1 / tmpq / 1000.0;
        double t1_over_q = t1 / tmpq;
        tmpCV   = R * t2 / tmpq - R * t1_over_q * t1_over_q;
    }

    void elecontri(const SystemData& sys, double& tmpq, double& tmpheat, double& tmpCV, double& tmpS)
    {
        elecontri(sys, sys.T, tmpq, tmpheat, tmpCV, tmpS);
    }


    /**
     * @brief Calculate vibrational contributions from a single mode
     *
     * Computes the thermodynamic contributions (ZPE, thermal energy, heat capacity,
     * entropy) for vibrational mode i, applying the configured low-frequency treatment
     * (Harmonic, Truhlar, Grimme, or Minenkov).
     *
     * @param sys SystemData with vibrational frequencies and treatment parameters
     * @param i Index of the vibrational mode (0-based)
     * @param T Temperature in Kelvin
     * @param tmpZPE [out] Zero-point energy contribution (kJ/mol)
     * @param tmpheat [out] Thermal energy contribution (kJ/mol)
     * @param tmpCV [out] Heat capacity contribution (J/mol/K)
     * @param tmpS [out] Entropy contribution (J/mol/K)
     */
    void getvibcontri(const SystemData& sys, int i, double T, double& tmpZPE, double& tmpheat, double& tmpCV, double& tmpS)
    {
        tmpZPE  = 0.0;
        tmpheat = 0.0;
        tmpCV   = 0.0;
        tmpS    = 0.0;

        if (i < 0 || i >= sys.nfreq || sys.freq[i] <= 0.0)
        {
            return;
        }

        double prefac_trunc = 0.0, term_trunc = 0.0;
        if (sys.lowVibTreatment == LowVibTreatment::Truhlar)
        {
            double freqtrunc = sys.ravib * wave2freq;
            prefac_trunc = h * freqtrunc / (kb * T);
            term_trunc   = std::exp(-h * freqtrunc / (kb * T));
        }

        tmpZPE = sys.wavenum[i] * sys.sclZPE / 2.0 / au2cm_1 * au2kJ_mol;

        if (T > 0.0)
        {
            double prefac = h * sys.freq[i] * sys.sclheat / (kb * T);
            double term   = std::exp(-h * sys.freq[i] * sys.sclheat / (kb * T));
            if (sys.lowVibTreatment == LowVibTreatment::Truhlar && sys.wavenum[i] < sys.ravib)
            {
                prefac = prefac_trunc;
                term   = term_trunc;
            }
            if (sys.lowVibTreatment == LowVibTreatment::Minenkov)
            {
                double UvRRHO = tmpZPE + R * T * prefac * term / (1.0 - term) / 1000.0;
                tmpZPE = 0.0;
                double Ufree = R * T / 2.0 / 1000.0;
                double intpvib_ratio = sys.intpvib / sys.wavenum[i];
                double intpvib_r2 = intpvib_ratio * intpvib_ratio;
                double tmpval = 1.0 + intpvib_r2 * intpvib_r2;
                tmpheat = (1.0 / tmpval) * UvRRHO + (1.0 - 1.0 / tmpval) * Ufree;
            }
            else
            {
                tmpheat = R * T * prefac * term / (1.0 - term) / 1000.0;
            }
        }

        double prefac = h * sys.freq[i] * sys.sclCV / (kb * T);
        double term   = std::exp(-h * sys.freq[i] * sys.sclCV / (kb * T));
        if (sys.lowVibTreatment == LowVibTreatment::Truhlar && sys.wavenum[i] < sys.ravib)
        {
            prefac = prefac_trunc;
            term   = term_trunc;
        }
        double one_minus_term = 1.0 - term;
        tmpCV = R * prefac * prefac * term / (one_minus_term * one_minus_term);

        prefac = h * sys.freq[i] * sys.sclS / (kb * T);
        term   = std::exp(-h * sys.freq[i] * sys.sclS / (kb * T));
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
            double Sfree = R * (0.5 + std::log(std::sqrt(8.0 * M_PI * M_PI * M_PI * miup * kb * T / (h * h))));
            double grim_ratio = sys.intpvib / sys.wavenum[i];
            double grim_r2 = grim_ratio * grim_ratio;
            double wei   = 1.0 / (1.0 + grim_r2 * grim_r2);
            tmpS         = wei * tmpS + (1.0 - wei) * Sfree;
        }
    }

    void getvibcontri(const SystemData& sys, int i, double& tmpZPE, double& tmpheat, double& tmpCV, double& tmpS)
    {
        getvibcontri(sys, i, sys.T, tmpZPE, tmpheat, tmpCV, tmpS);
    }


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
        const SystemData& csys = sys;
        calcthermo(csys, sys.T, sys.P, corrU, corrH, corrG, S, CV, CP, QV, Qbot);
    }

    /**
     * @brief Calculate total vibrational contributions from all modes at temperature T
     *
     * @param sys SystemData with vibrational frequency data
     * @param T Temperature in Kelvin
     * @param U_vib [out] Total vibrational internal energy (including ZPE)
     * @param CV_vib [out] Total vibrational heat capacity at constant volume
     * @param S_vib [out] Total vibrational entropy
     * @param QV [out] Total vibrational partition function
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
        const int natoms_inert = static_cast<int>(sys.a.size());
#ifdef _OPENMP
#pragma omp parallel for reduction(+:cenmassx) reduction(+:cenmassy) reduction(+:cenmassz) if(natoms_inert > 50)
#endif
        for (int ia = 0; ia < natoms_inert; ++ia)
        {
            cenmassx += sys.a[ia].x * sys.a[ia].mass;
            cenmassy += sys.a[ia].y * sys.a[ia].mass;
            cenmassz += sys.a[ia].z * sys.a[ia].mass;
        }
        cenmassx /= sys.totmass;
        cenmassy /= sys.totmass;
        cenmassz /= sys.totmass;

        // Build inertia matrix
        double sumxx = 0.0, sumyy = 0.0, sumzz = 0.0;
        double sumxy = 0.0, sumxz = 0.0, sumyz = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sumxx) reduction(+:sumyy) reduction(+:sumzz) reduction(+:sumxy) reduction(+:sumxz) reduction(+:sumyz) if(natoms_inert > 50)
#endif
        for (int ia = 0; ia < natoms_inert; ++ia)
        {
            double dx = sys.a[ia].x - cenmassx;
            double dy = sys.a[ia].y - cenmassy;
            double dz = sys.a[ia].z - cenmassz;
            sumxx += sys.a[ia].mass * (dy * dy + dz * dz);
            sumyy += sys.a[ia].mass * (dx * dx + dz * dz);
            sumzz += sys.a[ia].mass * (dx * dx + dy * dy);
            sumxy -= sys.a[ia].mass * dx * dy;
            sumxz -= sys.a[ia].mass * dx * dz;
            sumyz -= sys.a[ia].mass * dy * dz;
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
