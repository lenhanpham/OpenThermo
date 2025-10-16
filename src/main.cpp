/**
 * @file main.cpp
 * @brief Main entry point for OpenThermo molecular thermochemistry program
 * @author Le Nhan Pham
 * @date 2025
 *
 * This file contains the main function that orchestrates the entire OpenThermo
 * calculation workflow, including input parsing, molecular data Processing,
 * thermochemistry calculations, and result output.
 */

#include "atommass.h"
#include "calc.h"
#include "chemsys.h"
#include "help_utils.h"
#include "loadfile.h"
#include "symmetry.h"
#include "util.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


/**
 * @brief Extract basename without extension from a file path
 * @param filepath Full file path
 * @return Basename
 * without extension
 */
auto get_basename_without_extension(const std::string& filepath) -> std::string
{
    // Find last directory separator
    size_t      last_slash = filepath.find_last_of("/\\");
    std::string filename   = (last_slash != std::string::npos) ? filepath.substr(last_slash + 1) : filepath;

    // Find last dot for extension
    size_t last_dot = filename.rfind('.');
    return (last_dot != std::string::npos) ? filename.substr(0, last_dot) : filename;
}

/**
 * @brief Main entry point for OpenThermo
 *
 * This function initializes the program, processes command-line arguments,
 * loads molecular data from input files, performs thermochemistry calculations,
 * and outputs the results.
 *
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return Exit status (0 for success, non-zero for errors)
 */
auto main(int argc, char* argv[]) -> int
{
    try
    {
        SystemData            sys;                       // Main system data structure
        std::array<double, 3> rotcst = {0.0, 0.0, 0.0};  // Rotational constants

        // Print program information
        std::cout << "  " << "                                                                                     \n"
                  << "  " << "   ***********************************************************************     " << " \n"
                  << "  " << "                                OPENTHERMO                                     " << " \n"
                  << "  " << "   ***********************************************************************     " << " \n"
                  << "# " << "-------------------------------------------------------------------------------" << "#\n"
                  << "# " << "Version 0.001.1  Release date: 2025                                            " << "#\n"
                  << "# " << "Developer: Le Nhan Pham                                                        " << "#\n"
                  << "# " << "https://github.com/lenhanpham/openthermo                                       " << "#\n"
                  << "# " << "-------------------------------------------------------------------------------" << "#\n";

        std::cout << "  " << "                                                                                     \n"
                  << "  " << "                                                                               " << " \n"
                  << "  " << "Please cite this preprint if you use OpenThermo for your research              " << " \n"
                  << "  " << "                                                                               " << " \n"
                  << "# " << "-------------------------------------------------------------------------------" << "#\n"
                  << "# " << "L.N Pham, \"OpenThermo A Comprehensive C++ Program for Calculation of           " << "#\n"
                  << "# " << "Thermochemical Properties\" 2025, http://dx.doi.org/10.13140/RG.2.2.22380.63363 " << "#\n"
                  << "# " << "-------------------------------------------------------------------------------" << "#\n";


        // Handle help options before any other processing
        if (argc > 1)
        {
            std::string first_arg = argv[1];
            if (first_arg == "--help")
            {
                HelpUtils::print_help(argv[0]);
                return 0;
            }
            else if (first_arg == "--help-input")
            {
                HelpUtils::print_input_help();
                return 0;
            }
            else if (first_arg == "--help-output")
            {
                HelpUtils::print_output_help();
                return 0;
            }
            else if (first_arg == "--help-settings")
            {
                HelpUtils::print_settings_help();
                return 0;
            }
            else if (first_arg == "--create-config")
            {
                try
                {
                    util::create_default_settings_file();
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Error creating settings file: " << e.what() << "\n";
                    return 1;
                }
                return 0;
            }
            else if (first_arg.substr(0, 7) == "--help-")
            {
                std::string option = first_arg.substr(7);  // Remove "--help-"
                HelpUtils::print_option_help(option, argv[0]);
                return 0;
            }
        }

        // Initialize isotope mass table
        atommass::initmass(sys);

        // Load running parameters
        int narg = argc - 1;
        for (int iarg = 1; iarg <= narg; ++iarg)
        {
            std::string inputArgs = argv[iarg];
            if (inputArgs == "-noset")
            {
                sys.inoset = 1;
                break;
            }
        }
        if (sys.inoset == 1)
        {
            std::cout << "\"-noset\" is set: Setting parameters from settings.ini are ignored \n";
        }
        else
        {
            util::loadsettings(sys);
        }
        if (narg > 1)
        {
            std::vector<std::string> args(argv, argv + argc);
            util::loadarguments(sys, argc, args);
        }

        // Print running parameters
        std::cout << "\n                   --- Summary of Current Parameters ---\n\nRunning parameters:\n";
        if (sys.prtvib == 1)
        {
            std::cout << "Printing individual contribution of vibration modes: Yes\n";
        }
        else if (sys.prtvib == -1)
        {
            std::cout << "Printing individual contribution of vibration modes: Yes, to <basename>.vibcon file\n";
        }
        else
        {
            std::cout << "Printing individual contribution of vibration modes: No\n";
        }
        if (sys.Tstep == 0.0)
        {
            std::cout << " Temperature:     " << std::fixed << std::setprecision(3) << std::setw(12) << sys.T << " K\n";
        }
        else
        {
            std::cout << " Temperature scan, from " << std::fixed << std::setprecision(3) << std::setw(10) << sys.Tlow
                      << " to " << std::setw(10) << sys.Thigh << ", step: " << std::setw(8) << sys.Tstep << " K\n";
        }
        if (sys.Pstep == 0.0)
        {
            std::cout << " Pressure:      " << std::fixed << std::setprecision(3) << std::setw(12) << sys.P << " atm\n";
        }
        else
        {
            std::cout << " Pressure scan, from " << std::fixed << std::setprecision(3) << std::setw(10) << sys.Plow
                      << " to " << std::setw(10) << sys.Phigh << ", step: " << std::setw(8) << sys.Pstep << " atm\n";
        }
        if (sys.concstr != "0")
        {
            std::cout << " Concentration: " << std::fixed << std::setprecision(3) << std::setw(12)
                      << std::stod(sys.concstr) << " mol/L\n";
        }

        std::cout << " Scaling factor of vibrational frequencies for ZPE:       " << std::fixed << std::setprecision(4)
                  << std::setw(8) << sys.sclZPE << "\n"
                  << " Scaling factor of vibrational frequencies for U(T)-U(0): " << std::setw(8) << sys.sclheat << "\n"
                  << " Scaling factor of vibrational frequencies for S(T):      " << std::setw(8) << sys.sclS << "\n"
                  << " Scaling factor of vibrational frequencies for CV:        " << std::setw(8) << sys.sclCV << "\n";
        if (sys.lowVibTreatment == LowVibTreatment::Harmonic)
        {
            std::cout << "Low frequencies treatment: Harmonic approximation\n";
        }
        else if (sys.lowVibTreatment == LowVibTreatment::Truhlar)
        {
            std::cout << " Low frequencies treatment: Raising low frequencies (Truhlar's treatment)\n"
                      << " Lower frequencies will be raised to " << std::fixed << std::setprecision(2) << sys.ravib
                      << " cm^-1 during calculating S, U(T)-U(0), CV and q\n";
        }
        else if (sys.lowVibTreatment == LowVibTreatment::Grimme)
        {
            std::cout << " Low frequencies treatment: Grimme's interpolation for entropy\n";
        }
        else if (sys.lowVibTreatment == LowVibTreatment::Minenkov)
        {
            std::cout << " Low frequencies treatment: Minenkov's interpolation for entropy and internal energy\n";
        }
        if (sys.lowVibTreatment == LowVibTreatment::Grimme || sys.lowVibTreatment == LowVibTreatment::Minenkov)
        {
            std::cout << " Vibrational frequency threshold used in the interpolation is " << std::fixed
                      << std::setprecision(2) << sys.intpvib << " cm^-1\n";
        }
        if (sys.imagreal != 0.0)
        {
            std::cout << " Imaginary frequencies with norm < " << std::fixed << std::setprecision(2) << sys.imagreal
                      << " cm^-1 will be treated as real frequencies\n";
        }

        // Load input file
        if (narg >= 1)
        {
            sys.inputfile = argv[1];
        }
        if (sys.inputfile.empty())
        {
            std::cout << "\nInput file path, e.g. D:\\your_dir\\your_calc.log\n"
                      << " OpenThermo supports Gaussian, ORCA, GAMESS-US, NWChem, CP2K, and VASP "
                         "\n";
            while (true)
            {
                std::getline(std::cin, sys.inputfile);
                sys.inputfile.erase(std::remove(sys.inputfile.begin(), sys.inputfile.end(), '"'), sys.inputfile.end());
                std::ifstream check(sys.inputfile);
                sys.alive = check.good();
                check.close();
                if (sys.alive)
                    break;
                std::cout << "Cannot find the file, input again!\n";
            }
        }
        else
        {
            std::ifstream check(sys.inputfile);
            sys.alive = check.good();
            check.close();
            if (!sys.alive)
            {
                std::cout << "Error: Unable to find " << sys.inputfile << "\n";
                std::exit(1);
            }
        }

        // Print start message
        auto        start_now      = std::chrono::system_clock::now();
        std::time_t start_now_time = std::chrono::system_clock::to_time_t(start_now);
        // std::string basename       = get_basename_without_extension(sys.inputfile);
        std::cout << "                      -------- End of Summary --------\n";
        std::cout << "\n";
        std::cout << "OpenThermo started to process " << sys.inputfile << " at "
                  << std::ctime(&start_now_time);  // << "\n";

        // Process input file
        if (sys.inputfile.find(".otm") != std::string::npos)
        {
            std::cout << "\n Processing data from " << sys.inputfile << "\n";
            LoadFile::loadotm(sys);
        }
        else
        {
            auto qcprog = util::deterprog(sys);
            if (qcprog != util::QuantumChemistryProgram::Unknown)
            {
                std::cout << "\n";
                if (sys.massmod == 1)
                    std::cout << " Atomic masses used: Element\n";
                if (sys.massmod == 2)
                    std::cout << " Atomic masses used: Most abundant isotope\n";
                if (sys.massmod == 3)
                    std::cout << " Atomic masses used: Read from quantum chemical output\n";
                if (qcprog == util::QuantumChemistryProgram::Gaussian)
                {
                    std::cout << "Processing Gaussian output file...\n";
                    LoadFile::loadgau(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::Orca)
                {
                    std::cout << "Processing ORCA output file...\n";
                    LoadFile::loadorca(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::Gamess)
                {
                    std::cout << "Processing GAMESS-US output file...\n";
                    LoadFile::loadgms(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::Nwchem)
                {
                    std::cout << "Processing NWChem output file...\n";
                    LoadFile::loadnw(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::Cp2k)
                {
                    std::cout << "Processing CP2K output file...\n";
                    LoadFile::loadCP2K(sys);
                    if (sys.ipmode == 0)
                    {
                        std::cout << " Note: If your system is not isolated (periodic crystals, slabs or adsorbate on "
                                     "surface), \n"
                                     "you may want to set"
                                     "\"ipmode\" = 1 settings.ini in order to ignore translation and rotation "
                                     "contributions. \n"
                                  << "This is typical for condensed materials calculations with CP2K and VASP \n\n";
                    }
                }
                else if (qcprog == util::QuantumChemistryProgram::Xtb)
                {
                    std::cout << "Processing xtb g98.out file...\n";
                    LoadFile::loadxtb(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::Vasp)
                {
                    std::cout << "Processing VASP output file...\n";
                    LoadFile::loadvasp(sys);
                }
                // Debug modmass
                // std::cout << "Before modmass:\n";
                // for (int i = 0; i < sys.ncenter; ++i)
                //{
                //    std::cout << "Atom " << i + 1 << " mass: " << std::fixed << std::setprecision(6) << sys.a[i].mass
                //              << " amu\n";
                //}
                util::modmass(sys);
                // Debug modmass
                // std::cout << "After modmass:\n";
                // for (int i = 0; i < sys.ncenter; ++i)
                //{
                //    std::cout << "Atom " << i + 1 << " mass: " << std::fixed << std::setprecision(6) << sys.a[i].mass
                //              << " amu\n";
                //}
                // End Debug modmass
                sys.nelevel = 1;
                sys.elevel  = {0.0};
                sys.edegen  = {sys.spinmult};
                // Ensure degeneracy is positive
                if (!sys.edegen.empty() && sys.edegen[0] <= 0)
                {
                    sys.edegen[0] = 1;
                }
                if (sys.outotm == 1)
                    util::outotmfile(sys);
            }
            else
            {
                // Check if it's a list file (.list or .txt)
                bool is_list_file = (sys.inputfile.find(".list") != std::string::npos) ||
                                    (sys.inputfile.find(".txt") != std::string::npos);
                if (is_list_file)
                {
                    std::cout << "Processing list file...\n";
                    std::vector<std::string> filelist;
                    std::ifstream            listfile(sys.inputfile);
                    if (!listfile.is_open())
                    {
                        throw std::runtime_error("Unable to open list file: " + sys.inputfile);
                    }
                    std::string line;
                    while (std::getline(listfile, line))
                    {
                        // Trim whitespace
                        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch) {
                                       return !std::isspace(ch);
                                   }));
                        line.erase(std::find_if(line.rbegin(),
                                                line.rend(),
                                                [](unsigned char ch) {
                                                    return !std::isspace(ch);
                                                })
                                       .base(),
                                   line.end());
                        if (!line.empty())
                        {
                            filelist.push_back(line);
                        }
                    }
                    listfile.close();
                    if (filelist.empty())
                    {
                        throw std::runtime_error("List file is empty or contains no valid file paths");
                    }
                    // Process batch
                    size_t              nfile = filelist.size();
                    std::vector<double> Elist(nfile), Ulist(nfile), Hlist(nfile), Glist(nfile), Slist(nfile),
                        CVlist(nfile), CPlist(nfile), QVlist(nfile), Qbotlist(nfile);
                    calc::ensemble(sys, filelist, Elist, Ulist, Hlist, Glist, Slist, CVlist, CPlist, QVlist, Qbotlist);
                    // Batch processing complete, skip single file processing
                    return 0;
                }
                else
                {
                    std::cerr << "Error: Unable to identify the quantum chemical program that generated this file.\n";
                    std::cerr << "Supported programs: Gaussian, ORCA, GAMESS-US, NWChem, CP2K, VASP, xTB, and "
                                 "OpenThermo (.otm)\n";
                    std::cerr << "Maybe an old version or newly updated version of supported quantum chemical programs genereted this file \n";
                    std::cerr << "For batch processing, use a list file with .list or .txt extension containing file "
                                 "paths.\n";
                    throw std::runtime_error("Unknown file format");
                }
            }
            if (sys.Eexter != 0.0)
            {
                sys.E = sys.Eexter;
                std::cout << "Note: The electronic energy specified by \"E\" parameter will be used\n";
            }
            else if (sys.E != 0.0)
            {
                std::cout << "Note: The electronic energy extracted from input file will be used\n";
            }

            if (sys.imagreal != 0.0)
            {
                for (int ifreq = 0; ifreq < sys.nfreq; ++ifreq)
                {
                    if (sys.wavenum[ifreq] < 0 && std::abs(sys.wavenum[ifreq]) < sys.imagreal)
                    {
                        sys.wavenum[ifreq] = std::abs(sys.wavenum[ifreq]);
                        std::cout << " Note: Imaginary frequency " << std::fixed << std::setprecision(2)
                                  << sys.wavenum[ifreq] << " cm^-1 has been set to real frequency!\n";
                    }
                }
            }

            sys.totmass = 0.0;
            for (const auto& atom : sys.a)
            {
                sys.totmass += atom.mass;
            }

            calc::calcinertia(sys);

            sys.ilinear = 0;
            for (double in : sys.inert)
            {
                if (in < 0.001)
                {
                    sys.ilinear = 1;
                    break;
                }
            }
            // Debug GAMESS
            // std::cout << "Debug: After Processing, ncenter = " << sys.ncenter << ", a.size() = " << sys.a.size()
            //          << "\n";
            // End Debug GAMESS

            // Symmetry detection
            std::cout << "Number of atoms loaded: " << sys.a.size() << "\n";
            if (sys.a.empty())
            {
                std::cerr << "Error: No atoms loaded from input file!" << "\n";
                std::exit(1);
            }
            symmetry::SymmetryDetector symDetector;
            symDetector.PGnameinit = sys.PGnameinit;
            if (symDetector.PGnameinit == "?")
            {
                // else keep "?" for automatic detection
            }
            symDetector.ncenter = sys.a.size();
            symDetector.a       = sys.a;
            symDetector.a_index.resize(sys.a.size());
            for (size_t i = 0; i < sys.a.size(); ++i)
            {
                symDetector.a_index[i] = i;
            }
            symDetector.detectPG(1);
            sys.rotsym  = symDetector.rotsym;
            sys.PGname = symDetector.PGname;

            std::vector<double> freq(sys.nfreq);
            for (int j = 0; j < sys.nfreq; ++j)
            {
                freq[j] = sys.wavenum[j] * wave2freq;
            }
            sys.freq = freq;

            int nimag = 0;
            for (double f : sys.freq)
            {
                if (f < 0)
                    ++nimag;
            }
            if (nimag > 0)
            {
                std::cout << " Note: There are " << nimag
                          << " imaginary frequencies, they will be ignored in the calculation\n";
            }

            // Print molecular information
            sys.ncenter = sys.a.size();
            std::cout << "\n"
                      << "                      -------- Chemical System Data -------\n"
                      << "                      -------------------------------------\n"
                      << " Electronic energy: " << std::fixed << std::setprecision(8) << std::setw(18) << sys.E
                      << " a.u.\n";
            if (sys.spinmult != 0)
            {
                std::cout << " Spin multiplicity: " << std::setw(3) << sys.spinmult << "\n";
            }
            else
            {
                for (int ie = 0; ie < sys.nelevel; ++ie)
                {
                    std::cout << " Electronic energy level " << ie + 1 << "     E= " << std::fixed
                              << std::setprecision(6) << std::setw(12) << sys.elevel[ie]
                              << " eV     Degeneracy= " << std::setw(3) << sys.edegen[ie] << "\n";
                }
            }
            for (int iatm = 0; iatm < sys.ncenter; ++iatm)
            {
                std::cout << " Atom " << std::setw(5) << iatm + 1 << " (" << ind2name[sys.a[iatm].index]
                          << ")   Mass: " << std::fixed << std::setprecision(6) << std::setw(12) << sys.a[iatm].mass
                          << " amu\n";
                //<< " Index: " << sys.a[iatm].index << "\n";
            }
            // Debug mass of molecule
            // for (int iatm = 0; iatm < sys.ncenter; ++iatm)
            //{
            //    std::cout << "Debug: Atom " << iatm + 1 << " index: " << sys.a[iatm].index << " mass: " << std::fixed
            //              << std::setprecision(6) << sys.a[iatm].mass << " amu\n";
            //}
            // End Debug
            std::cout << " Total mass: " << std::fixed << std::setprecision(6) << std::setw(16) << sys.totmass
                      << " amu\n\n"
                      << " Point group: " << sys.PGname << "\n";
            if (sys.ipmode == 0)
            {
                std::cout << " Rotational symmetry number: " << std::setw(3) << sys.rotsym << "\n";

                // Optional: sort moments for consistent output
                std::array<double, 3> sorted_inert = sys.inert;
                std::sort(sorted_inert.begin(), sorted_inert.end());

                std::cout << " Principal moments of inertia (amu*Bohr^2):\n";
                for (double i : sorted_inert)
                {
                    std::cout << std::fixed << std::setprecision(6) << std::setw(16) << i << "\n";
                }

                double inert_sum = sorted_inert[0] + sorted_inert[1] + sorted_inert[2];
                if (inert_sum < 1e-10)
                {
                    std::cout << "This is a single atom system, rotational constant is zero\n";
                }
                else if (sys.ilinear == 1)
                {
                    // Use the largest (non-zero) moment for linear molecules
                    double largest_inert = sorted_inert[2];  // Since sorted, largest is last
                    double rotcst1       = h / (8.0 * pi * pi * largest_inert * amu2kg * (b2a * 1e-10) * (b2a * 1e-10));
                    std::cout << " Rotational constant (GHz): " << std::fixed << std::setprecision(6) << std::setw(14)
                              << rotcst1 / 1e9 << "\n"
                              << " Rotational temperature (K): " << std::fixed << std::setprecision(6) << std::setw(12)
                              << rotcst1 * h / kb << "\n"
                              << "This is a linear molecule\n";
                }
                else
                {
                    // Non-linear: compute all three
                    for (int i = 0; i < 3; ++i)
                    {
                        rotcst[i] = h / (8.0 * pi * pi * sys.inert[i] * amu2kg * (b2a * 1e-10) * (b2a * 1e-10));
                    }
                    std::cout << " Rotational constants relative to principal axes (GHz):\n";
                    for (double r : rotcst)
                    {
                        std::cout << std::fixed << std::setprecision(6) << std::setw(14) << r / 1e9 << "\n";
                    }
                    std::cout << " Rotational temperatures (K):";
                    for (double r : rotcst)
                    {
                        std::cout << std::fixed << std::setprecision(6) << std::setw(12) << r * h / kb;
                    }
                    std::cout << "\nThis is not a linear molecule\n";
                }
            }
            else
            {
                std::cout << "Rotation information is not shown here since ipmode=1\n";
            }
            if (sys.nfreq > 0)
            {
                std::cout << "\n There are " << sys.nfreq << " frequencies (cm^-1):\n";
                for (int ifreq = 0; ifreq < sys.nfreq; ++ifreq)
                {
                    std::cout << std::fixed << std::setprecision(1) << std::setw(8) << sys.wavenum[ifreq];
                    if ((ifreq + 1) % 9 == 0 || ifreq == sys.nfreq - 1)
                        std::cout << "\n";
                }
            }
            sys.freq.resize(sys.nfreq);
            for (int i = 0; i < sys.nfreq; ++i)
            {
                sys.freq[i] = sys.wavenum[i] * wave2freq;
            }
            if (nimag > 0)
            {
                std::cout << " Note: There are " << nimag
                          << " imaginary frequencies, they will be ignored in the calculation\n";
            }

            // Output thermochemistry results
            if (sys.Tstep == 0.0 && sys.Pstep == 0.0)
            {
                calc::showthermo(sys);
            }
            else
            {
                std::cout << "\nPerforming scan of temperature/pressure...\n";
                double P1 = sys.P, P2 = sys.P, Ps = 1.0;
                double T1 = sys.T, T2 = sys.T, Ts = 1.0;
                if (sys.Tstep != 0.0)
                {
                    T1 = sys.Tlow;
                    T2 = sys.Thigh;
                    Ts = sys.Tstep;
                }
                if (sys.Pstep != 0.0)
                {
                    P1 = sys.Plow;
                    P2 = sys.Phigh;
                    Ps = sys.Pstep;
                }

                // Generate dynamic filenames based on input file
                std::string basename     = get_basename_without_extension(sys.inputfile);
                std::string uhg_filename = basename + ".UHG";
                std::string scq_filename = basename + ".SCq";

                std::ofstream file_UHG(uhg_filename, std::ios::out);
                if (!file_UHG.is_open())
                {
                    throw std::runtime_error("Error: Could not open " + uhg_filename + " for writing");
                }
                file_UHG << "Ucorr, Hcorr and Gcorr are in kcal/mol; U, H and G are in a.u.\n\n"
                         << "     T(K)      P(atm)  Ucorr     Hcorr     Gcorr            U                H           "
                            "     G\n";
                std::ofstream file_SCq(scq_filename, std::ios::out);
                if (!file_SCq.is_open())
                {
                    file_UHG.close();
                    throw std::runtime_error("Error: Could not open " + scq_filename + " for writing");
                }
                file_SCq << "S, CV and CP are in cal/mol/K; q(V=0)/NA and q(bot)/NA are unitless\n\n"
                         << "    T(K)       P(atm)    S         CV        CP        q(V=0)/NA      q(bot)/NA\n";
                if (Ts > 0 && Ps > 0)
                {
                    // Calculate the number of temperature steps
                    // static_cast to int will truncate, so +1 makes it inclusive
                    const int num_step_T = static_cast<int>((T2 - T1) / Ts) + 1;

                    for (int i = 0; i < num_step_T; ++i)
                    {
                        const double T = T1 + i * Ts;  // Calculate T for this iteration

                        // Calculate the number of pressure steps for this temperature
                        const int num_step_P = static_cast<int>((P2 - P1) / Ps) + 1;

                        for (int j = 0; j < num_step_P; ++j)
                        {
                            const double P = P1 + j * Ps;  // Calculate P for this iteration

                            double corrU, corrH, corrG, S, CV, CP, QV, Qbot;
                            calc::calcthermo(sys, T, P, corrU, corrH, corrG, S, CV, CP, QV, Qbot);
                            file_UHG << std::fixed << std::setprecision(3) << std::setw(10) << T << std::setw(10) << P
                                     << std::setprecision(3) << std::setw(10) << corrU / cal2J << std::setw(10)
                                     << corrH / cal2J << std::setw(10) << corrG / cal2J << std::setprecision(6)
                                     << std::setw(17) << corrU / au2kJ_mol + sys.E << std::setw(17)
                                     << corrH / au2kJ_mol + sys.E << std::setw(17) << corrG / au2kJ_mol + sys.E << "\n";
                            file_SCq << std::fixed << std::setprecision(3) << std::setw(10) << T << std::setw(10) << P
                                     << std::setprecision(3) << std::setw(10) << S / cal2J << std::setw(10)
                                     << CV / cal2J << std::setw(10) << CP / cal2J << std::scientific
                                     << std::setprecision(6) << std::setw(16) << QV / NA << std::setw(16) << Qbot / NA
                                     << "\n";
                        }
                    }
                }
                file_UHG.close();
                file_SCq.close();
                std::cout << "\n Congratulation! Thermochemical properties at various temperatures/pressures were calculated " << "\n"
                          << " " << "All data were exported to " << uhg_filename << " and " << scq_filename << "\n" 
                          << " " << uhg_filename << " contains thermal correction to U, H and G, and sum of electronic energy and corresponding corrections\n"
                          << " " << scq_filename << " contains S, CV, CP, q(V=0) and q(bot)\n";
            }
        }

        if (narg == 0)
        {
            std::cout << "\nRunning finished! Press ENTER to exit\n";
            std::cin.get();
        }
        std::cout << "\n";
        auto        now      = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::cout << "Calculation completed at: " << std::ctime(&now_time);
        std::cout << "\n"
                  << "                    ---------- Happy calculation ----------" << "\n"
                  << "                    ---- OpenThermo normally terminated ---" << "\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nError: " << e.what() << "\n";
        std::cerr << "Program terminated due to an error." << "\n"
                  << "\n";
        if (argc == 1)
        {  // No arguments provided
            std::cout << "Press ENTER to exit" << "\n";
            std::cin.get();
        }
        return 1;
    }
}