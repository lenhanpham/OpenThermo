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
#include "omp_config.h"
#include "version.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
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
        // OpenMP thread detection and configuration is deferred until after
        // argument parsing, so that -omp-threads can be read first.

        SystemData            sys;                       // Main system data structure
        std::array<double, 3> rotcst = {0.0, 0.0, 0.0};  // Rotational constants

        // Handle early-exit options before any other processing
        if (argc > 1)
        {
            std::string first_arg = argv[1];
            if (first_arg == "--help")
            {
                HelpUtils::print_help(argv[0]);
                return 0;
            }
            else if (first_arg == "--version" || first_arg == "-v")
            {
                HelpUtils::print_version();
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

        // Print program banner
        HelpUtils::print_version();

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

        // --- Apply method-dependent Bav ---
        // Only HeadGordon supports the -bav option (grimme or qchem).
        // Grimme/Minenkov always use Grimme's Bav (1e-44 kg m^2).
        if (sys.lowVibTreatment == LowVibTreatment::HeadGordon)
        {
            if (!sys.bavUserOverride)
            {
                // HeadGordon defaults to Q-Chem's Bav
                sys.bavPreset = BavPreset::QChem;
                sys.Bav       = bavPresetValue(BavPreset::QChem);
            }
        }
        else
        {
            // Grimme/Minenkov: always use Grimme's Bav; warn if user tried to override
            if (sys.bavUserOverride && sys.bavPreset != BavPreset::Grimme)
            {
                std::cerr << "Warning: -bav option is only applicable to HeadGordon method. "
                          << "Ignoring -bav " << bavPresetName(sys.bavPreset)
                          << "; using grimme (1e-44 kg m^2).\n";
            }
            sys.bavPreset = BavPreset::Grimme;
            sys.Bav       = bavPresetValue(BavPreset::Grimme);
        }

        // --- OpenMP thread detection and configuration ---
        sys.exec.physical_cores_detected = detect_physical_cores();
        sys.exec.scheduler_cpus_detected = detect_scheduler_cpus();
        std::string thread_notification = validate_thread_count(
            sys.exec.omp_threads_requested,
            sys.exec.physical_cores_detected,
            sys.exec.scheduler_cpus_detected,
            sys.exec.omp_threads_actual,
            sys.exec.omp_user_override
        );
        configure_openmp(sys.exec.omp_threads_actual);

        if (sys.prtlevel >= 1 && !thread_notification.empty())
        {
            std::cout << "\n" << thread_notification << "\n";
        }

        // If prtlevel=3, auto-enable per-mode vibration output unless user explicitly set prtvib
        if (sys.prtlevel >= 3 && sys.prtvib == 0)
        {
            sys.prtvib = 1;
        }

        // Print running parameters
        if (sys.prtlevel >= 1)
        {
        std::cout << "\n                   --- Summary of Current Parameters ---\n\nRunning parameters:\n";
        std::cout << " Print level: " << sys.prtlevel << " (0=minimal, 1=default, 2=verbose, 3=full)\n";
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
        else if (sys.lowVibTreatment == LowVibTreatment::HeadGordon)
        {
            std::cout << " Low frequencies treatment: Head-Gordon's interpolation for energy";
            if (sys.hgEntropy)
                std::cout << " and entropy";
            std::cout << "\n";
        }
        if (sys.lowVibTreatment == LowVibTreatment::Grimme || sys.lowVibTreatment == LowVibTreatment::Minenkov
            || sys.lowVibTreatment == LowVibTreatment::HeadGordon)
        {
            std::cout << " Vibrational frequency threshold used in the interpolation is " << std::fixed
                      << std::setprecision(2) << sys.intpvib << " cm^-1\n";
        }
        if (sys.lowVibTreatment == LowVibTreatment::HeadGordon)
        {
            std::cout << " Average moment of inertia (Bav): " << bavPresetName(sys.bavPreset)
                      << " (" << std::scientific << std::setprecision(2) << sys.Bav << " kg m^2)\n"
                      << std::fixed;
        }
        if (sys.imagreal != 0.0)
        {
            std::cout << " Imaginary frequencies with norm < " << std::fixed << std::setprecision(2) << sys.imagreal
                      << " cm^-1 will be treated as real frequencies\n";
        }
        } // end if (sys.prtlevel >= 1) for parameter summary

        // Load input file
        if (narg >= 1)
        {
            sys.inputfile = argv[1];
        }
        if (sys.inputfile.empty())
        {
            std::cout << "\nInput file path, e.g. D:\\your_dir\\your_calc.log\n"
                      << " OpenThermo supports Gaussian, ORCA, GAMESS-US, NWChem, CP2K, VASP, and Q-Chem "
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
        if (sys.prtlevel >= 1)
        {
            std::cout << "                      -------- End of Summary --------\n";
            std::cout << "\n";
            std::cout << "OpenThermo started to process " << sys.inputfile << " at "
                      << std::ctime(&start_now_time);  // << "\n";
        }

        // Process input file
        if (sys.inputfile.find(".otm") != std::string::npos)
        {
            if (sys.prtlevel >= 2)
            {
                std::cout << "\n Processing data from " << sys.inputfile << "\n";
                std::cout << " Atomic masses used: Read from OTM file\n";
            }
            LoadFile::loadotm(sys);
        }
        else
        {
            auto qcprog = util::deterprog(sys);
            if (qcprog != util::QuantumChemistryProgram::Unknown)
            {
                if (sys.prtlevel >= 2)
                {
                std::cout << "\n";
                if (sys.massmod == 1)
                    std::cout << " Atomic masses used: Element\n";
                if (sys.massmod == 2)
                    std::cout << " Atomic masses used: Most abundant isotope\n";
                if (sys.massmod == 3)
                    std::cout << " Atomic masses used: Read from quantum chemical output\n";
                }
                if (qcprog == util::QuantumChemistryProgram::Gaussian)
                {
                    if (sys.prtlevel >= 2) std::cout << "Processing Gaussian output file...\n";
                    LoadFile::loadgau(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::Orca)
                {
                    if (sys.prtlevel >= 2) std::cout << "Processing ORCA output file...\n";
                    LoadFile::loadorca(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::Gamess)
                {
                    if (sys.prtlevel >= 2) std::cout << "Processing GAMESS-US output file...\n";
                    LoadFile::loadgms(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::Nwchem)
                {
                    if (sys.prtlevel >= 2) std::cout << "Processing NWChem output file...\n";
                    LoadFile::loadnw(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::Cp2k)
                {
                    if (sys.prtlevel >= 2) std::cout << "Processing CP2K output file...\n";
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
                    if (sys.prtlevel >= 2) std::cout << "Processing xtb g98.out file...\n";
                    LoadFile::loadxtb(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::Vasp)
                {
                    if (sys.prtlevel >= 2) std::cout << "Processing VASP output file...\n";
                    LoadFile::loadvasp(sys);
                }
                else if (qcprog == util::QuantumChemistryProgram::QChem)
                {
                    if (sys.prtlevel >= 2) std::cout << "Processing Q-Chem output file...\n";
                    LoadFile::loadqchem(sys);
                }
                util::modmass(sys);
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
                    std::cerr << "Supported programs: Gaussian, ORCA, GAMESS-US, NWChem, CP2K, VASP, Q-Chem, xTB, and "
                                 "OpenThermo (.otm)\n";
                    std::cerr << "Maybe an old version or newly updated version of supported quantum chemical programs genereted this file \n";
                    std::cerr << "For batch processing, use a list file with .list or .txt extension containing file "
                                 "paths.\n";
                    throw std::runtime_error("Unknown file format");
                }
            }
        }
        if (sys.Eexter != 0.0)
            {
                sys.E = sys.Eexter;
                if (sys.prtlevel >= 1)
                    std::cout << "Note: The electronic energy specified by \"E\" parameter will be used\n";
            }
            else if (sys.E != 0.0)
            {
                if (sys.prtlevel >= 2)
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
            // Symmetry detection
            if (sys.prtlevel >= 2)
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
            symDetector.detectPG((sys.prtlevel >= 2) ? 1 : 0);
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

            if (sys.prtlevel >= 1)
            {
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
                        std::cout << " Electronic energy level " << ie + 1 << "     E = " << std::fixed
                                  << std::setprecision(6) << std::setw(12) << sys.elevel[ie]
                                  << " eV     Degeneracy = " << std::setw(3) << sys.edegen[ie] << "\n";
                    }
                }
            }

            if (sys.prtlevel >= 2)
            {
                // Level 2+: full per-atom listing
                for (int iatm = 0; iatm < sys.ncenter; ++iatm)
                {
                    std::cout << " Atom " << std::setw(5) << iatm + 1 << " (" << ind2name[sys.a[iatm].index]
                              << ")   Mass: " << std::fixed << std::setprecision(6) << std::setw(12) << sys.a[iatm].mass
                              << " amu\n";
                }
                std::cout << " Total mass: " << std::fixed << std::setprecision(6) << std::setw(16) << sys.totmass
                          << " amu\n\n";
            }
            else if (sys.prtlevel == 1)
            {
                // Level 1: compact atom count summary
                std::map<int, int> elem_count;
                for (int iatm = 0; iatm < sys.ncenter; ++iatm)
                {
                    elem_count[sys.a[iatm].index]++;
                }
                std::cout << " Atoms: " << sys.ncenter << " (";
                bool first = true;
                for (const auto& [idx, count] : elem_count)
                {
                    if (!first) std::cout << ", ";
                    std::cout << count << " " << ind2name[idx];
                    first = false;
                }
                std::cout << ")  Total mass: " << std::fixed << std::setprecision(6) << sys.totmass << " amu\n";
            }

            if (sys.prtlevel >= 1)
            {
                std::cout << " Point group: " << sys.PGname;
                if (sys.ipmode == 0)
                    std::cout << "   Rotational symmetry number: " << std::setw(3) << sys.rotsym;
                std::cout << "\n";
            }

            if (sys.prtlevel >= 2 && sys.ipmode == 0)
            {
                // Level 2+: moments of inertia and rotational constants
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
                    double largest_inert = sorted_inert[2];
                    double rotcst1       = h / (8.0 * pi * pi * largest_inert * amu2kg * (b2a * 1e-10) * (b2a * 1e-10));
                    std::cout << " Rotational constant (GHz): " << std::fixed << std::setprecision(6) << std::setw(14)
                              << rotcst1 / 1e9 << "\n"
                              << " Rotational temperature (K): " << std::fixed << std::setprecision(6) << std::setw(12)
                              << rotcst1 * h / kb << "\n"
                              << "This is a linear molecule\n";
                }
                else
                {
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
            else if (sys.prtlevel >= 2 && sys.ipmode == 1)
            {
                std::cout << "Rotation information is not shown here since ipmode=1\n";
            }

            if (sys.nfreq > 0)
            {
                if (sys.prtlevel >= 2)
                {
                    // Level 2+: full frequency listing
                    std::cout << "\n There are " << sys.nfreq << " frequencies (cm^-1):\n";
                    for (int ifreq = 0; ifreq < sys.nfreq; ++ifreq)
                    {
                        std::cout << std::fixed << std::setprecision(1) << std::setw(8) << sys.wavenum[ifreq];
                        if ((ifreq + 1) % 9 == 0 || ifreq == sys.nfreq - 1)
                            std::cout << "\n";
                    }
                }
                else if (sys.prtlevel == 1)
                {
                    // Level 1: compact frequency count + range
                    double wmin = sys.wavenum[0], wmax = sys.wavenum[0];
                    for (int ifreq = 1; ifreq < sys.nfreq; ++ifreq)
                    {
                        if (sys.wavenum[ifreq] < wmin) wmin = sys.wavenum[ifreq];
                        if (sys.wavenum[ifreq] > wmax) wmax = sys.wavenum[ifreq];
                    }
                    std::cout << " Frequencies: " << sys.nfreq << " (range: " << std::fixed << std::setprecision(1)
                              << wmin << " -- " << wmax << " cm^-1)\n";
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
                // Single T/P point: always use inner strategy if beneficial
                sys.exec.omp_strategy = static_cast<int>(
                    select_strategy(1, sys.nfreq, sys.exec.omp_threads_actual));
                if (sys.prtlevel >= 2)
                {
                    std::cout << strategy_description(
                        static_cast<OMPStrategy>(sys.exec.omp_strategy), 1, sys.nfreq) << "\n";
                }
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
                    const int num_step_T = static_cast<int>((T2 - T1) / Ts) + 1;
                    const int num_step_P = static_cast<int>((P2 - P1) / Ps) + 1;
                    const int total_points = num_step_T * num_step_P;

                    // Auto-select parallelization strategy
                    OMPStrategy strategy = select_strategy(total_points, sys.nfreq, sys.exec.omp_threads_actual);
                    sys.exec.omp_strategy = static_cast<int>(strategy);

                    if (sys.prtlevel >= 2)
                    {
                        std::cout << strategy_description(strategy, total_points, sys.nfreq) << "\n";
                    }

                    struct ScanResult {
                        double T, P;
                        double corrU, corrH, corrG, S, CV, CP, QV, Qbot;
                    };
                    std::vector<ScanResult> results(total_points);

                    if (strategy == OMPStrategy::Outer)
                    {
                        // Outer strategy: parallelize T/P scan, calcthermo runs serially
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int idx = 0; idx < total_points; ++idx)
                        {
                            int i = idx / num_step_P;
                            int j = idx % num_step_P;
                            auto& r = results[idx];
                            r.T = T1 + i * Ts;
                            r.P = P1 + j * Ps;
                            calc::calcthermo(sys, r.T, r.P, r.corrU, r.corrH, r.corrG,
                                             r.S, r.CV, r.CP, r.QV, r.Qbot);
                        }
                    }
                    else
                    {
                        // Inner strategy: serial T/P scan, vibrational loop parallelized in calc.cpp
                        for (int idx = 0; idx < total_points; ++idx)
                        {
                            int i = idx / num_step_P;
                            int j = idx % num_step_P;
                            auto& r = results[idx];
                            r.T = T1 + i * Ts;
                            r.P = P1 + j * Ps;
                            calc::calcthermo(sys, r.T, r.P, r.corrU, r.corrH, r.corrG,
                                             r.S, r.CV, r.CP, r.QV, r.Qbot);
                        }
                    }

                    for (int idx = 0; idx < total_points; ++idx)
                    {
                        const auto& r = results[idx];
                        file_UHG << std::fixed << std::setprecision(3) << std::setw(10) << r.T << std::setw(10) << r.P
                                 << std::setprecision(3) << std::setw(10) << r.corrU / cal2J << std::setw(10)
                                 << r.corrH / cal2J << std::setw(10) << r.corrG / cal2J << std::setprecision(6)
                                 << std::setw(17) << r.corrU / au2kJ_mol + sys.E << std::setw(17)
                                 << r.corrH / au2kJ_mol + sys.E << std::setw(17) << r.corrG / au2kJ_mol + sys.E << "\n";
                        file_SCq << std::fixed << std::setprecision(3) << std::setw(10) << r.T << std::setw(10) << r.P
                                 << std::setprecision(3) << std::setw(10) << r.S / cal2J << std::setw(10)
                                 << r.CV / cal2J << std::setw(10) << r.CP / cal2J << std::scientific
                                 << std::setprecision(6) << std::setw(16) << r.QV / NA << std::setw(16) << r.Qbot / NA
                                 << "\n";
                    }
                }
                file_UHG.close();
                file_SCq.close();
                std::cout << "\n Congratulation! Thermochemical properties at various temperatures/pressures were calculated " << "\n"
                          << " " << "All data were exported to " << uhg_filename << " and " << scq_filename << "\n" 
                          << " " << uhg_filename << " contains thermal correction to U, H and G, and sum of electronic energy and corresponding corrections\n"
                          << " " << scq_filename << " contains S, CV, CP, q(V=0) and q(bot)\n";
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