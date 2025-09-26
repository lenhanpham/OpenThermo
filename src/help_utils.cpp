/**
 * @file help_utils.cpp
 * @brief Implementation of help utility functions for OpenThermo
 * @author Le Nhan Pham
 * @date 2025
 *
 * This file contains implementations for help-related utility functions
 * used throughout the OpenThermo molecular thermochemistry program.
 */

#include "help_utils.h"
#include <iostream>
#include <map>
#include <string>

namespace HelpUtils
{

    void print_help(const std::string& program_name)
    {
        std::cout << "OpenThermo: A Comprehensive C++ Program for Calculation of Thermochemical Properties\n"
                  << "Version 0.001.1\n"
                  << "Developer: Le Nhan Pham\n\n";
        std::cout << "Usage: " << program_name << " [input_file] [options]\n\n";
        std::cout << "Description:\n";
        std::cout << "  OpenThermo calculates thermochemical properties from quantum chemistry output files.\n";
        std::cout << "  It supports various input formats and provides comprehensive thermodynamic analysis\n";
        std::cout
            << "  including Gibbs free energy, enthalpy, entropy, heat capacity, and vibrational corrections.\n\n";
        std::cout << "Input Files:\n";
        std::cout << "  input_file    Path to input file (.otm format or quantum chemistry output)\n";
        std::cout << "                Supported formats: Gaussian, ORCA, GAMESS-US, NWChem, CP2K, VASP\n";
        std::cout << "                If no file specified, program will prompt for input\n\n";
        std::cout << "Options:\n";
        std::cout << "  -E <value>           Electronic energy in a.u. (overrides file value)\n";
        std::cout << "  -T <T>               Temperature in K (default: 298.15)\n";
        std::cout << "  -T <T1 T2 step>      Temperature scan from T1 to T2 with step\n";
        std::cout << "  -P <P>               Pressure in atm (default: 1.0)\n";
        std::cout << "  -P <P1 P2 step>      Pressure scan from P1 to P2 with step\n";
        std::cout << "  -sclZPE <factor>     Scale factor for ZPE frequencies (default: 1.0)\n";
        std::cout << "  -sclheat <factor>    Scale factor for thermal energy frequencies (default: 1.0)\n";
        std::cout << "  -sclS <factor>       Scale factor for entropy frequencies (default: 1.0)\n";
        std::cout << "  -sclCV <factor>      Scale factor for heat capacity frequencies (default: 1.0)\n";
        std::cout << "  -lowvibmeth <mode>     Low frequency treatment: 0/Harmonic, 1/Truhlar, 2/Grimme, 3/Minenkov\n";
        std::cout << "  -ravib <value>       Raising value for low frequencies in cm^-1 (default: 100.0)\n";
        std::cout << "  -ipmode <mode>       Calculation mode: 0=gas phase, 1=condensed phase\n";
        std::cout << "  -imagreal <value>    Treat imaginary freq < value as real (default: 0.0)\n";
        std::cout << "  -conc <string>       Concentration string for phase correction\n";
        std::cout << "  -massmod <type>      Default mass type: 1=element, 2=most abundant isotope, 3=file\n";
        std::cout << "  -PGlabel <label>     Force point group symmetry\n";
        std::cout << "  -prtvib <mode>       Print vibration contributions: 0=no, 1=yes, -1=to file\n";
        std::cout << "  -outotm <mode>       Output .otm file: 0=no, 1=yes\n";
        std::cout << "  -noset               Don't load settings from settings.ini\n";
        std::cout << "  --help               Show this help message\n";
        std::cout << "  --create-config      Create a default settings.ini file\n";
        std::cout << "  --help-input         Show input file format help\n";
        std::cout << "  --help-output        Show output format help\n";
        std::cout << "  --help-settings      Show settings file help\n";
        std::cout << "  --help-<option>      Show help for specific option (e.g., --help-T)\n\n";
        std::cout << "Settings File:\n";
        std::cout << "  Parameters can also be set in settings.ini file in current directory\n";
        std::cout << "  or in the directory specified by OPENTHERMOPATH environment variable\n\n";
        std::cout << "Output Files:\n";
        std::cout << "  <basename>.UHG       Thermodynamic quantities vs T/P (if scan performed)\n";
        std::cout << "  <basename>.SCq       Entropy, heat capacities, partition functions\n";
        std::cout << "  <basename>.vibcon    Individual vibration contributions (if requested)\n";
        std::cout << "  *.otm                OpenThermo format file (if -outotm 1)\n\n";
        std::cout << "Examples:\n";
        std::cout << "  " << program_name << " molecule.log\n";
        std::cout << "  " << program_name << " molecule.otm -T 300 -P 2.0\n";
        std::cout << "  " << program_name << " molecule.out -T 273 373 10 -lowvibmeth 2\n";
        std::cout << "  " << program_name << " --help-input\n";
        std::cout << "  " << program_name << " --help-T\n\n";
        std::cout << "For more detailed help on specific topics, use --help-<topic>\n";
    }

    void print_option_help(const std::string& option, const std::string& program_name)
    {
        std::cout << program_name << ": Help for option: -" << option << "\n\n";

        std::map<std::string, std::string> option_descriptions = {
            {"E",
             "Electronic Energy\n"
             "  -E <value>\n"
             "  Sets the electronic energy in atomic units (a.u.)\n"
             "  This overrides any electronic energy found in the input file\n"
             "  Example: -E -76.384729\n"
             "  Default: Use value from input file"},
            {"T",
             "Temperature\n"
             "  -T <T>          Single temperature in Kelvin\n"
             "  -T <T1 T2 step> Temperature scan from T1 to T2 with given step\n"
             "  Used for calculating temperature-dependent thermodynamic properties\n"
             "  Examples:\n"
             "    -T 298.15     (single temperature)\n"
             "    -T 200 400 25 (scan from 200K to 400K in 25K steps)\n"
             "  Default: 298.15 K"},
            {"P",
             "Pressure\n"
             "  -P <P>          Single pressure in atm\n"
             "  -P <P1 P2 step> Pressure scan from P1 to P2 with given step\n"
             "  Affects the calculation of Gibbs free energy through -RT ln(Q/N)\n"
             "  Examples:\n"
             "    -P 1.0        (standard pressure)\n"
             "    -P 0.5 2.0 0.2 (scan from 0.5 to 2.0 atm in 0.2 atm steps)\n"
             "  Default: 1.0 atm"},
            {"sclZPE",
             "ZPE Frequency Scaling\n"
             "  -sclZPE <factor>\n"
             "  Scaling factor applied to vibrational frequencies for ZPE calculation\n"
             "  ZPE = 0.5 * h * c * sum(sclZPE * ν_i)\n"
             "  Example: -sclZPE 0.98\n"
             "  Default: 1.0"},
            {"sclheat",
             "Thermal Energy Frequency Scaling\n"
             "  -sclheat <factor>\n"
             "  Scaling factor for vibrational frequencies in thermal energy calculation\n"
             "  U_vib = sum(hν/(exp(hν/kT)-1))\n"
             "  Example: -sclheat 0.99\n"
             "  Default: 1.0"},
            {"sclS",
             "Entropy Frequency Scaling\n"
             "  -sclS <factor>\n"
             "  Scaling factor for vibrational frequencies in entropy calculation\n"
             "  S_vib = R * sum(ln(1-exp(-hν/kT)) + hν/(kT(exp(hν/kT)-1)))\n"
             "  Example: -sclS 0.99\n"
             "  Default: 1.0"},
            {"sclCV",
             "Heat Capacity Frequency Scaling\n"
             "  -sclCV <factor>\n"
             "  Scaling factor for vibrational frequencies in heat capacity calculation\n"
             "  CV_vib = R * sum((hν/kT)^2 * exp(hν/kT) / (exp(hν/kT)-1)^2)\n"
             "  Example: -sclCV 0.99\n"
             "  Default: 1.0"},
            {"lowvibmeth",
             "Low Frequency Treatment\n"
             "  -lowvibmeth <mode>\n"
             "  Method for handling low vibrational frequencies:\n"
             "    0/Harmonic: Harmonic approximation (no special treatment)\n"
             "    1/Truhlar: Raise frequencies below threshold to ravib value\n"
             "    2/Grimme: Grimme's interpolation for entropy\n"
             "    3/Minenkov: Minenkov's interpolation for entropy and internal energy\n"
             "  Example: -lowvibmeth 2 or -lowvibmeth Grimme\n"
             "  Default: Grimme"},
            {"ravib",
             "Low Frequency Raising Value\n"
             "  -ravib <value>\n"
             "  Frequency value (cm^-1) to which low frequencies are raised when lowvibmeth=1\n"
             "  Example: -ravib 50.0\n"
             "  Default: 100.0"},
            {"ipmode",
             "Calculation Mode\n"
             "  -ipmode <mode>\n"
             "  Calculation mode:\n"
             "    0: Gas phase (include translational/rotational degrees of freedom)\n"
             "    1: Condensed phase (remove translational/rotational contributions)\n"
             "  Example: -ipmode 1\n"
             "  Default: 0"},
            {"imagreal",
             "Imaginary Frequency Treatment\n"
             "  -imagreal <value>\n"
             "  Treat imaginary frequencies with |ν| < value as real frequencies\n"
             "  Useful for transition states with small imaginary frequencies\n"
             "  Example: -imagreal 50.0\n"
             "  Default: 0.0 (no treatment)"},
            {"conc",
             "Concentration\n"
             "  -conc <value>\n"
             "  Concentration for phase correction (mol/L)\n"
             "  Format: numeric value representing concentration in mol/L\n"
             "  Example: -conc 2.5\n"
             "  Default: None"},
            {"massmod",
             "Default Atomic Masses\n"
             "  -massmod <type>\n"
             "  Type of atomic masses to use:\n"
             "    1: Element average mass\n"
             "    2: Most abundant isotope mass\n"
             "    3: Masses from input file\n"
             "  Example: -massmod 2\n"
             "  Default: 1"},
            {"PGlabel",
             "Point Group\n"
             "  -PGlabel <label>\n"
             "  Force specific point group for symmetry analysis\n"
             "  Overrides automatic symmetry detection\n"
             "  Example: -PGlabel C2v\n"
             "  Default: Auto-detect"},
            {"prtvib",
             "Print Vibration Contributions\n"
             "  -prtvib <mode>\n"
             "  Print individual vibrational mode contributions:\n"
             "    0: No (default)\n"
             "    1: Yes, to screen\n"
             "   -1: Yes, to <basename>.vibcon file\n"
             "  Example: -prtvib 1\n"
             "  Default: 0"},
            {"outotm",
             "Output OpenThermo File\n"
             "  -outotm <mode>\n"
             "  Generate .otm format file:\n"
             "    0: No (default)\n"
             "    1: Yes\n"
             "  Example: -outotm 1\n"
             "  Default: 0"},
            {"noset",
             "Skip Settings File\n"
             "  -noset\n"
             "  Don't load parameters from settings.ini file\n"
             "  Use this to ignore settings file and use defaults\n"
             "  Example: -noset\n"
             "  Default: Load settings if available"},
            {"create-config",
             "Create Default Settings File\n"
             "  --create-config\n"
             "  Generate a template settings.ini file with all default parameters\n"
             "  This creates a settings.ini file in the current directory that you can\n"
             "  edit to customize your default calculation parameters.\n"
             "  Example: ./OpenThermo --create-config\n"
             "  Note: This command exits immediately after creating the file"}};

        auto it = option_descriptions.find(option);
        if (it != option_descriptions.end())
        {
            std::cout << it->second << "\n\n";
        }
        else
        {
            std::cout << "Unknown option: -" << option << "\n";
            std::cout << "Use --help to see all available options\n\n";
        }
    }

    void print_input_help()
    {
        std::cout << "OpenThermo Input File Formats\n\n";
        std::cout << "Supported Input Formats:\n\n";
        std::cout << "1. OpenThermo Format (.otm):\n";
        std::cout << "   Native format containing all necessary molecular data\n";
        std::cout << "   Sections:\n";
        std::cout << "     *E          Electronic energy in a.u.\n";
        std::cout << "     *wavenum    Vibrational frequencies in cm^-1\n";
        std::cout << "     *atoms      Atomic coordinates and masses\n";
        std::cout << "     *elevel     Electronic energy levels and degeneracies\n\n";
        std::cout << "2. Quantum Chemistry Output Files:\n";
        std::cout << "   Gaussian (.log, .out):\n";
        std::cout << "     - Frequency analysis output\n";
        std::cout << "     - Extracts geometry, frequencies, and electronic energy\n\n";
        std::cout << "   ORCA (.out):\n";
        std::cout << "     - Frequency calculation output\n";
        std::cout << "     - Supports various ORCA output formats\n\n";
        std::cout << "   GAMESS-US (.log):\n";
        std::cout << "     - Hessian/frequency analysis output\n\n";
        std::cout << "   NWChem (.out):\n";
        std::cout << "     - Frequency analysis output\n\n";
        std::cout << "   CP2K (.out):\n";
        std::cout << "     - Vibrational analysis output\n";
        std::cout << "     - Supports both molecular and periodic systems\n\n";
        std::cout << "   VASP (OUTCAR):\n";
        std::cout << "     - Vibrational analysis output\n";
        std::cout << "     - Supports both molecular and periodic systems\n";  
        std::cout << "     - Both CONTCAR and OUTCAR are needed \n\n";        
        //std::cout << "   xtb (g98.out):\n";
        //std::cout << "     - xtb frequency analysis output in Gaussian format\n\n";
        std::cout << "3. List Files (text file):\n";
        std::cout << "   Text file containing paths to multiple input files (one per line)\n";
        std::cout << "   Useful for batch processing of multiple molecules\n\n";
        std::cout << "File Detection:\n";
        std::cout << "   The program automatically detects the input format based on file content\n";
        std::cout << "   No need to specify format explicitly\n\n";
        std::cout << "Requirements:\n";
        std::cout << "   - Geometry (atomic coordinates)\n";
        std::cout << "   - Vibrational frequencies (from frequency analysis)\n";
        std::cout << "   - Electronic energy (SCF or similar)\n";
        std::cout << "   - For multi-level calculations: electronic energy levels\n\n";
    }

    void print_output_help()
    {
        std::cout << "OpenThermo Output Information\n\n";
        std::cout << "Console Output:\n";
        std::cout << "  - Program version and developer information\n";
        std::cout << "  - Running parameters (temperature, pressure, scaling factors, etc.)\n";
        std::cout << "  - Molecular information (atoms, masses, point group, moments of inertia)\n";
        std::cout << "  - Vibrational frequencies\n";
        std::cout << "  - Thermodynamic properties at specified conditions\n\n";
        std::cout << "Output Files:\n\n";
        std::cout << "1. <basename>.UHG (Temperature/Pressure Scan):\n";
        std::cout << "   Contains thermal corrections and total energies\n";
        std::cout << "   Columns: T(K), P(atm), U_corr, H_corr, G_corr, U_total, H_total, G_total\n";
        std::cout << "   Units: Thermal corrections in kcal/mol, energies in a.u.\n\n";
        std::cout << "2. <basename>.SCq (Temperature/Pressure Scan):\n";
        std::cout << "   Contains entropy, heat capacities, and partition functions\n";
        std::cout << "   Columns: T(K), P(atm), S, CV, CP, q(V=0)/NA, q(bot)/NA\n";
        std::cout << "   Units: S, CV, CP in cal/mol/K, partition functions dimensionless\n\n";
        std::cout << "3. <basename>.vibcon (Vibration Contributions):\n";
        std::cout << "   Individual contributions of each vibrational mode\n";
        std::cout << "   Generated when prtvib = 1 or -1\n\n";
        std::cout << "4. *.otm (OpenThermo Format):\n";
        std::cout << "   Native format file containing all molecular data\n";
        std::cout << "   Generated when outotm = 1\n";
        std::cout << "   Can be used as input for subsequent calculations\n\n";
        std::cout << "Thermodynamic Properties:\n";
        std::cout << "  - Electronic energy (E)\n";
        std::cout << "  - Zero-point energy (ZPE)\n";
        std::cout << "  - Thermal corrections (U_corr, H_corr, G_corr)\n";
        std::cout << "  - Total energies (U, H, G)\n";
        std::cout << "  - Entropy (S)\n";
        std::cout << "  - Heat capacities (CV, CP)\n";
        std::cout << "  - Partition functions (translational, rotational, vibrational, electronic)\n\n";
        std::cout << "Units:\n";
        std::cout << "  - Energy: a.u. (atomic units) or kJ/mol\n";
        std::cout << "  - Temperature: K (Kelvin)\n";
        std::cout << "  - Pressure: atm (atmospheres)\n";
        std::cout << "  - Frequency: cm^-1\n";
        std::cout << "  - Mass: amu (atomic mass units)\n";
        std::cout << "  - Distance: Angstrom\n\n";
    }

    void print_settings_help()
    {
        std::cout << "OpenThermo Settings File (settings.ini)\n\n";
        std::cout << "Location:\n";
        std::cout << "  - Current working directory: ./settings.ini\n";
        std::cout << "  - Environment variable: $OPENTHERMOPATH/settings.ini\n";
        std::cout << "  - If both exist, current directory takes precedence\n\n";
        std::cout << "Format:\n";
        std::cout << "  # Lines starting with # are comments\n";
        std::cout << "  parameter = value\n";
        std::cout << "  # Values can be quoted: parameter = \"value with spaces\"\n\n";
        std::cout << "Available Parameters:\n\n";
        std::cout << "# Electronic energy\n";
        std::cout << "E = -76.384729\n\n";
        std::cout << "# Temperature (single value or T1 T2 step)\n";
        std::cout << "T = 298.15\n";
        std::cout << "# T = 200.0 400.0 25.0\n\n";
        std::cout << "# Pressure (single value or P1 P2 step)\n";
        std::cout << "P = 1.0\n";
        std::cout << "# P = 0.5 2.0 0.2\n\n";
        std::cout << "# Frequency scaling factors\n";
        std::cout << "sclZPE = 1.0\n";
        std::cout << "sclheat = 1.0\n";
        std::cout << "sclS = 1.0\n";
        std::cout << "sclCV = 1.0\n\n";
        std::cout << "# Low frequency treatment\n";
        std::cout << "lowvibmeth = Grimme\n";
        std::cout << "ravib = 100.0\n";
        std::cout << "intpvib = 100.0\n\n";
        std::cout << "# Calculation mode and options\n";
        std::cout << "ipmode = 0\n";
        std::cout << "imagreal = 0.0\n";
        std::cout << "conc = 1.0\n";
        std::cout << "massmod = 1\n";
        std::cout << "PGlabel = \"?\"\n\n";
        std::cout << "# Output options\n";
        std::cout << "prtvib = 0\n";
        std::cout << "outotm = 0\n\n";
        std::cout << "# VASP energy selection\n";
        std::cout << "extrape = false \n\n";
        std::cout << "# Mass modifications (optional section)\n";
        std::cout << "# modmass should match exactly the index order of atoms in quantum chemical outputs\n";
        std::cout << "# 1 1.007825  # Atom index, element, mass\n";
        std::cout << "# 2 12.0      # Use specific isotope mass\n\n";
        std::cout << "Notes:\n";
        std::cout << "  - Parameters set in command line override settings file values\n";
        std::cout << "  - Use -noset to skip loading settings file\n";
        std::cout << "  - Mass modifications section is optional and allows per-atom mass changes\n";
        std::cout << "  - Empty or missing settings file uses program defaults\n\n";
    }

}  // namespace HelpUtils