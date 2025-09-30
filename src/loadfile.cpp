/**
 * @file loadfile.cpp
 * @brief Implementation of input file parsing and loading functions
 * @author Le Nhan Pham
 * @date 2025
 *
 * This file contains the implementation of functions for reading and parsing
 * various computational chemistry output file formats including Gaussian,
 * ORCA, NWChem, GAMESS, CP2K, and xTB outputs, as well as OpenThermo's native
 * .otm format.
 */


#include "loadfile.h"
#include "chemsys.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


// Implemetation of LoadFile class
/**
 * @brief Locate a specific label/string in an input file stream
 *
 * Searches for a given label string in the file
 * stream and positions the stream
 * at the matching line, optionally skipping additional lines after the match.
 *
 *
 * @param file Input file stream to search
 * @param label String pattern to locate in the file
 * @param skip Number of
 * lines to skip after finding the label (default: 0)
 * @return true if label found and stream positioned successfully,
 * false otherwise
 *
 * @note Removes whitespace and carriage returns from lines before comparison
 * @note Positions
 * stream at the beginning of the matching line
 * @note Used extensively for parsing quantum chemistry output files
 */
bool LoadFile::loclabel(std::ifstream& file, const std::string& label, int skip)
{
    std::string    line;
    std::streampos pos;
    while (std::getline(file, line))
    {
        pos = file.tellg() - std::streamoff(line.length() + 1);  // Position at line start
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        if (line.find(label) != std::string::npos)
        {
            // std::cerr << "DEBUG: loclabel matched line: [" << line << "] for label: [" << label << "]" << "\n";
            file.clear();
            file.seekg(pos);  // Rewind to start of matching line
            for (int i = 0; i < skip; ++i)
            {
                if (!std::getline(file, line))
                    return false;
            }
            return true;
        }
    }
    // std::cerr << "Debug: Reached EOF or error while searching for: " << label << "\n";
    return false;
}

bool LoadFile::loclabelfinal(std::ifstream& file, const std::string& label, int& ncount)
{
    ncount = 0;
    std::streampos lastpos;
    std::string    line;
    file.clear();
    file.seekg(0);
    while (std::getline(file, line))
    {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        if (line.find(label) != std::string::npos)
        {
            ncount++;
            lastpos = file.tellg() - std::streamoff(line.length() + 1);  // Save start of line
        }
    }
    if (ncount > 0)
    {
        file.clear();
        file.seekg(lastpos);
        return true;
    }
    return false;
}

//void LoadFile::skiplines(std::ifstream& file, int n)
//{
//    std::string line;
//    for (int i = 0; i < n; ++i)
//    {
//        std::getline(file, line);
//    }
//}

void LoadFile::skiplines(std::ifstream& file, int n, bool print_debug)
{
    std::string line;
    for (int i = 0; i < n; ++i)
    {
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Failed to skip line " + std::to_string(i + 1));
        }
        if (print_debug)
        {
            // Remove trailing \r for Windows compatibility
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            std::cout << "Debug: Skipped line " << i + 1 << ": '" << line << "'" << "\n";
        }
    }
}

// double LoadFile::readaftersign(std::ifstream& file, const std::string& sign) {
//     std::string line;
//     // Get the current position
//     std::streampos current_pos = file.tellg();
//
//     if (!std::getline(file, line)) {
//         throw std::runtime_error("Could not read line for sign: " + sign);
//     }
//
//     // Go back to the beginning of the line we just read
//     file.seekg(current_pos);
//
//     if (line.find(sign) == std::string::npos) {
//         throw std::runtime_error("Sign '" + sign + "' not found in line: " + line);
//     }
//     // Find the last number in the line (handles ".... 4", "....   4", etc.)
//     std::istringstream iss(line);
//     std::string token;
//     double value = 0.0;
//     while (iss >> token) {
//         try {
//             value = std::stod(token); // Keep last number
//         } catch (const std::invalid_argument&) {
//             // Not a number, continue
//         }
//     }
//     if (value == 0.0) {
//         throw std::runtime_error("No number found after sign in line: " + line);
//     }
//     return value;
// }
// double LoadFile::readaftersign(std::ifstream& file, const std::string& sign) {
//     std::string line;
//     // Get the current position
//     std::streampos current_pos = file.tellg();
//
//     // Read the current line (not the next one)
//     if (!std::getline(file, line)) {
//         throw std::runtime_error("Could not read line for sign: " + sign);
//     }
//
//     // Go back to the beginning of the line we just read
//     file.seekg(current_pos);
//
//     // Find the last number in the line
//     std::istringstream iss(line);
//     std::string token;
//     double value = 0.0;
//     bool found_number = false;
//     while (iss >> token) {
//         try {
//             value = std::stod(token); // Keep last number
//             found_number = true;
//         } catch (const std::invalid_argument&) {
//             // Not a number, continue
//         }
//     }
//
//     if (!found_number) {
//         throw std::runtime_error("No number found after sign in line: " + line);
//     }
//     return value;
// }
//
double LoadFile::readaftersign(std::ifstream& file, const std::string& sign)
{
    std::string line;
    if (!std::getline(file, line))
    {
        throw std::runtime_error("Could not read line for sign: " + sign);
    }

    // Find position of sign
    size_t pos = line.find(sign);
    if (pos == std::string::npos)
    {
        // Try reading the next line if sign not found in current line
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Could not read next line for sign: " + sign);
        }
        pos = line.find(sign);
        if (pos == std::string::npos)
        {
            throw std::runtime_error("Sign '" + sign + "' not found in line: " + line);
        }
    }
    // std::cerr << "DEBUG: readaftersign read line: [" << line << "]" << "\n";  // ← ADD THIS

    // Read from after the sign
    std::istringstream iss(line.substr(pos + sign.length()));
    std::string        token;
    double             value = 0.0;
    bool               found = false;
    while (iss >> token)
    {
        try
        {
            value = std::stod(token);
            found = true;
            break;  // Take first number after sign (safer)
        }
        catch (const std::invalid_argument&)
        {
            continue;
        }
    }

    if (!found)
    {
        throw std::runtime_error("No number found after sign '" + sign + "' in line: " + line);
    }

    return value;
}

int LoadFile::readaftersign_int(std::ifstream& file, const std::string& sign)
{
    double value = readaftersign(file, sign);
    return static_cast<int>(value);
}

double LoadFile::readaftersign_from_line(const std::string& line, const std::string& sign)
{
    size_t pos = line.rfind(sign);  // Find last occurrence to match Fortran behavior
    if (pos != std::string::npos)
    {
        std::string value_str = line.substr(pos + sign.length());
        // Remove leading whitespace
        value_str.erase(0, value_str.find_first_not_of(" \t"));
        try
        {
            return std::stod(value_str);
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "ERROR: Invalid argument for stod on value_str: '" << value_str << "' - " << e.what() << '\n';
            throw std::runtime_error("Failed to parse numeric value from: " + value_str);
        }
        catch (const std::out_of_range& e)
        {
            std::cerr << "ERROR: Out of range for stod on value_str: '" << value_str << "' - " << e.what() << '\n';
            throw std::runtime_error("Numeric value out of range: " + value_str);
        }
    }
    throw std::runtime_error("Sign '" + sign + "' not found in line: " + line);
}

void LoadFile::elename2idx(const std::string& element, int& index)
{
    // Convert element name to atomic number using the ind2name array from chemsys.h
    std::string elem = element;
    // Pad single character elements with space for comparison
    if (elem.length() == 1)
        elem += " ";

    for (int i = 1; i <= nelesupp; ++i)
    {
        if (ind2name[i] == elem)
        {
            index = i;
            return;
        }
    }

    // If not found in standard list, try common variations
    if (element == "H")
        index = 1;
    else if (element == "He")
        index = 2;
    else if (element == "Li")
        index = 3;
    else if (element == "Be")
        index = 4;
    else if (element == "B")
        index = 5;
    else if (element == "C")
        index = 6;
    else if (element == "N")
        index = 7;
    else if (element == "O")
        index = 8;
    else if (element == "F")
        index = 9;
    else if (element == "Ne")
        index = 10;
    else if (element == "Na")
        index = 11;
    else if (element == "Mg")
        index = 12;
    else if (element == "Al")
        index = 13;
    else if (element == "Si")
        index = 14;
    else if (element == "P")
        index = 15;
    else if (element == "S")
        index = 16;
    else if (element == "Cl")
        index = 17;
    else if (element == "Ar")
        index = 18;
    else if (element == "K")
        index = 19;
    else if (element == "Ca")
        index = 20;
    else
        index = 0;  // Unknown element
}

void LoadFile::setatmmass(SystemData& sys)
{
    // Set atomic masses based on atomic numbers using elemass array from chemsys.h
    for (auto& atom : sys.a)
    {
        if (atom.index > 0 && atom.index <= nelesupp)
        {
            atom.mass = elemass[atom.index];
        }
        else
        {
            atom.mass = 1.0;  // Default for unknown elements
        }
    }
}

/**
 * @brief Load molecular data from OpenThermo's native .otm format file
 *
 * Parses OpenThermo's binary-like text
 * format containing molecular geometry,
 * vibrational frequencies, and electronic energy data. The .otm format uses
 *
 * simple keyword-based sections for different data types.
 *
 * @param sys SystemData structure to populate with loaded
 * molecular data
 *
 * @throws std::runtime_error if file cannot be opened or parsed
 * @note .otm files contain: *E
 * (energy), *wavenum (frequencies), *atoms (geometry)
 * @note This is OpenThermo's native format for storing processed
 * molecular data
 */
void LoadFile::loadotm(SystemData& sys)
{
    // Load file in binary mode to avoid OS-dependent line endings
    std::ifstream file(sys.inputfile, std::ios::binary);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open input file: " + sys.inputfile);
    }

    std::string line, strtmp;

    // Load energy
    if (loclabel(file, "*E"))
    {
        file >> sys.E;
    }

    // Load wave numbers
    if (loclabel(file, "*wavenum"))
    {
        sys.nfreq = 0;
        double         tmpval;
        std::streampos pos = file.tellg();
        while (file >> tmpval)
        {
            sys.nfreq++;
        }

        sys.wavenum.resize(sys.nfreq);
        sys.freq.resize(sys.nfreq);
        file.clear();
        file.seekg(pos);

        for (int i = 0; i < sys.nfreq; ++i)
        {
            file >> sys.wavenum[i];
            sys.freq[i] = sys.wavenum[i] * wave2freq;  // Convert cm^-1 to Hz
        }
    }

    // Load atoms
    file.clear();
    file.seekg(0);
    if (loclabel(file, "*atoms"))
    {
        sys.ncenter = 0;
        std::string    loadArgs;
        std::streampos pos = file.tellg();

        while (std::getline(file, loadArgs) && !loadArgs.empty() && loadArgs.find('*') == std::string::npos)
        {
            if (!loadArgs.empty() && loadArgs != " ")
                sys.ncenter++;
        }

        sys.a.resize(sys.ncenter);
        file.clear();
        file.seekg(pos);

        for (int i = 0; i < sys.ncenter; ++i)
        {
            file >> strtmp >> sys.a[i].mass >> sys.a[i].x >> sys.a[i].y >> sys.a[i].z;
            elename2idx(strtmp, sys.a[i].index);
        }
    }

    // Load energy levels
    file.clear();
    file.seekg(0);
    if (loclabel(file, "*elevel"))
    {
        sys.nelevel = 0;
        std::string    loadArgs;
        std::streampos pos = file.tellg();

        while (std::getline(file, loadArgs) && !loadArgs.empty() && loadArgs.find('*') == std::string::npos)
        {
            if (!loadArgs.empty() && loadArgs != " ")
                sys.nelevel++;
        }

        sys.elevel.resize(sys.nelevel);
        sys.edegen.resize(sys.nelevel);
        file.clear();
        file.seekg(pos);

        for (int i = 0; i < sys.nelevel; ++i)
        {
            std::getline(file, line);
            std::istringstream iss(line);
            if (!(iss >> sys.elevel[i] >> sys.edegen[i]))
            {
                iss.clear();
                iss.str(line);
                iss >> sys.elevel[i];
                sys.edegen[i] = 1;
            }
        }
    }

    sys.spinmult = 0;
    file.close();
}

void LoadFile::loadgau(SystemData& sys)
{
    // Load file in binary mode to avoid OS-dependent line endings
    std::ifstream file(sys.inputfile, std::ios::binary);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open input file: " + sys.inputfile);
    }

    // Load energy
    if (loclabel(file, "Sum of electronic and zero-point Energies=", 0))
    {
        double tmp1 = readaftersign(file, "=");
        file.clear();
        file.seekg(0);
        if (loclabel(file, "Zero-point correction=", 0))
        {
            double tmp2 = readaftersign(file, "=");
            sys.E       = tmp1 - tmp2;
        }
    }

    // Load multiplicity
    file.clear();
    file.seekg(0);
    if (loclabel(file, "Multiplicity =", 0))
    {
        sys.spinmult = readaftersign_int(file, "=");
    }
    else
    {
        sys.spinmult = 1;
        std::cout << "Note: \"Multiplicity =\" cannot be found, set spin multiplicity to 1" << '\n';
    }

    // Load geometry
    loadGaugeom(file, sys);

    // Set mass
    if (sys.massmod == 1 || sys.massmod == 2)
    {
        setatmmass(sys);
    }
    else if (sys.massmod == 3)
    {
        file.clear();
        file.seekg(0);
        if (loclabel(file, "has atomic number", 0))
        {
            std::string line;
            int         atom_count = 0;
            while (std::getline(file, line) && atom_count < sys.ncenter)
            {
                if (line.find("has atomic number") != std::string::npos)
                {
                    // Split line by spaces, take the last token as mass
                    std::istringstream       iss(line);
                    std::vector<std::string> tokens;
                    std::string              token;
                    while (iss >> token)
                    {
                        tokens.push_back(token);
                    }
                    if (!tokens.empty())
                    {
                        try
                        {
                            double mass            = std::stod(tokens.back());
                            sys.a[atom_count].mass = mass;
                            atom_count++;
                        }
                        catch (const std::invalid_argument&)
                        {
                            // Skip lines that don't have a valid mass at the end
                        }
                    }
                }
            }
            if (atom_count != sys.ncenter)
            {
                std::cerr << "Warning: Found " << atom_count << " mass entries, expected " << sys.ncenter
                          << ". Using default masses.\n";
                setatmmass(sys);
            }
        }
        else
        {
            std::cerr << "Error: Unable to find atomic mass data from quantum chemical output file!" << "\n";
            std::cerr
                << "Gaussian won't print atomic masses out if #T is used in input route. #, #N, #P should be used"
                << "\n"
                << "the parameter modmass can be set to 1 or 2 in settings.in if you do not want to recalculate your "
                   "system"
                << "\n"
                << "Mismatch in thermochemical data may happen due to different type of atomic masses are used \n";
            std::cerr << "Press ENTER button to exit program" << "\n";
            std::cin.get();
            exit(1);
        }
    }

    // Load frequencies
    loadGaufreq(file, sys);
    file.close();
}


void LoadFile::loadGaugeom(std::ifstream& file, SystemData& sys)
{
    std::string locstr;
    int         nskip = 0;

    for (int itime = 1; itime <= 2; ++itime)
    {
        if (itime == 1)
        {
            locstr = "Input orientation:";
        }
        else
        {
            locstr = "Standard orientation:";
        }

        nskip = 0;
        file.clear();
        file.seekg(0);

        // Count how many times the orientation label appears
        // This matches: do while(.true.) with call loclabel(ifileid,locstr,ifound,0)
        while (loclabel(file, locstr, 0))
        {
            nskip++;
            // This matches: read(ifileid,*)
            std::string dummy;
            std::getline(file, dummy);
        }

        if (nskip > 0)
        {
            // Found at least once
            // This matches: call loclabel(ifileid,locstr,ifound) - NOTE: no 0 parameter!
            file.clear();
            file.seekg(0);
            if (!loclabel(file, locstr, 0))
            {
                std::cerr << "Error: Could not relocate geometry section" << '\n';
                exit(1);
            }

            // This matches: call skiplines(ifileid,5)
            skiplines(file, 5);

            // Count atoms - this matches the Fortran counting loop exactly
            sys.ncenter = 0;
            std::string loadArgs;
            while (std::getline(file, loadArgs))
            {
                if (loadArgs.find("----") != std::string::npos)
                    break;
                if (!loadArgs.empty())
                {
                    std::istringstream iss(loadArgs);
                    int                inouse1, index, inouse2;
                    double             x, y, z;
                    if (iss >> inouse1 >> index >> inouse2 >> x >> y >> z)
                    {
                        sys.ncenter++;
                    }
                }
            }

            if (sys.ncenter == 0)
            {
                std::cerr << "Error: No atoms found in geometry section" << '\n';
                exit(1);
            }

            // Allocate memory for atoms
            sys.a.resize(sys.ncenter);

            // Now do the reading phase - this matches the Fortran reading sequence exactly
            // This matches: rewind(ifileid)
            file.clear();
            file.seekg(0);

            // This matches: do iload=1,nskip
            for (int iload = 1; iload <= nskip; ++iload)
            {
                if (!loclabel(file, locstr, 0))
                {
                    std::cerr << "Error: Could not navigate to geometry section " << iload << '\n';
                    exit(1);
                }
                // This matches: read(ifileid,*)
                std::string dummy;
                std::getline(file, dummy);
            }

            // CRITICAL: This matches: call skiplines(ifileid,4) - NOTE: 4, not 5!
            skiplines(file, 4);

            // Read the geometry data - this matches: do iatm=1,ncenter
            for (int iatm = 0; iatm < sys.ncenter; ++iatm)
            {
                int inouse1, inouse2;
                // This matches: read(ifileid,*) inouse,a(iatm)%index,inouse,a(iatm)%x,a(iatm)%y,a(iatm)%z
                if (!(file >> inouse1 >> sys.a[iatm].index >> inouse2 >> sys.a[iatm].x >> sys.a[iatm].y >>
                      sys.a[iatm].z))
                {
                    std::cerr << "Error: Failed to read atom " << (iatm + 1) << " coordinates" << '\n';
                    exit(1);
                }
            }

            // Successfully loaded geometry - exit the itime loop
            break;
        }
        else
        {
            // No geometry found with this orientation
            if (itime == 2)
            {
                std::cerr << "Error: Failed to load geometry from this file!" << '\n';
                std::cerr << "Press ENTER button to exit" << '\n';
                std::cin.get();
                exit(1);
            }
        }
    }
}


void LoadFile::loadGaufreq(std::ifstream& file, SystemData& sys)
{
    file.clear();
    file.seekg(0);

    // First pass: Count the actual number of frequencies
    int frequencyCount = 0;
    while (loclabel(file, "Frequencies -- ", 0))
    {
        std::string line;
        std::getline(file, line);
        size_t freqPos = line.find("Frequencies -- ");
        if (freqPos != std::string::npos)
        {
            std::string        freqData = line.substr(freqPos + 15);  // Skip "Frequencies -- "
            std::istringstream iss(freqData);
            double             temp;
            int                countOnThisLine = 0;

            // Count how many frequency values are on this line (max 3)
            while (iss >> temp && countOnThisLine < 3)
            {
                countOnThisLine++;
            }

            frequencyCount += countOnThisLine;

            // If we read fewer than 3 frequencies, this is the last line
            if (countOnThisLine < 3)
            {
                break;
            }
        }
    }

    sys.nfreq = frequencyCount;

    if (sys.nfreq == 0)
        return;

    // Allocate arrays
    sys.wavenum.resize(sys.nfreq);
    sys.freq.resize(sys.nfreq);

    // Second pass: Read the actual frequency values
    file.clear();
    file.seekg(0);

    int ilackdata = sys.nfreq;
    int inow      = 0;

    while (ilackdata > 0)
    {
        int iread = (ilackdata > 3) ? 3 : ilackdata;

        if (!loclabel(file, "Frequencies -- ", 0))
            break;

        std::string line;
        std::getline(file, line);
        size_t freqPos = line.find("Frequencies -- ");
        if (freqPos != std::string::npos)
        {
            std::string        freqData = line.substr(freqPos + 15);
            std::istringstream iss(freqData);

            if (iread == 1)
            {
                iss >> sys.wavenum[inow];
            }
            else if (iread == 2)
            {
                iss >> sys.wavenum[inow] >> sys.wavenum[inow + 1];
            }
            else if (iread == 3)
            {
                iss >> sys.wavenum[inow] >> sys.wavenum[inow + 1] >> sys.wavenum[inow + 2];
            }
        }

        ilackdata -= iread;
        inow += iread;
    }

    // Convert wavenumbers to frequencies
    for (int i = 0; i < sys.nfreq; ++i)
    {
        sys.freq[i] = sys.wavenum[i] * wave2freq;
    }
}


// CP2K
void LoadFile::loadCP2K(SystemData& sys)
{
    // Load file in binary mode to avoid OS-dependent line endings
    std::ifstream file(sys.inputfile, std::ios::binary);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open input file: " + sys.inputfile);
    }

    // ==================== Load Spin Multiplicity ====================
    file.clear();
    file.seekg(0);
    if (loclabel(file, "DFT| Multiplicity", 0))
    {
        std::string line;
        std::getline(file, line);
        if (line.length() > 19)
        {
            try
            {
                sys.spinmult = std::stoi(line.substr(19));
            }
            catch (...)
            {
                std::cout << "Warning: Failed to parse spin multiplicity, defaulting to 1." << '\n';
                sys.spinmult = 1;
            }
        }
    }
    else
    {
        std::cout << "Note: Unable to find spin multiplicity information; assume to be singlet" << '\n';
        sys.spinmult = 1;
    }

    // ==================== Load Electronic Energy ====================
    file.clear();
    file.seekg(0);
    if (loclabel(file, "Electronic energy (U)", 0))
    {
        sys.E = readaftersign(file, ":") / 2.62549961709828e3;  // Convert to Hartree
    }
    else
    {
        if (sys.Eexter == 0)
        {
            std::cout << "Warning: Unable to find \"Electronic energy (U)\" from the input file, electronic energy is "
                         "thus set to zero. "
                      << "You should directly specify it via \"E\" parameter in settings.ini" << '\n';
            std::cout << "Press ENTER button to continue" << '\n';
            std::cin.get();
        }
        sys.E = 0;
    }

    // ==================== Load Number of Atoms ====================
    file.clear();
    file.seekg(0);
    if (loclabel(file, "- Atoms:", 0))
    {
        sys.ncenter = readaftersign_int(file, ":");
        // std::cerr << "Debug: Found " << sys.ncenter << " atoms." << "\n";
    }
    else
    {
        throw std::runtime_error("'- Atoms:' section not found in CP2K output.");
    }

    sys.a.resize(sys.ncenter);

    // ==================== Load Atom Coordinates and Masses ====================
    file.clear();
    file.seekg(0);
    bool foundAtoms = false;
    if (loclabel(file, "Atom  Kind  Element ", 0) || loclabel(file, "Atom Kind Element ", 0))
    {
        foundAtoms = true;
    }

    if (!foundAtoms)
    {
        std::cerr << "Error: Unable to find atom information! Please make sure that PRINT_LEVEL has been set to MEDIUM "
                     "or higher."
                  << '\n';
        std::cerr << "Press ENTER to exit..." << '\n';
        std::cin.get();
        exit(1);
    }

    // Skip header line
    std::string line;
    std::getline(file, line);

    // Skip any blank lines
    while (std::getline(file, line) && line.find_first_not_of(" \t\r\n") == std::string::npos)
    {
        continue;
    }

    // Parse atoms
    for (int i = 0; i < sys.ncenter; ++i)
    {
        if (line.empty())
        {
            if (!std::getline(file, line))
            {
                throw std::runtime_error("Unexpected end of file while reading atom " + std::to_string(i + 1));
            }
            // Skip empty lines
            while (line.find_first_not_of(" \t\r\n") == std::string::npos)
            {
                if (!std::getline(file, line))
                    break;
            }
        }

        std::istringstream iss(line);
        int                atom_idx, kind_idx, atomic_number;
        std::string        element_symbol;
        double             x, y, z, zeff, mass;

        if (!(iss >> atom_idx >> kind_idx >> element_symbol >> atomic_number >> x >> y >> z >> zeff >> mass))
        {
            throw std::runtime_error("Failed to parse atom line " + std::to_string(i + 1) + ": " + line);
        }

        // Set element index from symbol (optional: you could also use atomic_number)
        elename2idx(element_symbol, sys.a[i].index);
        sys.a[i].x    = x;
        sys.a[i].y    = y;
        sys.a[i].z    = z;
        sys.a[i].mass = mass;  // ← CORRECT: Read 9th field as mass

        // std::cerr << "Debug: Atom " << (i + 1) << " (" << element_symbol << ") index: " << sys.a[i].index
        //           << " mass: " << mass << " amu" << "\n";

        // Prepare next line
        line = "";
        if (i < sys.ncenter - 1)
        {
            while (line.empty() && std::getline(file, line))
            {
                // Skip empty lines
                if (line.find_first_not_of(" \t\r\n") == std::string::npos)
                {
                    line = "";
                }
            }
        }
    }

    // Apply default masses if requested
    if (sys.massmod == 1 || sys.massmod == 2)
    {
        setatmmass(sys);
    }

    // ==================== Load Vibrational Frequencies ====================
    file.clear();
    file.seekg(0);

    std::vector<double> allFreqs;
    std::string         freqLine;

    while (std::getline(file, freqLine))
    {
        if (freqLine.find("VIB|Frequency (cm^-1)") != std::string::npos)
        {
            size_t pos = freqLine.find("VIB|Frequency (cm^-1)");
            if (pos == std::string::npos)
                continue;

            // Extract substring after the label
            std::string        freqPart = freqLine.substr(pos + 22);  // "VIB|Frequency (cm^-1)" is 22 chars
            std::istringstream iss(freqPart);
            double             freq;

            // Read all frequencies on this line
            while (iss >> freq)
            {
                // Optional: Skip near-zero frequencies (translations/rotations)
                // if (std::abs(freq) > 10.0) {
                allFreqs.push_back(freq);
                //}
            }
            // std::cerr << "Debug: Parsed frequency line: " << freqPart << " → found " << allFreqs.size()
            //           << " total so far." << "\n";
        }
    }

    sys.nfreq = allFreqs.size();
    if (sys.nfreq > 0)
    {
        sys.wavenum.resize(sys.nfreq);
        sys.freq.resize(sys.nfreq);

        for (int i = 0; i < sys.nfreq; ++i)
        {
            sys.wavenum[i] = allFreqs[i];
            sys.freq[i]    = allFreqs[i] * wave2freq;  // Convert cm⁻¹ to Hz
        }

        // std::cerr << "Debug: Loaded " << sys.nfreq << " vibrational frequencies." << "\n";
        // for (int i = 0; i < std::min(5, sys.nfreq); ++i)
        //{
        //     std::cerr << "Debug: Frequency " << (i + 1) << ": " << sys.wavenum[i] << " cm⁻¹" << "\n";
        // }
    }
    else
    {
        std::cerr << "No vibrational frequencies found in CP2K output." << '\n';
    }

    // ==================== Set Point Group ====================
    if (sys.PGnameinit == "?")
    {
        if (sys.ipmode == 1)
        {
            sys.PGnameinit = "C1";
            std::cout
                << "Note: When using CP2K to treat periodic systems or solid states (ipmode=1), OpenThermo does not "
                   "automatically detect point group and simply set it to C1."
                << "You might want to use other point groups, and it can be set manually via \"PGname\" in settings.ini"
                << '\n';
        }
    }

    file.close();
}

// ORCA
void LoadFile::loadorca(SystemData& sys)
{
    // Load file in binary mode to avoid OS-dependent line endings
    std::ifstream file(sys.inputfile, std::ios::binary);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open input file: " + sys.inputfile);
    }

    // Load energy
    // Load energy
    int ncount;
    if (!loclabelfinal(file, "FINAL SINGLE POINT ENERGY", ncount))
    {
        std::cerr << "Error: FINAL SINGLE POINT ENERGY not found in ORCA file" << '\n';
        throw std::runtime_error("Energy section not found");
    }
    std::string line;
    if (!std::getline(file, line))
    {
        throw std::runtime_error("Could not read energy line");
    }
    size_t pos = line.find("FINAL SINGLE POINT ENERGY");
    if (pos != std::string::npos)
    {
        std::istringstream iss(line);
        std::string        token;
        while (iss >> token)
            ;
        try
        {
            sys.E = std::stod(token);
        }
        catch (const std::invalid_argument& e)
        {
            throw std::runtime_error("Failed to parse energy value from: " + line);
        }
    }
    else
    {
        throw std::runtime_error("Energy line format not recognized: " + line);
    }

    // Load multiplicity
    file.clear();
    file.seekg(0);  // Reset to start
    if (loclabel(file, "Ideal value S", 0))
    {
        try
        {
            double tmp   = readaftersign(file, "=");
            sys.spinmult = static_cast<int>(2 * tmp + 1);
        }
        catch (const std::runtime_error& e)
        {
            std::cerr << "Warning: Failed to parse spin multiplicity: " << e.what() << '\n';
            sys.spinmult = 1;
        }
    }
    else
    {
        std::cerr << "Note: Ideal value S not found, assuming singlet (spin multiplicity = 1)" << '\n';
        sys.spinmult = 1;
    }

    // Load geometry - CRITICAL FIX: Pass file handle, don't reset
    loadORCAgeom(file, sys);

    // Set mass - CRITICAL FIX: Exact pattern matching and positioning
    if (sys.massmod == 1 || sys.massmod == 2)
    {
        setatmmass(sys);
    }
    else if (sys.massmod == 3)
    {
        file.clear();
        file.seekg(0);

        std::string line;
        bool        found = false;
        while (std::getline(file, line))
        {
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            line.erase(0, line.find_first_not_of(" \t"));
            line.erase(line.find_last_not_of(" \t") + 1);
            if (line == "CARTESIAN COORDINATES (A.U.)")
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            throw std::runtime_error("CARTESIAN COORDINATES (A.U.) section not found for mass reading");
        }

        skiplines(file, 2);  // Skip header lines

        for (int i = 0; i < sys.ncenter; ++i)
        {
            if (!std::getline(file, line))
            {
                throw std::runtime_error("Failed to read data for atom " + std::to_string(i + 1));
            }
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            line.erase(0, line.find_first_not_of(" \t"));
            line.erase(line.find_last_not_of(" \t") + 1);

            std::istringstream iss(line);
            int                no, frag;
            std::string        lb;
            double             za, mass, x, y, z;
            if (!(iss >> no >> lb >> za >> frag >> mass >> x >> y >> z))
            {
                throw std::runtime_error("Failed to parse mass for atom " + std::to_string(i + 1) +
                                         " from line: " + line);
            }
            sys.a[i].mass = mass;
        }

        // Validate masses
        for (int i = 0; i < sys.ncenter; ++i)
        {
            if (sys.a[i].mass < 0.1)
            {
                throw std::runtime_error("Invalid mass for atom " + std::to_string(i + 1) + ": " +
                                         std::to_string(sys.a[i].mass) + " amu");
            }
        }
    }

    // Load frequencies - CRITICAL FIX: Pass file handle, don't reset
    loadORCAfreq(file, sys);
    file.close();
}

void LoadFile::loadORCAgeom(std::ifstream& file, SystemData& sys)
{
    // Reset file to start to ensure "Number of atoms" is found
    file.clear();
    file.seekg(0);

    if (!loclabel(file, "Number of atoms", 0))
    {
        // std::cerr << "Debug: Current file position: " << file.tellg() << "\n";
        throw std::runtime_error("Number of atoms not found in ORCA file");
    }

    try
    {
        sys.ncenter = readaftersign_int(file, ". ");
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "Error: Failed to parse number of atoms: " << e.what() << '\n';
        throw;
    }

    sys.a.resize(sys.ncenter);

    int ncount;
    if (!loclabelfinal(file, "CARTESIAN COORDINATES (ANGSTROEM)", ncount))
    {
        throw std::runtime_error("CARTESIAN COORDINATES (ANGSTROEM) section not found");
    }

    std::string dummy;
    std::getline(file, dummy);  // Skip first header line
    std::getline(file, dummy);  // Skip second header line

    for (int i = 0; i < sys.ncenter; ++i)
    {
        std::string line;
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Failed to read coordinates for atom " + std::to_string(i + 1));
        }
        // Remove trailing \r for Windows compatibility
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        std::istringstream       iss(line);
        std::vector<std::string> tokens;
        std::string              token;
        while (iss >> token)
        {
            tokens.push_back(token);
        }
        if (tokens.size() < 4)
        {
            throw std::runtime_error("Insufficient tokens in coordinate line for atom " + std::to_string(i + 1) + ": " +
                                     line);
        }
        try
        {
            sys.a[i].x = std::stod(tokens[1]);
            sys.a[i].y = std::stod(tokens[2]);
            sys.a[i].z = std::stod(tokens[3]);
            elename2idx(tokens[0], sys.a[i].index);
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error("Failed to parse coordinates for atom " + std::to_string(i + 1) + ": " + e.what() +
                                     " from: " + line);
        }
    }
}

void LoadFile::loadORCAfreq(std::ifstream& file, SystemData& sys)
{
    int ncount;
    if (!loclabelfinal(file, "Scaling factor for frequencies =", ncount) || ncount == 0)
    {
        std::cerr << "\nError: Unable to load frequencies from this file! Please check keywords in the ORCA input file"
                  << '\n';
        std::cerr << "\nOr your ORCA is maybe too old (2.x). Please use the recent versions" << '\n';
        // std::cerr << "\n";
        // std::cerr << "Press ENTER button to exit" << "\n";
        // std::cin.get();
        throw std::runtime_error("Unable to load frequencies from ORCA file: Scaling factor section not found");
        std::cerr << "\n" << '\n';
    }

    std::string dummy;
    std::getline(file, dummy);  // Skip first header line
    std::getline(file, dummy);  // Skip second header line

    sys.nfreq                 = 0;
    std::streampos countStart = file.tellg();
    std::string    loadArgs;
    while (std::getline(file, loadArgs))
    {
        if (loadArgs.find_first_not_of(" \t\r\n") == std::string::npos)
            break;
        if (loadArgs.find(" 0.00 cm") != std::string::npos)
            continue;
        sys.nfreq++;
    }

    sys.wavenum.resize(sys.nfreq);
    sys.freq.resize(sys.nfreq);

    file.clear();
    file.seekg(countStart);

    int ifreq = 0;
    while (ifreq < sys.nfreq && std::getline(file, loadArgs))
    {
        if (loadArgs.find_first_not_of(" \t\r\n") == std::string::npos)
            break;
        if (loadArgs.find(" 0.00 cm") != std::string::npos)
            continue;

        std::istringstream iss(loadArgs);
        std::string        dummy_str;
        double             freq_val;
        if (!(iss >> dummy_str >> freq_val))
        {
            std::cerr << "Error: Failed to parse frequency from line: " << loadArgs << '\n';
            throw std::runtime_error("Invalid frequency format");
        }

        sys.wavenum[ifreq] = freq_val;
        sys.freq[ifreq]    = freq_val * wave2freq;
        ifreq++;
    }
}

// GAMESS
void LoadFile::loadgms(SystemData& sys)
{
    // Load file in binary mode to avoid OS-dependent line endings
    std::ifstream file(sys.inputfile, std::ios::binary);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open input file: " + sys.inputfile);
    }

    // Load energy - robust method: find last "FINAL" line, read next line, split and take last value
    file.clear();
    file.seekg(0);  // Rewind to start
    std::string    line;
    std::streampos last_final_pos = -1;
    while (std::getline(file, line))
    {
        if (line.find("FINAL") != std::string::npos)
        {
            last_final_pos = file.tellg();
        }
    }
    if (last_final_pos == -1)
    {
        throw std::runtime_error("FINAL energy section not found in GAMESS file");
    }
    // Seek to the position after the last "FINAL" line
    file.clear();
    file.seekg(last_final_pos);
    if (!std::getline(file, line))
    {
        throw std::runtime_error("Could not read energy line after FINAL");
    }
    // Split the line by spaces and take the last token as energy
    std::istringstream       iss(line);
    std::vector<std::string> tokens;
    std::string              token;
    while (iss >> token)
    {
        tokens.push_back(token);
    }
    if (tokens.empty())
    {
        throw std::runtime_error("Energy line is empty: " + line);
    }
    // The last token should be the energy value
    std::istringstream energy_iss(tokens.back());
    if (!(energy_iss >> sys.E))
    {
        throw std::runtime_error("Failed to parse energy from last token: " + tokens.back());
    }
    // Try to extract method name from the line
    std::string method = line;
    // Remove trailing spaces and the energy value
    method.erase(method.find_last_not_of(" \t") + 1);
    size_t last_space = method.find_last_of(" \t");
    if (last_space != std::string::npos)
    {
        method = method.substr(0, last_space);
    }
    // Clean up method name
    method.erase(0, method.find_first_not_of(" \t"));
    method.erase(method.find_last_not_of(" \t") + 1);
    if (!method.empty())
    {
        std::cout << "Note: " << method << " energy (" << std::fixed << std::setprecision(8) << sys.E
                  << " a.u.) is loaded. If this is not the intended method, the final U, H, G may be misleading"
                  << '\n';
    }
    else
    {
        std::cout << "Note: Energy (" << std::fixed << std::setprecision(8) << sys.E
                  << " a.u.) is loaded from GAMESS file." << '\n';
    }

    // Load multiplicity
    file.clear();
    file.seekg(0);  // Rewind to start
    if (loclabel(file, "SPIN MULTIPLICITY", 0))
    {
        sys.spinmult = readaftersign_int(file, "=");
    }
    else
    {
        std::cerr << "Warning: SPIN MULTIPLICITY not found, assuming singlet (spin multiplicity = 1)" << '\n';
        sys.spinmult = 1;
    }

    // Load geometry
    loadGmsgeom(file, sys);

    // Set mass
    file.clear();
    file.seekg(0);  // Rewind to start
    if (sys.massmod == 1 || sys.massmod == 2)
    {
        setatmmass(sys);
    }
    else if (sys.massmod == 3)
    {
        file.clear();
        file.seekg(0);  // Rewind to start

        // Use loclabel to find the "ATOMIC WEIGHTS (AMU)" label.
        // The 'skip' parameter in loclabel will be 0, meaning it leaves the
        // file pointer at the *beginning* of the line containing the label.
        if (loclabel(file, "ATOMIC WEIGHTS (AMU)", 0))
        {
            std::string current_line;

            // Consume the "ATOMIC WEIGHTS (AMU)" label line itself.
            // After this, the file pointer is at the beginning of the next line (the blank one).
            if (!std::getline(file, current_line))
            {
                throw std::runtime_error("Failed to read 'ATOMIC WEIGHTS (AMU)' label line.");
            }
            // std::cerr << "Debug (Mass): Consumed label line: '" << current_line << "'" << "\n";

            // Consume the blank line after the label.
            // After this, the file pointer is at the beginning of the first atom's mass data.
            if (!std::getline(file, current_line))
            {
                throw std::runtime_error("Failed to read blank line after 'ATOMIC WEIGHTS (AMU)'.");
            }
            // std::cerr << "Debug (Mass): Consumed blank line: '" << current_line << "'" << "\n";


            std::vector<std::pair<std::string, double>> mass_data;
            // Now, the loop should start reading from the first atom's mass data.
            // We use `std::getline(file, current_line)` directly in the loop condition
            // to fetch the next line for processing.
            while (std::getline(file, current_line))
            {
                // Remove carriage returns and trim leading/trailing whitespace
                current_line.erase(std::remove(current_line.begin(), current_line.end(), '\r'), current_line.end());
                std::string trimmed_line = current_line;
                trimmed_line.erase(0, trimmed_line.find_first_not_of(" \t"));
                trimmed_line.erase(trimmed_line.find_last_not_of(" \t") + 1);

                if (trimmed_line.empty())
                {
                    // std::cerr << "Debug (Mass): Encountered an empty/whitespace-only line, continuing." << "\n";
                    continue;  // Skip truly empty or whitespace-only lines.
                }

                // Check if the line starts with an integer, which indicates an atom entry.
                // This is a robust way to identify the end of the mass data section.
                std::istringstream test_iss(trimmed_line);
                int                atom_num_check;
                if (!(test_iss >> atom_num_check))
                {
                    // std::cerr << "Debug (Mass): Stopping mass reading at non-atom line (e.g., 'MODES...'): '"
                    //           << trimmed_line << "'" << "\n";
                    break;  // If it doesn't start with an integer, it's not atom data, so we've reached the end of the
                            // section.
                }

                // Parse the mass data from the trimmed line.
                // Example line: "    1     C                12.00000"
                std::istringstream iss(trimmed_line);
                int                inouse;          // Atom number (e.g., 1)
                std::string        element_symbol;  // Element symbol (e.g., "C")
                double             mass_value;      // Mass value (e.g., 12.00000)
                if (!(iss >> inouse >> element_symbol >> mass_value))
                {
                    // This should ideally not happen if test_iss succeeded, but good for robust error handling.
                    throw std::runtime_error("Failed to parse mass data from line: '" + current_line + "'");
                }
                mass_data.emplace_back(element_symbol, mass_value);
                // std::cerr << "Debug (Mass): Read mass for element " << element_symbol << ": " << mass_value
                //           << " amu from line: '" << trimmed_line << "'" << "\n";
            }

            // After the loop, verify that the number of masses read matches the expected number of atoms.
            if (mass_data.size() != static_cast<size_t>(sys.ncenter))
            {
                throw std::runtime_error("Mismatch in number of atoms read (" + std::to_string(mass_data.size()) +
                                         ") vs expected (" + std::to_string(sys.ncenter) + ") for mass data.");
            }

            // Assign the read masses to the SystemData structure based on element symbol.
            for (int i = 0; i < sys.ncenter; ++i)
            {
                std::string expected_element_symbol =
                    ind2name[sys.a[i].index];  // Get element symbol from atom index (e.g., 6 -> "C")

                // Trim any trailing spaces from the expected element symbol
                expected_element_symbol.erase(expected_element_symbol.find_last_not_of(" \t") + 1);

                bool found = false;
                for (const auto& [element_symbol_read, mass_read] : mass_data)
                {
                    // Also trim the read element symbol to ensure clean comparison
                    std::string trimmed_read_symbol = element_symbol_read;
                    trimmed_read_symbol.erase(trimmed_read_symbol.find_last_not_of(" \t") + 1);

                    if (trimmed_read_symbol == expected_element_symbol)
                    {
                        sys.a[i].mass = mass_read;
                        found         = true;
                        // std::cerr << "Debug (Mass): Assigned mass " << mass_read << " amu to atom " << i + 1 << " ("
                        //           << expected_element_symbol << ")" << "\n";
                        break;
                    }
                }
                if (!found)
                {
                    // Enhanced error message with both trimmed and original symbols for debugging
                    // std::cerr << "Debug (Mass): Available elements in mass_data:" << "\n";
                    for (const auto& [sym, mass] : mass_data)
                    {
                        std::cerr << "  '" << sym << "'" << '\n';
                    }
                    // std::cerr << "Debug (Mass): Looking for: '" << expected_element_symbol << "' (original: '"
                    //           << ind2name[sys.a[i].index] << "')" << "\n";
                    throw std::runtime_error("No mass found for element '" + expected_element_symbol + "' for atom " +
                                             std::to_string(i + 1));
                }
            }
        }
        else
        {
            throw std::runtime_error("'ATOMIC WEIGHTS (AMU)' section not found in GAMESS file.");
        }
    }

    // Load frequencies
    file.clear();
    file.seekg(0);  // Rewind to start
    loadGmsfreq(file, sys);
    file.close();
}

void LoadFile::loadGmsgeom(std::ifstream& file, SystemData& sys)
{
    file.clear();
    file.seekg(0);  // Rewind to start
    if (loclabel(file, "TOTAL NUMBER OF ATOMS", 0))
    {
        sys.ncenter = readaftersign_int(file, "=");
    }
    else
    {
        throw std::runtime_error("TOTAL NUMBER OF ATOMS not found in GAMESS file");
    }

    sys.a.resize(sys.ncenter);

    int ncount;
    if (loclabelfinal(file, "EQUILIBRIUM GEOMETRY LOCATED", ncount) && ncount > 0)
    {
        //Debug
        std::cout << "Debug: Found 'EQUILIBRIUM GEOMETRY LOCATED' section" << "\n"
                  << "Debug: ncount = " << ncount << "\n";
        skiplines(file, 4, true);

        for (int i = 0; i < sys.ncenter;)
        {
            std::string line;
            if (!std::getline(file, line))
            {
                throw std::runtime_error("Failed to read coordinates for atom " + std::to_string(i + 1));
            }
            // Remove trailing \r for Windows compatibility
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
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
            if (line.empty())
                continue;  // Skip blank lines
            // std::cerr << "Debug: Geometry line for atom " << i + 1 << ": " << line << "\n";
            std::istringstream       iss(line);
            std::vector<std::string> tokens;
            std::string              token;
            while (iss >> token)
            {
                tokens.push_back(token);
            }
            if (tokens.size() < 5)
            {
                throw std::runtime_error("Insufficient tokens in coordinate line for atom (ANGS)" + std::to_string(i + 1) +
                                         ": " + line);
            }
            try
            {
                sys.a[i].x     = std::stod(tokens[2]);
                sys.a[i].y     = std::stod(tokens[3]);
                sys.a[i].z     = std::stod(tokens[4]);
                sys.a[i].index = static_cast<int>(std::stod(tokens[1]));  // Atomic number from charge
            }
            catch (const std::exception& e)
            {
                throw std::runtime_error("Failed to parse coordinates for atom " + std::to_string(i + 1) + ": " +
                                         e.what() + " from: " + line);
            }
            // std::cerr << "Debug: Atom " << i + 1 << " index = " << sys.a[i].index << "\n";
            ++i;
        }
    }
    else
    {
        file.clear();
        file.seekg(0);
        if (loclabel(file, "COORDINATES (BOHR)"))
        {
            skiplines(file, 2);
            for (int i = 0; i < sys.ncenter;)
            {
                std::string line;
                if (!std::getline(file, line))
                {
                    throw std::runtime_error("Failed to read coordinates for atom " + std::to_string(i + 1));
                }
                // Remove trailing \r for Windows compatibility
                line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
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
                if (line.empty())
                    continue;  // Skip blank lines
                std::istringstream       iss(line);
                std::vector<std::string> tokens;
                std::string              token;
                while (iss >> token)
                {
                    tokens.push_back(token);
                }
                if (tokens.size() < 5)
                {
                    throw std::runtime_error("Insufficient tokens in coordinate line for atom (BOHR)" +
                                             std::to_string(i + 1) + ": " + line);
                }
                try
                {
                    sys.a[i].x     = std::stod(tokens[2]);
                    sys.a[i].y     = std::stod(tokens[3]);
                    sys.a[i].z     = std::stod(tokens[4]);
                    sys.a[i].index = static_cast<int>(std::stod(tokens[1]));  // Atomic number from charge
                }
                catch (const std::exception& e)
                {
                    throw std::runtime_error("Failed to parse coordinates for atom " + std::to_string(i + 1) + ": " +
                                             e.what() + " from: " + line);
                }
                sys.a[i].x *= b2a;
                sys.a[i].y *= b2a;
                sys.a[i].z *= b2a;
                ++i;
            }
        }
        else
        {
            throw std::runtime_error("No valid coordinate section found in GAMESS file");
        }
    }
}

void LoadFile::loadGmsfreq(std::ifstream& file, SystemData& sys)
{
    file.clear();
    file.seekg(0);  // Rewind to start
    if (loclabel(file, "VIBRATIONAL MODES ARE USED IN THERMOCHEMISTRY", 0))
    {
        std::string line;
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Failed to read frequency range line");
        }
        std::istringstream iss(line);
        int                istart, iend;
        std::string        dummy;
        if (!(iss >> istart >> dummy >> iend))
        {
            throw std::runtime_error("Failed to parse frequency range: " + line);
        }
        sys.nfreq = iend - istart + 1;

        sys.wavenum.resize(sys.nfreq);
        sys.freq.resize(sys.nfreq);

        file.clear();
        file.seekg(0);
        if (loclabel(file, "MODE FREQ(CM**-1)", 0))
        {
            skiplines(file, 1);  // Skip header line
            int ifreq = 0;
            for (int idx = 1; idx <= iend; ++idx)
            {
                std::string line;
                if (!std::getline(file, line))
                {
                    throw std::runtime_error("Failed to read frequency for mode " + std::to_string(idx));
                }
                std::istringstream iss(line);
                int                inouse;
                double             tmpval;
                if (!(iss >> inouse >> tmpval))
                {
                    throw std::runtime_error("Failed to parse frequency for mode " + std::to_string(idx) +
                                             " from: " + line);
                }
                if (idx >= istart)
                {
                    sys.wavenum[ifreq] = tmpval;
                    sys.freq[ifreq]    = tmpval * wave2freq;
                    ifreq++;
                }
            }
        }
        else
        {
            throw std::runtime_error("MODE FREQ(CM**-1) section not found");
        }
    }
    else
    {
        throw std::runtime_error("VIBRATIONAL MODES ARE USED IN THERMOCHEMISTRY section not found");
    }
}

// NWCHEM
void LoadFile::loadnw(SystemData& sys)
{
    // Load file in binary mode to avoid OS-dependent line endings
    std::ifstream file(sys.inputfile, std::ios::binary);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open input file: " + sys.inputfile);
    }

    // Load energy
    int ncount;
    if (loclabelfinal(file, "Total DFT energy =", ncount) && ncount > 0)
    {
        sys.E = readaftersign(file, "=");
    }
    else
    {
        file.clear();
        file.seekg(0);
        if (loclabelfinal(file, "Total SCF energy =", ncount) && ncount > 0)
        {
            sys.E = readaftersign(file, "=");
            std::cout << "Note: SCF energy is loaded. If the theoretical method presently used is other one, "
                      << "the finally printed total U, H, G will be misleading" << '\n';
        }
        else
        {
            std::cout << "Warning: Unable to load electronic energy, thus it is set to zero" << '\n';
            sys.E = 0;
        }
    }

    // Load multiplicity
    file.clear();
    file.seekg(0);
    if (loclabel(file, "Spin multiplicity:", 0))
    {
        sys.spinmult = readaftersign_int(file, ":");
    }

    // Load geometry
    loadNwgeom(file, sys);

    // Set mass
    if (sys.massmod == 1 || sys.massmod == 2)
    {
        setatmmass(sys);
    }
    else if (sys.massmod == 3)
    {
        file.clear();
        file.seekg(0);
        if (loclabel(file, "- Atom information -"))
        {
            skiplines(file, 3);
            // std::cerr << "Debug: Found 'Atom information' section, reading " << sys.ncenter << " atoms" << "\n";
            int         atoms_read = 0;
            std::string line;
            bool        parse_failed = false;
            while (std::getline(file, line) && atoms_read < sys.ncenter)
            {
                if (line.find("------") != std::string::npos)
                    break;
                if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos)
                    continue;

                std::replace(line.begin(), line.end(), 'D', 'E');  // ← FIX: Handle D notation

                std::istringstream iss(line);
                std::string        element_symbol;
                int                atom_number;
                double             x, y, z, mass;

                if (!(iss >> element_symbol >> atom_number >> x >> y >> z >> mass))
                {
                    std::cerr << "Warning: Failed to parse mass line " << (atoms_read + 1) << ": " << line << '\n';
                    parse_failed = true;
                    break;
                }

                sys.a[atoms_read].mass = mass;
                // std::cerr << "Debug: Atom " << atoms_read + 1 << " (" << element_symbol << ") mass: " << mass << "
                // amu"
                //           << "\n";
                atoms_read++;
            }

            if (atoms_read != sys.ncenter || parse_failed)
            {
                std::cerr << "Warning: Mass section incomplete or corrupt. Falling back to default atomic masses."
                          << '\n';
                setatmmass(sys);  // ← SAFETY: Don’t lose your atoms!
            }
        }
        else
        {
            throw std::runtime_error("'- Atom information -' section not found for mass");
        }
    }

    // Load frequencies
    loadNwfreq(file, sys);
    file.close();
}

void LoadFile::loadNwgeom(std::ifstream& file, SystemData& sys)
{
    file.clear();
    file.seekg(0);
    if (loclabel(file, "No. of atoms     :"))
    {
        // std::cerr << "DEBUG: About to read after sign..." << "\n";
        sys.ncenter = readaftersign_int(file, ":");
        // std::cerr << "Debug: ncenter set to " << sys.ncenter << "\n";
    }
    else
    {
        throw std::runtime_error("'No. of atoms     :' not found in NWChem file");
    }

    sys.a.resize(sys.ncenter);

    file.clear();
    file.seekg(0);
    if (loclabel(file, "- Atom information -"))
    {
        // Skip exactly three header lines
        std::string line;
        // 1. Header line: "---------------------------- Atom information ----------------------------"
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Failed to read header line after '- Atom information -'");
        }
        // std::cerr << "Debug: Skipped header line: " << line << "\n";
        //  2. Column labels: "atom    #        X              Y              Z            mass"
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Failed to read column labels line");
        }
        // std::cerr << "Debug: Skipped column labels: " << line << "\n";
        //  3. Separator: "--------------------------------------------------------------------------"
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Failed to read separator line");
        }
        // std::cerr << "Debug: Skipped separator: " << line << "\n";

        // Now read geometry data
        for (int i = 0; i < sys.ncenter; ++i)
        {
            if (!std::getline(file, line))
            {
                throw std::runtime_error("Failed to read geometry for atom " + std::to_string(i + 1));
            }
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            line.erase(0, line.find_first_not_of(" \t"));
            line.erase(line.find_last_not_of(" \t") + 1);
            if (line.empty())
            {
                // std::cerr << "Debug: Skipping empty line before atom " << i + 1 << "\n";
                --i;  // Retry this atom
                continue;
            }
            if (line.find("------") != std::string::npos)
            {
                throw std::runtime_error("Reached end of Atom information section prematurely at atom " +
                                         std::to_string(i + 1));
            }

            // CRITICAL: Convert Fortran D-exponent to C++ E-exponent
            std::replace(line.begin(), line.end(), 'D', 'E');
            std::replace(line.begin(), line.end(), 'd', 'e');  // or 'E' if you prefer consistency

            // std::cerr << "Debug: Geometry line " << i + 1 << ": " << line << "\n";
            std::istringstream iss(line);
            std::string        strtmp;
            int                inouse;
            double             x, y, z;
            std::string        mass_str;  // Skip mass
            if (!(iss >> strtmp >> inouse >> x >> y >> z >> mass_str))
            {
                throw std::runtime_error("Failed to parse geometry for atom " + std::to_string(i + 1) +
                                         " from: " + line);
            }
            elename2idx(strtmp, sys.a[i].index);
            sys.a[i].x = x * b2a;
            sys.a[i].y = y * b2a;
            sys.a[i].z = z * b2a;
            // std::cerr << "Debug: Atom " << i + 1 << " (" << strtmp << ") index: " << sys.a[i].index << "\n";
        }
    }
    else
    {
        throw std::runtime_error("'- Atom information -' section not found");
    }
}

void LoadFile::loadNwfreq(std::ifstream& file, SystemData& sys)
{
    if (loclabel(file, "Projected Derivative Dipole Moments"))
    {
        skiplines(file, 3);  // Skip: header, units, separator line

        sys.nfreq = 0;
        std::vector<double> tempFreq;

        std::string line;
        while (std::getline(file, line))
        {
            line.erase(0, line.find_first_not_of(" \t"));
            line.erase(line.find_last_not_of(" \t\n\r") + 1);
            // Stop at separator line
            if (line.find("------") != std::string::npos && line.find("Mode") == std::string::npos)
            {
                break;
            }

            // Skip if line doesn't start with a digit
            if (line.empty() || !std::isdigit(line[0]))
            {
                continue;
            }
            // Skip empty lines
            if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos)
            {
                continue;
            }

            // We expect lines like:
            //    7       78.253 ||      -0.000               0.000             2.508

            std::istringstream iss(line);
            int                mode;
            double             freq_cm;
            std::string        sep;

            // Read mode number and frequency
            if (!(iss >> mode >> freq_cm))
            {
                continue;  // Skip malformed lines
            }

            // Optional: read "||" to ensure we're parsing correctly
            iss >> sep;
            if (sep != "||")
            {
                // Not a problem — maybe extra space, but we already have freq
            }

            // Only store non-zero frequencies: (original logic)
            // But note: NWChem prints 0.000 for translations/rotations — you may want to SKIP them
            // Since you said "if (tmp != 0)" — skipping zeros.
            // But in vibrational analysis, we usually want ALL real modes (positive frequencies)
            // Negative frequencies = imaginary = transition states — you might want to keep them too.

            // Decide: Still want to skip zero frequencies?
            // In the example, modes 1-6 are ~0 (translations/rotations) — usually excluded in thermochemistry.
            // Let's keep only POSITIVE frequencies (normal vibrations)

            if (freq_cm > 1e-3)
            {  // Skip near-zero (translations/rotations)
                tempFreq.push_back(freq_cm);
                sys.nfreq++;
            }
            // If want to include negative (imaginary) frequencies, use: if (std::abs(freq_cm) > 1e-3)
        }

        // Resize and convert
        sys.wavenum.resize(sys.nfreq);
        sys.freq.resize(sys.nfreq);

        for (int i = 0; i < sys.nfreq; ++i)
        {
            sys.wavenum[i] = tempFreq[i];
            sys.freq[i]    = tempFreq[i] * wave2freq;  // Convert cm⁻¹ to Hz
        }

        // std::cerr << "Debug: Loaded " << sys.nfreq << " vibrational frequencies." << "\n";
    }
    else
    {
        std::cerr << "'Projected Derivative Dipole Moments' section not found — no frequencies loaded." << "\n";
    }
}

// XTB
void LoadFile::loadxtb(SystemData& sys)
{
    // Load file in binary mode to avoid OS-dependent line endings
    std::ifstream file(sys.inputfile, std::ios::binary);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open input file: " + sys.inputfile);
    }

    // Load multiplicity
    if (loclabel(file, "alpha electrons", 0))
    {
        int         naelec, nbelec;
        std::string dummy;
        file >> naelec >> dummy >> dummy >> nbelec;
        sys.spinmult = naelec - nbelec + 1;
    }

    // Load energy
    if (sys.Eexter == 0)
    {
        std::cout << "\nNOTE: This file does not contain electronic energy, and you also did not explicitly set \"E\" "
                     "in settings.ini. "
                  << "If you want to let OpenThermo load electronic energy from a xtb output file, "
                  << "input its path now, e.g. D:\\ltwd\\xtb.out. If you press ENTER button directly, then electronic "
                     "energy will simply be set to 0"
                  << '\n';

        while (true)
        {
            std::string c200tmp;
            std::getline(std::cin, c200tmp);

            if (c200tmp.empty())
            {
                sys.E = 0;
                break;
            }
            else
            {
                std::ifstream xtbfile(c200tmp);
                if (xtbfile.is_open())
                {
                    int ncount;
                    if (loclabelfinal(xtbfile, "total energy", ncount) && ncount > 0)
                    {
                        std::string line;
                        std::getline(xtbfile, line);
                        if (line.length() > 36)
                        {
                            sys.E = std::stod(line.substr(36));
                            std::cout << "Loaded electronic energy: " << std::fixed << std::setprecision(10) << sys.E
                                      << " Hartree" << '\n';
                        }
                        xtbfile.close();
                        break;
                    }
                    else
                    {
                        std::cout << "Error: Unable to locate \"total energy\"! Input again" << '\n';
                        xtbfile.close();
                        continue;
                    }
                }
                else
                {
                    std::cout << "Cannot find the file, input again!" << '\n';
                }
            }
        }
    }
    else
    {
        sys.E = sys.Eexter;
    }

    // Load geometry
    if (loclabel(file, "Coordinates (Angstroms)", 0))
    {
        skiplines(file, 3);

        sys.ncenter = 0;
        std::string              line;
        std::vector<std::string> atomLines;

        while (std::getline(file, line))
        {
            if (line.find("--") != std::string::npos)
                break;
            if (!line.empty())
            {
                atomLines.push_back(line);
                sys.ncenter++;
            }
        }

        sys.a.resize(sys.ncenter);

        for (int i = 0; i < sys.ncenter; ++i)
        {
            std::istringstream iss(atomLines[i]);
            int                inouse;
            iss >> inouse >> sys.a[i].index >> inouse >> sys.a[i].x >> sys.a[i].y >> sys.a[i].z;
        }
    }

    if (sys.massmod == 3)
    {
        std::cout << "Note: massmod=3 is meaningless for present case because input file does not record atomic mass. "
                     "Now set atomic masses to element mass (massmod=1)"
                  << '\n';
        sys.massmod = 1;
    }
    setatmmass(sys);

    loadGaufreq(file, sys);
    file.close();
}

// VASP
void LoadFile::loadvasp(SystemData& sys)
{
    // Load file in binary mode to avoid OS-dependent line endings
    std::ifstream file(sys.inputfile, std::ios::binary);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open input file: " + sys.inputfile);
    }

    // Determine if it's CONTCAR or OUTCAR based on content
    // Check for OUTCAR-specific patterns first
    bool           isOUTCAR = false;
    std::string    checkLine;
    std::streampos originalPos = file.tellg();
    while (std::getline(file, checkLine))
    {
        if (checkLine.find("VRHFIN") != std::string::npos || checkLine.find("POMASS") != std::string::npos ||
            checkLine.find("FREE ENERGIE OF THE ION-ELECTRON SYSTEM") != std::string::npos)
        {
            isOUTCAR = true;
            break;
        }
    }
    file.clear();
    file.seekg(originalPos);

    if (isOUTCAR)
    {
        // OUTCAR - load geometry, energy, and frequencies
        loadVASPgeom(file, sys, true);
        file.clear();
        file.seekg(0);
        loadVASPEnergy(file, sys);
        file.clear();
        file.seekg(0);
        loadVASPfreq(file, sys);
    }
    else
    {
        // Assume CONTCAR - load geometry
        loadVASPgeom(file, sys, false);
        // For CONTCAR, try to load energy and frequencies from OUTCAR in the same directory
        std::filesystem::path contcarPath = sys.inputfile;
        std::filesystem::path outcarPath  = contcarPath.parent_path() / "OUTCAR";
        if (std::filesystem::exists(outcarPath))
        {
            std::ifstream outcarFile(outcarPath);
            if (outcarFile.is_open())
            {
                loadVASPEnergy(outcarFile, sys);
                outcarFile.clear();
                outcarFile.seekg(0);
                loadVASPfreq(outcarFile, sys);
                outcarFile.close();
            }
        }
        else
        {
            // For CONTCAR, set default atomic masses if no OUTCAR found
            setatmmass(sys);
        }
    }


    // Set default multiplicity for VASP (usually singlet)
    if (sys.spinmult == 0)
    {
        sys.spinmult = 1;
    }

    // ==================== Set Point Group ====================
    if (sys.PGnameinit == "?")
    {
        if (sys.ipmode == 1)
        {
            sys.PGnameinit = "C1";
            std::cout
                << "Note: When using VASP to treat periodic systems or solid states (ipmode=1), OpenThermo does not "
                   "automatically detect point group and simply set it to C1."
                << "You might want to use other point group, and it can be set manually via \"PGname\" in settings.ini"
                << '\n';
        }
    }

    file.close();
}

void LoadFile::loadVASPgeom(std::ifstream& file, SystemData& sys, bool isOUTCAR)
{
    if (isOUTCAR)
    {
        // OUTCAR case: read geometry from OUTCAR, elements from CONTCAR, masses from OUTCAR
        std::filesystem::path outcarPath  = sys.inputfile;
        std::filesystem::path contcarPath = outcarPath.parent_path() / "CONTCAR";

        if (!std::filesystem::exists(contcarPath))
        {
            throw std::runtime_error("CONTCAR file not found in the same directory as OUTCAR: " + contcarPath.string());
        }

        // Read CONTCAR for elements and atom counts
        std::ifstream contcarFile(contcarPath);
        if (!contcarFile.is_open())
        {
            throw std::runtime_error("Cannot open CONTCAR file: " + contcarPath.string());
        }

        // Skip first 5 lines
        for (int i = 0; i < 5; ++i)
        {
            std::string dummy;
            std::getline(contcarFile, dummy);
        }

        // Read element names (line 6)
        std::string elemLine;
        std::getline(contcarFile, elemLine);
        std::istringstream       elemIss(elemLine);
        std::vector<std::string> elements;
        std::string              elem;
        while (elemIss >> elem)
        {
            // Remove trailing slash if present
            if (!elem.empty() && elem.back() == '/')
                elem.pop_back();
            elements.push_back(elem);
        }

        // Read atom counts (line 7)
        std::string countLine;
        std::getline(contcarFile, countLine);
        std::istringstream countIss(countLine);
        std::vector<int>   counts;
        int                count;
        sys.ncenter = 0;
        while (countIss >> count)
        {
            counts.push_back(count);
            sys.ncenter += count;
        }

        if (elements.size() != counts.size())
        {
            throw std::runtime_error("Mismatch between number of elements and counts in CONTCAR");
        }

        contcarFile.close();

        // Read masses from OUTCAR
        std::map<std::string, double> massMap;
        file.clear();
        file.seekg(0);
        std::string line;
        while (std::getline(file, line))
        {
            if (line.find("POMASS") != std::string::npos)
            {
                size_t pos = line.find('=');
                if (pos != std::string::npos)
                {
                    std::string massStr = line.substr(pos + 1);
                    // Find the mass value before ';'
                    size_t semiPos = massStr.find(';');
                    if (semiPos != std::string::npos)
                    {
                        massStr = massStr.substr(0, semiPos);
                    }
                    try
                    {
                        double mass = std::stod(massStr);
                        // Find the element name before POMASS
                        size_t      pomassPos    = line.find("POMASS");
                        std::string beforePomass = line.substr(0, pomassPos);
                        // Find the element symbol (usually 2 characters before POMASS)
                        std::istringstream iss(beforePomass);
                        std::string        token;
                        std::string        elemSymbol;
                        while (iss >> token)
                        {
                            elemSymbol = token;
                        }
                        // Remove trailing colon if present
                        if (!elemSymbol.empty() && elemSymbol.back() == ':')
                            elemSymbol.pop_back();
                        massMap[elemSymbol] = mass;
                    }
                    catch (const std::invalid_argument&)
                    {
                        // Skip invalid mass
                    }
                }
            }
        }

        // Read geometry from OUTCAR
        file.clear();
        file.seekg(0);
        int ncount;
        if (!loclabelfinal(file, "POSITION", ncount) || ncount == 0)
        {
            throw std::runtime_error("POSITION section not found in OUTCAR");
        }

        // Skip header lines
        skiplines(file, 2);  // Skip "POSITION" line and "TOTAL-FORCE" line

        std::vector<std::array<double, 3>> coordinates;
        while (std::getline(file, line))
        {
            if (line.find("----") != std::string::npos)
                break;
            if (line.empty())
                continue;

            std::istringstream iss(line);
            double             x, y, z;
            if (iss >> x >> y >> z)
            {
                coordinates.push_back({x, y, z});
            }
        }

        if (coordinates.size() != static_cast<size_t>(sys.ncenter))
        {
            throw std::runtime_error("Number of coordinates (" + std::to_string(coordinates.size()) +
                                     ") does not match number of atoms (" + std::to_string(sys.ncenter) + ")");
        }

        // Assign elements and masses
        sys.a.resize(sys.ncenter);
        int atomIdx = 0;
        for (size_t i = 0; i < elements.size(); ++i)
        {
            for (int j = 0; j < counts[i]; ++j)
            {
                sys.a[atomIdx].x = coordinates[atomIdx][0];
                sys.a[atomIdx].y = coordinates[atomIdx][1];
                sys.a[atomIdx].z = coordinates[atomIdx][2];

                elename2idx(elements[i], sys.a[atomIdx].index);

                auto massIt = massMap.find(elements[i]);
                if (massIt != massMap.end())
                {
                    sys.a[atomIdx].mass = massIt->second;
                }
                else
                {
                    // Use default mass if not found
                    if (sys.a[atomIdx].index > 0 && sys.a[atomIdx].index <= nelesupp)
                    {
                        sys.a[atomIdx].mass = elemass[sys.a[atomIdx].index];
                    }
                    else
                    {
                        sys.a[atomIdx].mass = 1.0;
                    }
                }

                atomIdx++;
            }
        }
    }
    else
    {
        // CONTCAR case: read geometry from CONTCAR, use default masses
        file.clear();
        file.seekg(0);

        // Skip comment
        std::string dummy;
        std::getline(file, dummy);

        // Skip scaling
        std::getline(file, dummy);

        // Read lattice vectors
        std::array<std::array<double, 3>, 3> lattice;
        for (auto& row : lattice)
        {
            std::getline(file, dummy);
            std::istringstream iss(dummy);
            iss >> row[0] >> row[1] >> row[2];
        }

        // Read element names
        std::string elemLine;
        std::getline(file, elemLine);
        std::istringstream       elemIss(elemLine);
        std::vector<std::string> elements;
        std::string              elem;
        while (elemIss >> elem)
        {
            if (!elem.empty() && elem.back() == '/')
                elem.pop_back();
            elements.push_back(elem);
        }

        // Read atom counts
        std::string countLine;
        std::getline(file, countLine);
        std::istringstream countIss(countLine);
        std::vector<int>   counts;
        int                count;
        sys.ncenter = 0;
        while (countIss >> count)
        {
            counts.push_back(count);
            sys.ncenter += count;
        }

        if (elements.size() != counts.size())
        {
            throw std::runtime_error("Mismatch between number of elements and counts in CONTCAR");
        }

        // Read coordinate type
        std::string coordType;
        std::getline(file, coordType);

        sys.a.resize(sys.ncenter);
        int atomIdx = 0;
        for (size_t i = 0; i < elements.size(); ++i)
        {
            for (int j = 0; j < counts[i]; ++j)
            {
                std::string coordLine;
                std::getline(file, coordLine);
                std::istringstream iss(coordLine);
                double             fx, fy, fz;
                iss >> fx >> fy >> fz;

                elename2idx(elements[i], sys.a[atomIdx].index);

                // Convert Direct to Cartesian
                sys.a[atomIdx].x = fx * lattice[0][0] + fy * lattice[1][0] + fz * lattice[2][0];
                sys.a[atomIdx].y = fx * lattice[0][1] + fy * lattice[1][1] + fz * lattice[2][1];
                sys.a[atomIdx].z = fx * lattice[0][2] + fy * lattice[1][2] + fz * lattice[2][2];

                atomIdx++;
            }
        }

        // Use default masses
        setatmmass(sys);
    }
}

void LoadFile::loadVASPEnergy(std::ifstream& file, SystemData& sys)
{
    int ncount;
    if (loclabelfinal(file, "energy  without entropy", ncount) && ncount > 0)
    {
        // The file pointer is at the beginning of the "energy  without entropy" line
        std::string energy_line;
        std::getline(file, energy_line);

        // Read the whole line and split it to get the value at the fourth position
        // Line format: "  energy  without entropy=      -27.39346935  energy(sigma->0) =      -27.39346935"
        std::istringstream       iss(energy_line);
        std::string              token;
        std::vector<std::string> tokens;

        while (iss >> token)
        {
            tokens.push_back(token);
        }

        if (tokens.size() >= 4)
        {
            try
            {
                size_t energy_index = (sys.vasp_energy_select == 1) ? tokens.size() - 1
                                                                    : 3;  // Last token if extrape=true, else 4th token
                double energy_eV    = std::stod(tokens[energy_index]);
                sys.E               = energy_eV / 27.2114;  // Convert eV to Hartree
            }
            catch (const std::invalid_argument& e)
            {
                std::cerr << "Warning: Invalid energy value '" << tokens[3] << "'. Setting energy to 0." << "\n";
                sys.E = 0.0;
            }
            catch (const std::out_of_range& e)
            {
                std::cerr << "Warning: Energy value out of range. Setting energy to 0." << "\n";
                sys.E = 0.0;
            }
        }
        else
        {
            std::cerr << "Warning: Unexpected format in energy line. Expected at least 4 tokens, got " << tokens.size()
                      << ". Setting energy to 0." << "\n";
            sys.E = 0.0;
        }
    }
    else
    {
        std::cerr << "Warning: Unable to find energy section in VASP file, setting to 0" << '\n';
        sys.E = 0.0;
    }
}

void LoadFile::loadVASPfreq(std::ifstream& file, SystemData& sys)
{
    std::vector<double> allFreqs;
    std::string         freqLine;

    while (std::getline(file, freqLine))
    {
        if (freqLine.find("cm-1") != std::string::npos)
        {
            // Split line into tokens
            std::vector<std::string> tokens;
            std::istringstream       iss(freqLine);
            std::string              token;
            while (iss >> token)
            {
                tokens.push_back(token);
            }

            // Find "cm-1" token and extract frequency from previous token
            for (size_t i = 0; i < tokens.size(); ++i)
            {
                if (tokens[i] == "cm-1" && i > 0)
                {
                    try
                    {
                        double freq = std::stod(tokens[i - 1]);
                        if (freq > 1e-3)  // Skip near-zero frequencies
                            allFreqs.push_back(freq);
                    }
                    catch (const std::invalid_argument&)
                    {
                        // Skip invalid numbers
                    }
                    break;  // Only one frequency per line
                }
            }
        }
    }

    sys.nfreq = allFreqs.size();
    if (sys.nfreq > 0)
    {
        sys.wavenum.resize(sys.nfreq);
        sys.freq.resize(sys.nfreq);

        for (int i = 0; i < sys.nfreq; ++i)
        {
            sys.wavenum[i] = allFreqs[i];
            sys.freq[i]    = allFreqs[i] * wave2freq;  // Convert cm⁻¹ to Hz
        }
    }
    else
    {
        std::cerr << "No vibrational frequencies found in VASP file." << '\n';
    }
}
