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
#include "defvar.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
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
            // std::cerr << "DEBUG: loclabel matched line: [" << line << "] for label: [" << label << "]" << std::endl;
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
    // std::cerr << "Debug: Reached EOF or error while searching for: " << label << std::endl;
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

void LoadFile::skiplines(std::ifstream& file, int n)
{
    std::string line;
    for (int i = 0; i < n; ++i)
    {
        std::getline(file, line);
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
        throw std::runtime_error("Sign '" + sign + "' not found in line: " + line);
    }
    // std::cerr << "DEBUG: readaftersign read line: [" << line << "]" << std::endl;  // ← ADD THIS

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
            std::cerr << "ERROR: Invalid argument for stod on value_str: '" << value_str << "' - " << e.what()
                      << std::endl;
            throw std::runtime_error("Failed to parse numeric value from: " + value_str);
        }
        catch (const std::out_of_range& e)
        {
            std::cerr << "ERROR: Out of range for stod on value_str: '" << value_str << "' - " << e.what() << std::endl;
            throw std::runtime_error("Numeric value out of range: " + value_str);
        }
    }
    throw std::runtime_error("Sign '" + sign + "' not found in line: " + line);
}

void LoadFile::elename2idx(const std::string& element, int& index)
{
    // Convert element name to atomic number using the ind2name array from defvar.h
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
    // Set atomic masses based on atomic numbers using elemass array from defvar.h
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
    std::ifstream file(sys.inputfile);
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
    std::ifstream file(sys.inputfile);
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
        std::cout << "Note: \"Multiplicity =\" cannot be found, set spin multiplicity to 1" << std::endl;
    }

    // Load geometry
    loadGaugeom(file, sys);

    // Set mass
    if (sys.defmass == 1 || sys.defmass == 2)
    {
        setatmmass(sys);
    }
    else if (sys.defmass == 3)
    {
        file.clear();
        file.seekg(0);
        if (loclabel(file, "has atomic number", 0))
        {
            for (int i = 0; i < sys.ncenter; ++i)
            {
                sys.a[i].mass = readaftersign(file, "mass");
            }
        }
        else
        {
            std::cerr << "Error: Unable to find atomic mass information!" << std::endl;
            std::cerr << "If you have used #T for your freq task, you should use # instead" << std::endl;
            std::cerr << "Press ENTER button to exit program" << std::endl;
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
                std::cerr << "Error: Could not relocate geometry section" << std::endl;
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
                std::cerr << "Error: No atoms found in geometry section" << std::endl;
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
                    std::cerr << "Error: Could not navigate to geometry section " << iload << std::endl;
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
                    std::cerr << "Error: Failed to read atom " << (iatm + 1) << " coordinates" << std::endl;
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
                std::cerr << "Error: Failed to load geometry from this file!" << std::endl;
                std::cerr << "Press ENTER button to exit" << std::endl;
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
    std::ifstream file(sys.inputfile);
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
                std::cout << "Warning: Failed to parse spin multiplicity, defaulting to 1." << std::endl;
                sys.spinmult = 1;
            }
        }
    }
    else
    {
        std::cout << "Note: Unable to find spin multiplicity information, assume to be singlet" << std::endl;
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
                      << "You should directly specify it via \"E\" parameter in settings.ini" << std::endl;
            std::cout << "Press ENTER button to continue" << std::endl;
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
        // std::cerr << "Debug: Found " << sys.ncenter << " atoms." << std::endl;
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
                  << std::endl;
        std::cerr << "Press ENTER to exit..." << std::endl;
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
        //           << " mass: " << mass << " amu" << std::endl;

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
    if (sys.defmass == 1 || sys.defmass == 2)
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
            //           << " total so far." << std::endl;
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

        // std::cerr << "Debug: Loaded " << sys.nfreq << " vibrational frequencies." << std::endl;
        // for (int i = 0; i < std::min(5, sys.nfreq); ++i)
        //{
        //     std::cerr << "Debug: Frequency " << (i + 1) << ": " << sys.wavenum[i] << " cm⁻¹" << std::endl;
        // }
    }
    else
    {
        std::cerr << "No vibrational frequencies found in CP2K output." << std::endl;
    }

    // ==================== Set Point Group ====================
    if (sys.PGlabelinit == "?")
    {
        if (sys.imode == 1)
        {
            sys.PGlabelinit = "C1";
            std::cout << "Note: In the case of dealing with CP2K output file with imode=1, point group is not "
                         "automatically detected but simply set to C1. "
                      << "If you need to use other point group, please manually set \"PGlabel\" in settings.ini"
                      << std::endl;
        }
    }

    file.close();
}

// ORCA

void LoadFile::loadorca(SystemData& sys)
{
    std::ifstream file(sys.inputfile);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open input file: " + sys.inputfile);
    }

    // Load energy
    // Load energy
    int ncount;
    if (!loclabelfinal(file, "FINAL SINGLE POINT ENERGY", ncount))
    {
        std::cerr << "Error: FINAL SINGLE POINT ENERGY not found in ORCA file" << std::endl;
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
            std::cerr << "Warning: Failed to parse spin multiplicity: " << e.what() << std::endl;
            sys.spinmult = 1;
        }
    }
    else
    {
        std::cerr << "Note: Ideal value S not found, assuming singlet (spin multiplicity = 1)" << std::endl;
        sys.spinmult = 1;
    }

    // Load geometry - CRITICAL FIX: Pass file handle, don't reset
    loadORCAgeom(file, sys);

    // Set mass - CRITICAL FIX: Exact pattern matching and positioning
    if (sys.defmass == 1 || sys.defmass == 2)
    {
        setatmmass(sys);
    }
    else if (sys.defmass == 3)
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
        // std::cerr << "Debug: Current file position: " << file.tellg() << std::endl;
        throw std::runtime_error("Number of atoms not found in ORCA file");
    }

    try
    {
        sys.ncenter = readaftersign_int(file, ". ");
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "Error: Failed to parse number of atoms: " << e.what() << std::endl;
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
        std::string loadArgs;
        if (!(file >> loadArgs >> sys.a[i].x >> sys.a[i].y >> sys.a[i].z))
        {
            std::cerr << "Error: Failed to read coordinates for atom " << (i + 1) << std::endl;
            throw std::runtime_error("Incomplete geometry data");
        }
        elename2idx(loadArgs, sys.a[i].index);
    }
}

void LoadFile::loadORCAfreq(std::ifstream& file, SystemData& sys)
{
    int ncount;
    if (!loclabelfinal(file, "Scaling factor for frequencies =", ncount) || ncount == 0)
    {
        std::cerr << "\nError: Unable to load frequencies from this file! Please check keywords in the ORCA input file"
                  << std::endl;
        std::cerr << "\nOr your ORCA is maybe too old (2.x). Please use the recent versions" << std::endl;
        // std::cerr << std::endl;
        // std::cerr << "Press ENTER button to exit" << std::endl;
        // std::cin.get();
        throw std::runtime_error("Unable to load frequencies from ORCA file: Scaling factor section not found");
        std::cerr << "\n" << std::endl;
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
            std::cerr << "Error: Failed to parse frequency from line: " << loadArgs << std::endl;
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
    std::ifstream file(sys.inputfile);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open input file: " + sys.inputfile);
    }

    // Load energy
    int ncount;
    if (loclabelfinal(file, "FINAL ", ncount))
    {
        std::string line;
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Could not read energy line");
        }
        size_t pos = line.find("IS ");
        if (pos != std::string::npos)
        {
            std::istringstream iss(line.substr(pos + 3));
            if (!(iss >> sys.E))
            {
                throw std::runtime_error("Failed to parse energy from: " + line);
            }
            // Extract method name (e.g., "RHF ENERGY")
            std::string method = line.substr(0, pos - 1);  // Before "IS "
            method.erase(0, method.find_first_not_of(" \t"));
            method.erase(method.find_last_not_of(" \t") + 1);
            std::cout << "Note: " << method << " energy (" << std::fixed << std::setprecision(8) << sys.E
                      << " a.u.) is loaded. If this is not the intended method, the final U, H, G may be misleading"
                      << std::endl;
        }
        else
        {
            throw std::runtime_error("Energy line format not recognized: " + line);
        }
    }
    else
    {
        throw std::runtime_error("FINAL energy section not found in GAMESS file");
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
        std::cerr << "Warning: SPIN MULTIPLICITY not found, assuming singlet (spin multiplicity = 1)" << std::endl;
        sys.spinmult = 1;
    }

    // Load geometry
    loadgmsgeom(file, sys);

    // Set mass
    file.clear();
    file.seekg(0);  // Rewind to start
    if (sys.defmass == 1 || sys.defmass == 2)
    {
        setatmmass(sys);
    }
    else if (sys.defmass == 3)
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
            // std::cerr << "Debug (Mass): Consumed label line: '" << current_line << "'" << std::endl;

            // Consume the blank line after the label.
            // After this, the file pointer is at the beginning of the first atom's mass data.
            if (!std::getline(file, current_line))
            {
                throw std::runtime_error("Failed to read blank line after 'ATOMIC WEIGHTS (AMU)'.");
            }
            // std::cerr << "Debug (Mass): Consumed blank line: '" << current_line << "'" << std::endl;


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
                    // std::cerr << "Debug (Mass): Encountered an empty/whitespace-only line, continuing." << std::endl;
                    continue;  // Skip truly empty or whitespace-only lines.
                }

                // Check if the line starts with an integer, which indicates an atom entry.
                // This is a robust way to identify the end of the mass data section.
                std::istringstream test_iss(trimmed_line);
                int                atom_num_check;
                if (!(test_iss >> atom_num_check))
                {
                    // std::cerr << "Debug (Mass): Stopping mass reading at non-atom line (e.g., 'MODES...'): '"
                    //           << trimmed_line << "'" << std::endl;
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
                //           << " amu from line: '" << trimmed_line << "'" << std::endl;
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
                        //           << expected_element_symbol << ")" << std::endl;
                        break;
                    }
                }
                if (!found)
                {
                    // Enhanced error message with both trimmed and original symbols for debugging
                    // std::cerr << "Debug (Mass): Available elements in mass_data:" << std::endl;
                    for (const auto& [sym, mass] : mass_data)
                    {
                        std::cerr << "  '" << sym << "'" << std::endl;
                    }
                    // std::cerr << "Debug (Mass): Looking for: '" << expected_element_symbol << "' (original: '"
                    //           << ind2name[sys.a[i].index] << "')" << std::endl;
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
    loadgmsfreq(file, sys);
    file.close();
}

void LoadFile::loadgmsgeom(std::ifstream& file, SystemData& sys)
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
    if (loclabelfinal(file, "COORDINATES OF ALL ATOMS ARE (ANGS)", ncount) && ncount > 0)
    {
        skiplines(file, 3);
        for (int i = 0; i < sys.ncenter; ++i)
        {
            std::string line;
            if (!std::getline(file, line))
            {
                throw std::runtime_error("Failed to read coordinates for atom " + std::to_string(i + 1));
            }
            // std::cerr << "Debug: Geometry line for atom " << i + 1 << ": " << line << std::endl;
            std::istringstream iss(line);
            std::string        loadArgs;
            double             index;  // GAMESS uses float for index (e.g., 6.0)
            if (!(iss >> loadArgs >> index >> sys.a[i].x >> sys.a[i].y >> sys.a[i].z))
            {
                throw std::runtime_error("Failed to parse coordinates for atom " + std::to_string(i + 1) +
                                         " from: " + line);
            }
            sys.a[i].index = static_cast<int>(index);  // Convert float to int
            // std::cerr << "Debug: Atom " << i + 1 << " index = " << sys.a[i].index << std::endl;
        }
    }
    else
    {
        file.clear();
        file.seekg(0);
        if (loclabel(file, "ATOMIC                      COORDINATES (BOHR)"))
        {
            skiplines(file, 2);
            for (int i = 0; i < sys.ncenter; ++i)
            {
                std::string line;
                if (!std::getline(file, line))
                {
                    throw std::runtime_error("Failed to read coordinates for atom " + std::to_string(i + 1));
                }
                std::istringstream iss(line);
                std::string        loadArgs;
                if (!(iss >> loadArgs >> sys.a[i].index >> sys.a[i].x >> sys.a[i].y >> sys.a[i].z))
                {
                    throw std::runtime_error("Failed to parse coordinates for atom " + std::to_string(i + 1) +
                                             " from: " + line);
                }
                sys.a[i].x *= b2a;
                sys.a[i].y *= b2a;
                sys.a[i].z *= b2a;
            }
        }
        else
        {
            throw std::runtime_error("No valid coordinate section found in GAMESS file");
        }
    }
}

void LoadFile::loadgmsfreq(std::ifstream& file, SystemData& sys)
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
    std::ifstream file(sys.inputfile);
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
                      << "the finally printed total U, H, G will be misleading" << std::endl;
        }
        else
        {
            std::cout << "Warning: Unable to load electronic energy, thus it is set to zero" << std::endl;
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
    loadnwgeom(file, sys);

    // Set mass
    if (sys.defmass == 1 || sys.defmass == 2)
    {
        setatmmass(sys);
    }
    else if (sys.defmass == 3)
    {
        file.clear();
        file.seekg(0);
        if (loclabel(file, "- Atom information -"))
        {
            skiplines(file, 3);
            // std::cerr << "Debug: Found 'Atom information' section, reading " << sys.ncenter << " atoms" << std::endl;
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
                    std::cerr << "Warning: Failed to parse mass line " << (atoms_read + 1) << ": " << line << std::endl;
                    parse_failed = true;
                    break;
                }

                sys.a[atoms_read].mass = mass;
                // std::cerr << "Debug: Atom " << atoms_read + 1 << " (" << element_symbol << ") mass: " << mass << "
                // amu"
                //           << std::endl;
                atoms_read++;
            }

            if (atoms_read != sys.ncenter || parse_failed)
            {
                std::cerr << "Warning: Mass section incomplete or corrupt. Falling back to default atomic masses."
                          << std::endl;
                setatmmass(sys);  // ← SAFETY: Don’t lose your atoms!
            }
        }
        else
        {
            throw std::runtime_error("'- Atom information -' section not found for mass");
        }
    }

    // Load frequencies
    loadnwfreq(file, sys);
    file.close();
}

void LoadFile::loadnwgeom(std::ifstream& file, SystemData& sys)
{
    file.clear();
    file.seekg(0);
    if (loclabel(file, "No. of atoms     :"))
    {
        // std::cerr << "DEBUG: About to read after sign..." << std::endl;
        sys.ncenter = readaftersign_int(file, ":");
        // std::cerr << "Debug: ncenter set to " << sys.ncenter << std::endl;
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
        // std::cerr << "Debug: Skipped header line: " << line << std::endl;
        //  2. Column labels: "atom    #        X              Y              Z            mass"
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Failed to read column labels line");
        }
        // std::cerr << "Debug: Skipped column labels: " << line << std::endl;
        //  3. Separator: "--------------------------------------------------------------------------"
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Failed to read separator line");
        }
        // std::cerr << "Debug: Skipped separator: " << line << std::endl;

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
                // std::cerr << "Debug: Skipping empty line before atom " << i + 1 << std::endl;
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

            // std::cerr << "Debug: Geometry line " << i + 1 << ": " << line << std::endl;
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
            // std::cerr << "Debug: Atom " << i + 1 << " (" << strtmp << ") index: " << sys.a[i].index << std::endl;
        }
    }
    else
    {
        throw std::runtime_error("'- Atom information -' section not found");
    }
}

void LoadFile::loadnwfreq(std::ifstream& file, SystemData& sys)
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

        // std::cerr << "Debug: Loaded " << sys.nfreq << " vibrational frequencies." << std::endl;
    }
    else
    {
        std::cerr << "'Projected Derivative Dipole Moments' section not found — no frequencies loaded." << std::endl;
    }
}

// XTB
void LoadFile::loadxtb(SystemData& sys)
{
    std::ifstream file(sys.inputfile);
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
                  << std::endl;

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
                                      << " Hartree" << std::endl;
                        }
                        xtbfile.close();
                        break;
                    }
                    else
                    {
                        std::cout << "Error: Unable to locate \"total energy\"! Input again" << std::endl;
                        xtbfile.close();
                        continue;
                    }
                }
                else
                {
                    std::cout << "Cannot find the file, input again!" << std::endl;
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

    if (sys.defmass == 3)
    {
        std::cout << "Note: defmass=3 is meaningless for present case because input file does not record atomic mass. "
                     "Now set atomic masses to element mass (defmass=1)"
                  << std::endl;
        sys.defmass = 1;
    }
    setatmmass(sys);

    loadGaufreq(file, sys);
    file.close();
}
