/**
 * @file loadfile.h
 * @brief Header for input file parsing and loading class
 * @author Le Nhan Pham
 * @date 2025
 *
 * This header file declares functions for reading and parsing various
 * computational chemistry output file formats including Gaussian, ORCA,
 * NWChem, GAMESS, CP2K, and xTB outputs, as well as OpenThermo's native
 * .otm format.
 */

#ifndef LOADFILE_H
#define LOADFILE_H

#include "chemsys.h"
#include <string>
#include <istream>


/**
 * @brief Class for loading quantum chemistry output files
 *
 * This class provides functionality to read and parse various quantum chemistry
 * output file formats including Gaussian, ORCA, CP2K, GAMESS-US, NWChem, and XTB.
 * All data is loaded into a SystemData structure to maintain modularity.
 */
class LoadFile
{
private:
    // Utility functions
    static auto loclabel(std::istream& file, const std::string& label, int skip = 0) -> bool;
    static auto loclabelfinal(std::istream& file, const std::string& label, int& ncount) -> bool;
    static void skiplines(std::istream& file, int n);
    static auto readaftersign(std::istream& file, const std::string& sign) -> double;
    static auto readaftersign_int(std::istream& file, const std::string& sign) -> int;
    static auto readaftersign_from_line(const std::string& line, const std::string& sign) -> double;
    static void elename2idx(const std::string& element, int& index);
    static void setatmmass(SystemData& sys);

    // Geometry loading functions
    static void loadGaugeom(std::istream& file, SystemData& sys);
    static void loadORCAgeom(std::istream& file, SystemData& sys);
    static void loadGmsgeom(std::istream& file, SystemData& sys);
    static void loadNwgeom(std::istream& file, SystemData& sys);

    // Frequency loading functions
    static void loadGaufreq(std::istream& file, SystemData& sys);
    static void loadORCAfreq(std::istream& file, SystemData& sys);
    static void loadGmsfreq(std::istream& file, SystemData& sys);
    static void loadNwfreq(std::istream& file, SystemData& sys);

    // VASP loading functions
    static void loadVASPgeom(std::istream& file, SystemData& sys, bool isOUTCAR);
    static void loadVASPEnergy(std::istream& file, SystemData& sys);
    static void loadVASPfreq(std::istream& file, SystemData& sys);


public:
    // Main loading functions
    static void loadotm(SystemData& sys);
    static void loadgau(SystemData& sys);
    static void loadCP2K(SystemData& sys);
    static void loadorca(SystemData& sys);
    static void loadgms(SystemData& sys);
    static void loadnw(SystemData& sys);
    static void loadxtb(SystemData& sys);
    static void loadvasp(SystemData& sys);
};


#endif  // LOADFILE_H
