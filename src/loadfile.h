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

//#include <vector>
#include <string>
#include "defvar.h"


/**
 * @brief Class for loading quantum chemistry output files
 * 
 * This class provides functionality to read and parse various quantum chemistry
 * output file formats including Gaussian, ORCA, CP2K, GAMESS-US, NWChem, and XTB.
 * All data is loaded into a SystemData structure to maintain modularity.
 */
class LoadFile {
private:
    // Utility functions
    static auto loclabel(std::ifstream& file, const std::string& label, int skip = 0) -> bool;
    static auto loclabelfinal(std::ifstream& file, const std::string& label, int& ncount) -> bool;
    static void skiplines(std::ifstream& file, int n);
    static auto readaftersign(std::ifstream& file, const std::string& sign) -> double;
    static auto readaftersign_int(std::ifstream& file, const std::string& sign) -> int;
    static auto readaftersign_from_line(const std::string& line, const std::string& sign) -> double;
    static void elename2idx(const std::string& element, int& index);
    static void setatmmass(SystemData& sys);
    
    // Geometry loading functions
    static void loadGaugeom(std::ifstream& file, SystemData& sys);
    static void loadORCAgeom(std::ifstream& file, SystemData& sys);
    static void loadgmsgeom(std::ifstream& file, SystemData& sys);
    static void loadnwgeom(std::ifstream& file, SystemData& sys);
    
    // Frequency loading functions
    static void loadGaufreq(std::ifstream& file, SystemData& sys);
    static void loadORCAfreq(std::ifstream& file, SystemData& sys);
    static void loadgmsfreq(std::ifstream& file, SystemData& sys);
    static void loadnwfreq(std::ifstream& file, SystemData& sys);

public:
    // Main loading functions
    static void loadotm(SystemData& sys);
    static void loadgau(SystemData& sys);
    static void loadCP2K(SystemData& sys);
    static void loadorca(SystemData& sys);
    static void loadgms(SystemData& sys);
    static void loadnw(SystemData& sys);
    static void loadxtb(SystemData& sys);
};


#endif // LOADFILE_H