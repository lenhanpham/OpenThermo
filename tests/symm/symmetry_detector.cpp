#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <map>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cctype>

// Include symmetry headers
#include "../../cpp/symmetry.h"
#include "../../cpp/defvar.h"

namespace fs = std::filesystem;

/**
 * @brief Map of element symbols to atomic numbers
 */
std::map<std::string, int> element_to_atomic_number = {
    {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
    {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20},
    {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
    {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40},
    {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50},
    {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60},
    {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70},
    {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80},
    {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
    {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95}, {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100}
};

/**
 * @brief Trim whitespace from string
 */
std::string trim(const std::string& str) {
    auto start = std::find_if(str.begin(), str.end(), [](unsigned char c) {
        return !std::isspace(c);
    });
    auto end = std::find_if(str.rbegin(), str.rend(), [](unsigned char c) {
        return !std::isspace(c);
    }).base();
    return (start < end) ? std::string(start, end) : std::string();
}

/**
 * @brief Parse geometry file and extract atomic numbers and coordinates
 * @param filename Path to geometry file
 * @param natoms Output: number of atoms
 * @param nat Output: vector of atomic numbers
 * @param coord Output: vector of coordinates (natoms x 3)
 * @return true if parsing successful, false otherwise
 */
bool parse_geometry_file(const std::string& filename, int& natoms,
                        std::vector<int>& nat, std::vector<std::vector<double>>& coord) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }

    std::string line;
    std::vector<std::string> lines;

    // Read all lines
    while (std::getline(file, line)) {
        lines.push_back(line);
    }

    if (lines.size() < 2) {
        std::cerr << "Error: File " << filename << " has insufficient lines" << std::endl;
        return false;
    }

    // Determine format based on file extension
    bool is_xyz = filename.find(".xyz") != std::string::npos;

    std::string title;
    int num_atoms;

    if (is_xyz) {
        // XYZ format: first line is number of atoms, second line is title
        try {
            num_atoms = std::stoi(trim(lines[0]));
            title = trim(lines[1]);
        } catch (const std::exception&) {
            std::cerr << "Error: Invalid XYZ format in " << filename << std::endl;
            return false;
        }
    } else {
        // Custom format: first line is title, second line is number of atoms
        title = trim(lines[0]);
        try {
            num_atoms = std::stoi(trim(lines[1]));
        } catch (const std::exception&) {
            std::cerr << "Error: Invalid format in " << filename << std::endl;
            return false;
        }
    }

    if (num_atoms <= 0 || lines.size() < 2 + num_atoms) {
        std::cerr << "Error: Invalid number of atoms or insufficient data in " << filename << std::endl;
        return false;
    }

    natoms = num_atoms;
    nat.resize(natoms);
    coord.assign(natoms, std::vector<double>(3, 0.0));

    // Parse atom lines
    for (int i = 0; i < natoms; ++i) {
        std::string atom_line = lines[2 + i];
        std::istringstream iss(atom_line);
        std::string token;

        if (is_xyz) {
            // XYZ format: element_symbol x y z
            std::string element;
            double x, y, z;
            if (!(iss >> element >> x >> y >> z)) {
                std::cerr << "Error: Invalid atom line " << i + 1 << " in " << filename << std::endl;
                return false;
            }

            // Convert element symbol to atomic number
            auto it = element_to_atomic_number.find(element);
            if (it == element_to_atomic_number.end()) {
                std::cerr << "Error: Unknown element " << element << " in " << filename << std::endl;
                return false;
            }
            nat[i] = it->second;
            coord[i][0] = x;
            coord[i][1] = y;
            coord[i][2] = z;
        } else {
            // Custom format: atomic_number x y z
            int atomic_num;
            double x, y, z;
            if (!(iss >> atomic_num >> x >> y >> z)) {
                std::cerr << "Error: Invalid atom line " << i + 1 << " in " << filename << std::endl;
                return false;
            }
            nat[i] = atomic_num;
            coord[i][0] = x;
            coord[i][1] = y;
            coord[i][2] = z;
        }
    }

    return true;
}

/**
 * @brief Main function
 */
int main() {
    // Get current directory (should be tests/symm/)
    fs::path current_dir = fs::current_path();

    std::cout << "Symmetry Detector - Processing geometry files in: " << current_dir << std::endl;
    std::cout << "============================================================" << std::endl;

    // Iterate over all regular files in current directory
    for (const auto& entry : fs::directory_iterator(current_dir)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();

            try {
                // Parse geometry file
                int natoms;
                std::vector<int> nat;
                std::vector<std::vector<double>> coord;

                if (!parse_geometry_file(entry.path().string(), natoms, nat, coord)) {
                    std::cout << std::left << std::setw(20) << filename << ": Error parsing file" << std::endl;
                    continue;
                }

                // Create SymmetryDetector instance
                symmetry::SymmetryDetector detector;

                // Prepare output variables
                int ng, ni, nsg, ncr, nsr, np;

                // Prepare working arrays
                std::vector<std::vector<double>> symn(3, std::vector<double>(150, 0.0));
                std::vector<std::vector<int>> nsym(150, std::vector<int>(5, 0));
                int nprm;
                std::vector<std::vector<int>> nper(natoms, std::vector<int>(250, 0));
                int nseq;
                std::vector<int> nccl(natoms, 0);
                std::vector<std::vector<int>> nscl(natoms, std::vector<int>(natoms, 0));

                // Call sym_elements
                double delta = 0.01; // tolerance
                int nout = 0; // silent mode
                std::vector<std::string> symb; // empty vector for symbols

                detector.sym_elements(natoms, nat, coord, symb, delta, nout, ng, ni, nsg, ncr, nsr, np,
                                    symn, nsym, nprm, nper, nseq, nccl, nscl);

                // Get point group
                std::string point_group = detector.symm_point_group(ng, ni, nsg, ncr, nsr, np, nout);

                // Output result
                std::cout << std::left << std::setw(20) << filename << ": " << point_group << std::endl;

            } catch (const std::exception& e) {
                std::cout << std::left << std::setw(20) << filename << ": Exception - " << e.what() << std::endl;
            }
        }
    }

    std::cout << "============================================================" << std::endl;
    std::cout << "Processing complete." << std::endl;

    return 0;
}