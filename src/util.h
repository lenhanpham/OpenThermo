// util.h
/**
 * @file util.h
 * @brief Header for mathematical and utility functions
 * @author Le Nhan Pham
 * @date 2025
 *
 * This header file declares various mathematical utility functions including
 * matrix operations, file parsing utilities, sorting algorithms, and string
 * manipulation functions used throughout the OpenThermo program.
 */

#ifndef UTIL_H
#define UTIL_H

#include "chemsys.h"  // For SystemData
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>


namespace util
{
    /**
     * @brief Enumeration of supported quantum chemistry programs
     */
    enum class QuantumChemistryProgram : std::uint8_t
    {
        Unknown  = 0,
        Gaussian = 1,
        Orca     = 2,
        Gamess   = 3,
        Nwchem   = 4,
        Cp2k     = 5,
        Xtb      = 6,
        Vasp     = 7,
        Otm      = 8
    };

    /**
     * @brief Multiply two matrices A (m x p) and B (p x n).
     * @param A First matrix.
     * @param B Second matrix.
     * @return Resulting matrix C (m x n).
     */
    auto matmul(const std::vector<std::vector<double>>& A,
                const std::vector<std::vector<double>>& B) -> std::vector<std::vector<double>>;

    /**
     * @brief Compute the transpose of a matrix.
     * @param M Input matrix.
     * @return Transposed matrix.
     */
    auto transpose(const std::vector<std::vector<double>>& M) -> std::vector<std::vector<double>>;

    /**
     * @brief Diagonalize a symmetric matrix using the Jacobi method.
     * @param mat Input matrix (will be modified to diagonal form).
     * @param S Output eigenvector matrix (columns are eigenvectors).
     * @param eigval Output eigenvalue vector.
     * @param maxcyc Maximum number of cycles (default: 200).
     * @param thres Convergence threshold (default: 1e-9).
     */
    void diagmat(std::vector<std::vector<double>>& mat,
                 std::vector<std::vector<double>>& S,
                 std::array<double, 3>&            eigval,
                 int                               maxcyc = 200,
                 double                            thres  = 1e-9);

    ///**
    // * @brief Skip a specific number of lines in the input stream.
    // * @param file Input stream.
    // * @param nskip Number of lines to skip.
    // */
    // void skiplines(std::istream& file, int nskip);

    ///**
    // * @brief Read a double value after the last occurrence of a sign in the current line.
    // * @param file Input stream.
    // * @param sign Sign to search for.
    // * @param val Output value.
    // */
    // void readaftersign(std::istream& file, const std::string& sign, double& val);

    ///**
    // * @brief Read an integer value after the last occurrence of a sign in the current line.
    // * @param file Input stream.
    // * @param sign Sign to search for.
    // * @param val Output value.
    // */
    // void readaftersign_int(std::istream& file, const std::string& sign, int& val);


    /**
     * @brief Sort a double array from small to large (bubble sort).
     * @param array Array to sort.
     * @param mode "val" to sort by value, "abs" to sort by absolute value (default: "val").
     * @param list Optional index list to permute.
     * @param list2 Optional second index list to permute.
     */
    // void sortr8(std::vector<double>& array,
    //             const std::string&   mode  = "val",
    //             std::vector<int>*    list  = nullptr,
    //             std::vector<int>*    list2 = nullptr);

    /**
     * @brief Parse command-line arguments and update system parameters
     */
    void loadarguments(SystemData& sys, int argc, std::vector<std::string>& argv);

    /**
     * @brief Extract string value for a command-line option
     */
    void get_option_str(std::ifstream& file, const std::string& key, std::string& value);

    /**
     * @brief Load calculation settings from settings.ini file
     */
    void loadsettings(SystemData& sys);

    /**
     * @brief Locate a specific label/string in an input file
     */
    auto loclabel(std::ifstream&     file,
                  const std::string& label,
                  int&               nskip,
                  bool               rewind    = true,
                  bool               find_last = false,
                  int                maxline   = 0) -> bool;

    /**
     * @brief Modify atomic masses based on user settings
     */
    void modmass(SystemData& sys);

    /**
     * @brief Determine the computational chemistry program that generated the input
     */
    auto deterprog(SystemData& sys) -> QuantumChemistryProgram;

    /**
     * @brief Output molecular data to .otm format file
     */
    void outotmfile(SystemData& sys);

    /**
     * @brief Create a default settings.ini file with all default parameters
     */
    void create_default_settings_file();

}  // namespace util

#endif  // UTIL_H