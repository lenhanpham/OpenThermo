// util.cpp
/**
 * @file util.cpp
 * @brief Implementation of mathematical and utility functions
 * @author Le Nhan Pham
 * @date 2025
 *
 * This file contains the implementation of various mathematical utility
 * functions including matrix operations, file parsing utilities, sorting
 * algorithms, and string manipulation functions used throughout OpenThermo.
 */

#include "util.h"
#include "chemsys.h"  // For SystemData and related types
#include <algorithm>  // for std::swap
#include <array>
#include <cmath>  // for atan, sin, cos, abs
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace util
{

    /**
     * @brief Perform matrix multiplication A * B (UNUSED - consider removal)
     *
     * Computes the matrix
     * product of two matrices A (m
     * x p) and B (p x n),
     * resulting in matrix C (m x n).
     *
     *
     * @param A First matrix (m x p)
     *
     * @param B Second matrix (p x n)
     * @return Result matrix C (m x
     * n)
     *
     * @throws std::runtime_error if
     * matrices have incompatible dimensions or are empty
     *
     * @note Uses standard matrix multiplication algorithm
     * with O(m*n*p) complexity
     * @note UNUSED: Only
     * called internally within util.cpp - no external usage
     */
    auto matmul(const std::vector<std::vector<double>>& A,
                const std::vector<std::vector<double>>& B) -> std::vector<std::vector<double>>
    {
        size_t m = A.size();
        if (m == 0)
            throw std::runtime_error("matmul: Empty matrix A");
        size_t p = A[0].size();
        size_t n = B[0].size();
        if (B.size() != p)
            throw std::runtime_error("matmul: Dimension mismatch");

        std::vector<std::vector<double>> C(m, std::vector<double>(n, 0.0));
        for (size_t i = 0; i < m; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                for (size_t k = 0; k < p; ++k)
                {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }

    /**
     * @brief Compute the transpose of a matrix (UNUSED - consider removal)
     *
     * Returns the transpose
     * of the input matrix M, where
     * element M[i][j]
     * becomes element M[j][i] in the result.
     *
     *
     * @param M Input matrix to transpose

     * * @return Transposed matrix
     *
     * @note Returns empty matrix
     * if input is empty
     * @note Assumes
     * rectangular matrix (all rows have same number of columns)
     *
     * @note UNUSED: Only called internally within util.cpp - no external usage
     */
    auto transpose(const std::vector<std::vector<double>>& M) -> std::vector<std::vector<double>>
    {
        if (M.empty())
            return {};
        size_t                           rows = M.size();
        size_t                           cols = M[0].size();
        std::vector<std::vector<double>> T(cols, std::vector<double>(rows, 0.0));
        for (size_t i = 0; i < rows; ++i)
        {
            for (size_t j = 0; j < cols; ++j)
            {
                T[j][i] = M[i][j];
            }
        }
        return T;
    }

    /**
     * @brief Diagonalize a symmetric matrix using the Jacobi method
     *
     * Computes eigenvalues and
     * eigenvectors of a symmetric matrix using the
     * classical Jacobi eigenvalue algorithm. The input matrix is
     * transformed
     * to diagonal form, with eigenvalues stored in eigval and eigenvectors
     * in the columns of
     * matrix S.
     *
     * @param mat [in/out] Input symmetric matrix (modified to diagonal form)
     * @param S
     * [out] Eigenvector matrix (columns are eigenvectors)
     * @param eigval [out] Eigenvalue array (first 3 elements
     * used)
     * @param maxcyc Maximum number of Jacobi rotations (default: 200)
     * @param thres Convergence
     * threshold for off-diagonal elements (default: 1e-9)
     *
     * @throws std::runtime_error if matrix is not
     * square or empty
     * @note Uses Jacobi rotations to iteratively reduce off-diagonal elements
     * @note
     * Particularly efficient for small matrices (3x3 in OpenThermo's case)
     */
    void diagmat(std::vector<std::vector<double>>& mat,
                 std::vector<std::vector<double>>& S,
                 std::array<double, 3>&            eigval,
                 int                               maxcyc,
                 double                            thres)
    {
        size_t n = mat.size();
        if (n == 0 || mat[0].size() != n)
            throw std::runtime_error("diagmat: Matrix must be square and non-empty");

        S.assign(n, std::vector<double>(n, 0.0));
        for (size_t i = 0; i < n; ++i)
            S[i][i] = 1.0;

        for (int k = 0; k <= maxcyc; ++k)
        {
            std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));
            for (size_t i = 0; i < n; ++i)
                R[i][i] = 1.0;

            size_t max_i = 0, max_j = 1;
            double max_offdiag = 0.0;
            for (size_t ii = 0; ii < n; ++ii)
            {
                for (size_t jj = ii + 1; jj < n; ++jj)
                {
                    double abs_val = std::abs(mat[ii][jj]);
                    if (abs_val > max_offdiag)
                    {
                        max_offdiag = abs_val;
                        max_i       = ii;
                        max_j       = jj;
                    }
                }
            }

            if (max_offdiag < thres)
                break;
            if (k == maxcyc)
                std::cerr << "Note: Matrix diagonalization exceeded max cycle before convergence" << "\n";

            double phi      = std::atan(2.0 * mat[max_i][max_j] / (mat[max_i][max_i] - mat[max_j][max_j])) / 2.0;
            double c        = std::cos(phi);
            double s        = std::sin(phi);
            R[max_i][max_i] = c;
            R[max_j][max_j] = c;
            R[max_i][max_j] = -s;
            R[max_j][max_i] = s;

            auto Rt = transpose(R);
            mat     = matmul(matmul(Rt, mat), R);
            S       = matmul(S, R);
        }

        for (size_t i = 0; i < 3; ++i)
        {
            eigval[i] = mat[i][i];
        }
    }


    // void loclabel(std::istream& file, const std::string& label, bool& found, bool rewind, int maxline)
    // {
    //     if (rewind)
    //     {
    //         file.clear();
    //         file.seekg(0);
    //     }
    //
    //     std::string line;
    //     int         linecount = 0;
    //     found                 = false;
    //
    //     while (std::getline(file, line))
    //     {
    //         if (line.find(label) != std::string::npos)
    //         {
    //             // Backspace: seek back to start of this line
    //             file.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
    //             found = true;
    //             return;
    //         }
    //         ++linecount;
    //         if (maxline > 0 && linecount >= maxline)
    //             break;
    //     }
    // }

    // void loclabelfinal(std::istream& file, const std::string& label, int& nfound)
    // {
    //     nfound = 0;
    //     file.clear();
    //     file.seekg(0);
    //
    //     bool found;
    //     while (true)
    //     {
    //         loclabel(file, label, found, false);
    //         if (!found)
    //             break;
    //         ++nfound;
    //         std::string dummy;
    //         std::getline(file, dummy);  // Skip the line
    //     }
    //
    //     file.clear();
    //     file.seekg(0);
    //     for (int i = 0; i < nfound; ++i)
    //     {
    //         loclabel(file, label, found, false);
    //         if (i < nfound - 1)
    //         {
    //             std::string dummy;
    //             std::getline(file, dummy);
    //         }
    //     }
    // }

    // void get_option_str(std::istream& file, const std::string& label, std::string& str)
    // {
    //     bool found;
    //     loclabel(file, label, found);
    //     str.clear();
    //     if (found)
    //     {
    //         std::string line;
    //         std::getline(file, line);
    //
    //         // Remove all '"' characters
    //         line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
    //
    //         size_t ibeg = line.find('=');
    //         if (ibeg == std::string::npos)
    //             return;
    //
    //         size_t iend = line.find("//");
    //         if (iend == std::string::npos)
    //             iend = line.length();
    //
    //         str = line.substr(ibeg + 1, iend - ibeg - 1);
    //
    //         // Trim leading whitespace
    //         size_t start = str.find_first_not_of(" \t");
    //         if (start != std::string::npos)
    //         {
    //             str = str.substr(start);
    //         }
    //         else
    //         {
    //             str.clear();
    //         }
    //     }
    // }

    // void sortr8(std::vector<double>& array, const std::string& mode, std::vector<int>* list, std::vector<int>* list2)
    // {
    //     size_t N        = array.size();
    //     bool   byabs    = (mode == "abs");
    //     bool   haslist  = (list != nullptr && list->size() == N);
    //     bool   haslist2 = (list2 != nullptr && list2->size() == N);
    //
    //     for (size_t i = 0; i < N; ++i)
    //     {
    //         for (size_t j = i + 1; j < N; ++j)
    //         {
    //             bool swap_needed = false;
    //             if (byabs)
    //             {
    //                 if (std::abs(array[i]) > std::abs(array[j]))
    //                     swap_needed = true;
    //             }
    //             else
    //             {
    //                 if (array[i] > array[j])
    //                     swap_needed = true;
    //             }
    //             if (swap_needed)
    //             {
    //                 std::swap(array[i], array[j]);
    //                 if (haslist)
    //                     std::swap((*list)[i], (*list)[j]);
    //                 if (haslist2)
    //                     std::swap((*list2)[i], (*list2)[j]);
    //             }
    //         }
    //     }
    // }

    /**
     * @brief Parse low vibrational frequency treatment method
     *
     * Accepts both integer values (0-3)
     * and string names for backward compatibility
     * and user convenience.
     *
     * @param str Input string to
     * parse
     * @return LowVibTreatment enum value
     * @throws std::runtime_error if input is invalid
     */
    LowVibTreatment parseLowVibTreatment(const std::string& str)
    {
        // Try to parse as integer first
        std::istringstream iss(str);
        int                intValue;
        if (iss >> intValue)
        {
            switch (intValue)
            {
                case 0:
                    return LowVibTreatment::Harmonic;
                case 1:
                    return LowVibTreatment::Truhlar;
                case 2:
                    return LowVibTreatment::Grimme;
                case 3:
                    return LowVibTreatment::Minenkov;
                default:
                    throw std::runtime_error("Invalid low frequency treatment value: " + str +
                                             ". Must be 0-3 or method name.");
            }
        }

        // Try to parse as string (case-insensitive)
        std::string lowVibMth = str;
        std::transform(lowVibMth.begin(), lowVibMth.end(), lowVibMth.begin(), ::tolower);

        if (lowVibMth == "harmonic")
            return LowVibTreatment::Harmonic;
        if (lowVibMth == "truhlar")
            return LowVibTreatment::Truhlar;
        if (lowVibMth == "grimme")
            return LowVibTreatment::Grimme;
        if (lowVibMth == "minenkov")
            return LowVibTreatment::Minenkov;

        throw std::runtime_error("Invalid low frequency treatment method: " + str +
                                 ". Valid options: 0/Harmonic, 1/Truhlar, 2/Grimme, 3/Minenkov");
    }

    /**
     * @brief Parse command-line arguments and update system parameters
     *
     * Processes command-line
     * arguments passed to
     * OpenThermo and updates the
     * SystemData structure with user-specified parameters.
     * Supports various
     * options for
     * temperature, pressure, scaling factors, and other calculation
     *
     * settings.
     *
     * @param sys SystemData structure to
     * update with parsed arguments
     * @param argc
     * Number of command-line arguments
     * @param argv Array of command-line
     * argument strings
     *
     *
     * @throws std::runtime_error if invalid arguments or missing values are encountered
     * @note
     * Arguments
     * start from index 2 (after program name and input file)
     * @note Supports options like -T, -P, -E, -prtvib,

     * * -sclZPE, etc.
     */
    void loadarguments(SystemData& sys, int argc, std::vector<std::string>& argv)
    {
        // Print full command line
        std::string command = argv[0];
        for (int i = 1; i < argc; ++i)
        {
            command += " " + std::string(argv[i]);
        }
        std::cout << "Command of invoking OpenThermo:\n " << command << "\n";

        // Check if additional arguments override parameters
        if (argc > 2)
        {  // argv[0] = program, argv[1] = inputfile
            std::cout << "Note: One or more running parameters are overridden by arguments\n";
        }

        // Parse arguments starting from index 2
        for (int iarg = 2; iarg < argc; ++iarg)
        {
            const std::string& inputArgs = argv[iarg];
            if (inputArgs == "-E")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -E");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.Eexter))
                    throw std::runtime_error("Error: Invalid value for -E");
            }
            else if (inputArgs == "-prtvib")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -prtvib");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.prtvib))
                    throw std::runtime_error("Error: Invalid value for -prtvib");
            }
            else if (inputArgs == "-T")
            {
                // Check if next 3 arguments exist for scanning
                if (iarg + 3 < argc)
                {
                    std::istringstream iss_low(argv[iarg + 1]);
                    std::istringstream iss_high(argv[iarg + 2]);
                    std::istringstream iss_step(argv[iarg + 3]);
                    if (iss_low >> sys.Tlow && iss_high >> sys.Thigh && iss_step >> sys.Tstep)
                    {
                        iarg += 3;  // Skip the consumed arguments
                        continue;
                    }
                }
                // Fallback to single temperature
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -T");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.T))
                    throw std::runtime_error("Error: Invalid value for -T");
            }
            else if (inputArgs == "-P")
            {
                // Check if next 3 arguments exist for scanning
                if (iarg + 3 < argc)
                {
                    std::istringstream iss_low(argv[iarg + 1]);
                    std::istringstream iss_high(argv[iarg + 2]);
                    std::istringstream iss_step(argv[iarg + 3]);
                    if (iss_low >> sys.Plow && iss_high >> sys.Phigh && iss_step >> sys.Pstep)
                    {
                        iarg += 3;  // Skip the consumed arguments
                        continue;
                    }
                }
                // Fallback to single pressure
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -P");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.P))
                    throw std::runtime_error("Error: Invalid value for -P");
            }
            else if (inputArgs == "-sclZPE")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -sclZPE");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.sclZPE))
                    throw std::runtime_error("Error: Invalid value for -sclZPE");
            }
            else if (inputArgs == "-sclheat")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -sclheat");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.sclheat))
                    throw std::runtime_error("Error: Invalid value for -sclheat");
            }
            else if (inputArgs == "-sclS")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -sclS");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.sclS))
                    throw std::runtime_error("Error: Invalid value for -sclS");
            }
            else if (inputArgs == "-sclCV")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -sclCV");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.sclCV))
                    throw std::runtime_error("Error: Invalid value for -sclCV");
            }
            else if (inputArgs == "-lowvibmeth")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -lowvibmeth");
                try
                {
                    sys.lowVibTreatment = parseLowVibTreatment(argv[iarg]);
                }
                catch (const std::runtime_error& e)
                {
                    throw std::runtime_error("Error: " + std::string(e.what()));
                }
            }
            else if (inputArgs == "-ravib")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -ravib");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.ravib))
                    throw std::runtime_error("Error: Invalid value for -ravib");
            }
            else if (inputArgs == "-ipmode")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -ipmode");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.ipmode))
                    throw std::runtime_error("Error: Invalid value for -ipmode");
            }
            else if (inputArgs == "-imagreal")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -imagreal");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.imagreal))
                    throw std::runtime_error("Error: Invalid value for -imagreal");
            }
            else if (inputArgs == "-conc")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -conc");
                sys.concstr = argv[iarg];
            }
            else if (inputArgs == "-outotm")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -outotm");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.outotm))
                    throw std::runtime_error("Error: Invalid value for -outotm");
            }
            else if (inputArgs == "-massmod")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -massmod");
                std::istringstream iss(argv[iarg]);
                if (!(iss >> sys.massmod))
                    throw std::runtime_error("Error: Invalid value for -massmod");
            }
            else if (inputArgs == "-PGname")
            {
                if (++iarg >= argc)
                    throw std::runtime_error("Error: Missing value for -PGname");
                sys.PGnameinit = argv[iarg];
            }
            else if (inputArgs == "-noset")
            {
                continue;
            }
            else
            {
                std::cerr << "Error: Unable to recognize argument " << inputArgs << "\n";
                std::exit(1);
            }
        }
    }

    // Utility function to read a string after a key in settings.ini
    void get_option_str(std::ifstream& file, const std::string& key, std::string& value)
    {
        file.clear();
        file.seekg(0);
        std::string line;
        value = "";
        while (std::getline(file, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            // Find the key at the beginning of the line or after whitespace
            size_t pos = 0;
            // Skip leading whitespace
            while (pos < line.length() && std::isspace(line[pos]))
                pos++;
            // Check if the key matches at this position
            if (pos + key.length() <= line.length() && line.substr(pos, key.length()) == key)
            {
                // Found the key, extract everything after it
                value = line.substr(pos + key.length());
                return;
            }
        }
    }

    /**
     * @brief Load calculation settings from settings.ini file
     *
     * Reads the settings.ini
     * configuration file and
     * updates the SystemData
     * structure with user-defined parameters and
     * preferences. The file is
     * searched in the
     * current directory first, then in the openthermopath
     *
     * environment variable if set.
     *
     * @param sys SystemData
     * structure to update with settings
     *

     * * @note If settings.ini is not found, default values are used
     * @note Supports
     * parameters like
     * temperature, pressure, scaling factors, etc.
     * @note Environment variable OPENTHERMOPATH can specify
     *
     * settings location
     */
    void loadsettings(SystemData& sys)
    {
        std::string   settingpath = "settings.ini";
        std::ifstream file(settingpath);
        if (!file.is_open())
        {
            char* openthermopath = std::getenv("openthermopath");
            if (openthermopath)
            {
                settingpath = std::string(openthermopath) +
                              (std::filesystem::path::preferred_separator == '\\' ? "\\settings.ini" : "/settings.ini");
                file.open(settingpath);
            }
        }
        if (file.is_open())
        {
            std::cout << "\nLoading running parameters from settings.ini...\n";
            std::string inputArgs;
            get_option_str(file, "E", inputArgs);
            if (!inputArgs.empty())
            {
                std::istringstream iss(inputArgs);
                iss >> std::ws;
                iss.ignore(1, '=');
                if (!(iss >> sys.Eexter))
                    throw std::runtime_error("Error: Invalid value for E in settings.ini");
            }
            get_option_str(file, "prtvib", inputArgs);
            if (!inputArgs.empty())
            {
                std::istringstream iss(inputArgs);
                iss >> std::ws;
                iss.ignore(1, '=');
                if (!(iss >> sys.prtvib))
                    throw std::runtime_error("Error: Invalid value for prtvib in settings.ini");
            }
            get_option_str(file, "T", inputArgs);
            if (!inputArgs.empty())
            {
                std::istringstream iss(inputArgs);
                iss.ignore(1, '=');
                iss >> std::ws;
                if (!(iss >> sys.Tlow >> sys.Thigh >> sys.Tstep))
                {
                    // Try to parse as single temperature value
                    std::string tempStr = inputArgs.substr(inputArgs.find('=') + 1);
                    // Remove leading/trailing whitespace
                    size_t start = tempStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        tempStr = tempStr.substr(start);
                    size_t end = tempStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        tempStr = tempStr.substr(0, end + 1);

                    std::istringstream tempIss(tempStr);
                    if (!(tempIss >> sys.T))
                        throw std::runtime_error("Error: Invalid value for T in settings.ini");
                }
            }
            get_option_str(file, "P", inputArgs);
            if (!inputArgs.empty())
            {
                std::istringstream iss(inputArgs);
                iss >> std::ws;
                iss.ignore(1, '=');
                if (!(iss >> sys.Plow >> sys.Phigh >> sys.Pstep))
                {
                    iss.clear();
                    iss.seekg(0);
                    iss >> std::ws;
                    iss.ignore(1, '=');
                    if (!(iss >> sys.P))
                        throw std::runtime_error("Error: Invalid value for P in settings.ini");
                }
            }
            get_option_str(file, "sclZPE", inputArgs);
            if (!inputArgs.empty())
            {
                std::istringstream iss(inputArgs);
                iss >> std::ws;
                iss.ignore(1, '=');
                if (!(iss >> sys.sclZPE))
                    throw std::runtime_error("Error: Invalid value for sclZPE in settings.ini");
            }
            get_option_str(file, "sclheat", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    std::istringstream iss(valueStr);
                    if (!(iss >> sys.sclheat))
                        throw std::runtime_error("Error: Invalid value for sclheat in settings.ini");
                }
            }
            get_option_str(file, "sclS", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    std::istringstream iss(valueStr);
                    if (!(iss >> sys.sclS))
                        throw std::runtime_error("Error: Invalid value for sclS in settings.ini");
                }
            }
            get_option_str(file, "sclCV", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    std::istringstream iss(valueStr);
                    if (!(iss >> sys.sclCV))
                        throw std::runtime_error("Error: Invalid value for sclCV in settings.ini");
                }
            }
            get_option_str(file, "lowvibmeth", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    try
                    {
                        sys.lowVibTreatment = parseLowVibTreatment(valueStr);
                    }
                    catch (const std::runtime_error& e)
                    {
                        throw std::runtime_error("Error: " + std::string(e.what()) + " in settings.ini");
                    }
                }
            }
            get_option_str(file, "ravib", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    std::istringstream iss(valueStr);
                    if (!(iss >> sys.ravib))
                        throw std::runtime_error("Error: Invalid value for ravib in settings.ini");
                }
            }
            get_option_str(file, "intpvib", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    std::istringstream iss(valueStr);
                    if (!(iss >> sys.intpvib))
                        throw std::runtime_error("Error: Invalid value for intpvib in settings.ini");
                }
            }
            get_option_str(file, "imagreal", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    std::istringstream iss(valueStr);
                    if (!(iss >> sys.imagreal))
                        throw std::runtime_error("Error: Invalid value for imagreal in settings.ini");
                }
            }
            get_option_str(file, "ipmode", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    std::istringstream iss(valueStr);
                    if (!(iss >> sys.ipmode))
                        throw std::runtime_error("Error: Invalid value for ipmode in settings.ini");
                }
            }
            get_option_str(file, "conc", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    sys.concstr = valueStr;
                }
            }
            get_option_str(file, "outotm", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    std::istringstream iss(valueStr);
                    if (!(iss >> sys.outotm))
                        throw std::runtime_error("Error: Invalid value for outotm in settings.ini");
                }
            }
            get_option_str(file, "PGname", inputArgs);
            if (!inputArgs.empty())
                sys.PGnameinit = inputArgs;
            get_option_str(file, "massmod", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    std::istringstream iss(valueStr);
                    if (!(iss >> sys.massmod))
                        throw std::runtime_error("Error: Invalid value for massmod in settings.ini");
                }
            }
            get_option_str(file, "extrape", inputArgs);
            if (!inputArgs.empty())
            {
                // Find the position of '='
                size_t eqPos = inputArgs.find('=');
                if (eqPos != std::string::npos)
                {
                    std::string valueStr = inputArgs.substr(eqPos + 1);
                    // Remove leading/trailing whitespace
                    size_t start = valueStr.find_first_not_of(" \t");
                    if (start != std::string::npos)
                        valueStr = valueStr.substr(start);
                    size_t end = valueStr.find_last_not_of(" \t");
                    if (end != std::string::npos)
                        valueStr = valueStr.substr(0, end + 1);

                    // Convert to lowercase for case-insensitive comparison
                    std::transform(valueStr.begin(), valueStr.end(), valueStr.begin(), ::tolower);

                    if (valueStr == "true" || valueStr == "yes" || valueStr == "1")
                    {
                        sys.vasp_energy_select = 1;
                    }
                    else if (valueStr == "false" || valueStr == "no" || valueStr == "0")
                    {
                        sys.vasp_energy_select = 0;
                    }
                    else
                    {
                        throw std::runtime_error(
                            "Error: Invalid value for extrape in settings.ini. Use true/yes/1 or false/no/0");
                    }
                }
            }
            file.close();
        }
        else
        {
            std::cout << "\nWarning: settings.ini could not be found in either current directory or the directory defined by "
                         "openthermopath environment, "
                      << "thus default parameters are used!\n";
        }
    }

    /**
     * @brief Advanced label locator with positioning control (PRIMARY FUNCTION)
     *
     * This is the main
     * label-finding
     * function used throughout OpenThermo for parsing
     * various quantum chemistry output file
     * formats (Gaussian, ORCA,
     * NWChem, etc.).
     *
     * Key Features:
     * - Flexible search options
     * (rewind, find_last, maxline)
     * - Returns count of
     * skipped occurrences (nskip)
     * - Positions
     * stream at BEGINNING of target line
     * - Supports both first and last
     * occurrence finding
     * - Used
     * in program detection (deterprog) and file parsing
     *
     * Positioning Behavior:
     * - After
     *
     * finding label, calculates position at line START
     * - Useful for re-reading the same line for data
     * extraction
     * - Enables subsequent getline() calls to parse line content
     *
     * Usage Patterns:
     *
     * // Find first occurrence within
     * 200 lines
     *   if (loclabel(file, "Energy ", nskip, true, false,
     * 200)) {
     *       double energy;
     *       file >>
     * energy;  // Read value from same line
     *   }

     * *
     *   // Find last occurrence of section header
     *   if
     * (loclabel(file, "Geometry", nskip, true,
     * true, 0)) {
     *       // Process geometry section
     *   }
     *
     * Parameters:
     *
     * - rewind:
     * If true, reset to file beginning before search
     * - find_last: If true, find last occurrence instead of

     * * first
     * - maxline: Maximum lines to search (0 = unlimited)
     * - nskip: Returns count of occurrences
     * found (useful
     * for find_last)
     *
     * Status: HEAVILY USED - 8 calls in sub.cpp alone
     */
    auto loclabel(std::ifstream& file, const std::string& label, int& nskip, bool rewind, bool find_last, int maxline)
        -> bool
    {
        std::string line;
        if (rewind)
        {
            file.clear();
            file.seekg(0);
        }
        nskip                     = 0;
        std::streampos last_pos   = file.tellg();
        int            line_count = 0;
        while (std::getline(file, line) && (maxline == 0 || line_count < maxline))
        {
            if (line.find(label) != std::string::npos)
            {
                nskip++;
                // Calculate position at beginning of found line
                std::streampos current   = file.tellg();
                std::streampos begin_pos = current - static_cast<std::streamoff>(line.length() + 1);  // assuming \n
                if (!find_last)
                {
                    file.seekg(begin_pos);  // Position at beginning of found line
                    return true;
                }
                last_pos = begin_pos;
            }
            line_count++;
        }
        if (nskip > 0 && find_last)
        {
            file.clear();
            file.seekg(last_pos);  // Return to last occurrence
            return true;
        }
        if (rewind)
        {
            file.clear();
            file.seekg(0);
        }
        return false;
    }

    /**
     * @brief Apply user-specified atomic mass modifications
     *
     * Reads mass modification settings from
     * settings.ini
     * and applies them to
     * specific atoms. Supports both explicit mass values and isotope
     * specifications.
     * Used for
     * isotopic substitution studies and custom mass assignments.
     *
     *
     * @param sys SystemData structure containing atomic
     * data to modify
     *
     * @note Reads from "modmass"
     * section in settings.ini
     * @note Format: atom_index isotope_number or
     * atom_index explicit_mass
     *
     * @note Isotope numbers are element-specific (1-based indexing)
     * @note Only modifies
     * atoms specified
     * in the settings file
     */
    void modmass(SystemData& sys)
    {
        std::ifstream file("settings.ini");
        if (file.is_open())
        {
            int nskip;
            if (loclabel(file, "modmass", nskip))
            {
                std::string line;
                std::getline(file, line);  // Skip modmass line
                while (std::getline(file, line))
                {
                    if (line.empty() || line[0] == '#')
                        continue;
                    std::istringstream iss(line);
                    int                iatm;
                    std::string        inputArgs;
                    if (!(iss >> iatm >> inputArgs))
                    {
                        break;  // Exit on read failure (end of section)
                    }
                    if (iatm < 1 || iatm > sys.ncenter)
                    {
                        throw std::runtime_error("Error: Invalid atom index " + std::to_string(iatm));
                    }
                    if (inputArgs.find('.') == std::string::npos)
                    {
                        // Integer (isotope)
                        int                isotmp;
                        std::istringstream iss_iso(inputArgs);
                        if (!(iss_iso >> isotmp))
                        {
                            throw std::runtime_error("Error: Invalid isotope for atom " + std::to_string(iatm));
                        }
                        sys.a[iatm - 1].mass = isomass[sys.a[iatm - 1].index][isotmp];
                        if (sys.a[iatm - 1].mass == 0.0)
                        {
                            std::cerr << "Error: No isotope " << isotmp << " for atom " << iatm << "!\n";
                            std::cerr << "Exiting program\n";
                            file.close();
                            std::exit(1);
                        }
                    }
                    else
                    {
                        // Explicit mass
                        std::istringstream iss_mass(inputArgs);
                        if (!(iss_mass >> sys.a[iatm - 1].mass))
                        {
                            throw std::runtime_error("Error: Invalid mass for atom " + std::to_string(iatm));
                        }
                    }
                    std::cout << "Mass of atom " << iatm << " (" << ind2name[sys.a[iatm - 1].index]
                              << ") has been modified to " << std::fixed << std::setprecision(6) << sys.a[iatm - 1].mass
                              << " amu\n";
                }
            }
            file.close();
        }
    }

    // Determine the program that generated the input file
    auto deterprog(SystemData& sys) -> QuantumChemistryProgram
    {
        std::ifstream file(sys.inputfile);
        if (!file.is_open())
        {
            throw std::runtime_error("Error: Could not open file " + sys.inputfile);
        }

        int nskip;
        if (loclabel(file, "generated by the xtb code", nskip, true, false, 200))
        {
            file.close();
            return QuantumChemistryProgram::Xtb;  // xtb g98.out
        }
        if (loclabel(file, "Gaussian, Inc", nskip, true, false, 200) ||
            loclabel(file, "Entering Gaussian System", nskip, true, false, 200))
        {
            file.close();
            return QuantumChemistryProgram::Gaussian;  // Gaussian
        }
        if (loclabel(file, "O   R   C   A", nskip, true, false, 200))
        {
            file.close();
            return QuantumChemistryProgram::Orca;  // ORCA
        }
        if (loclabel(file, "GAMESS", nskip, true, false, 200))
        {
            file.close();
            return QuantumChemistryProgram::Gamess;  // GAMESS-US
        }
        if (loclabel(file, "Northwest Computational Chemistry Package", nskip, true, false, 200))
        {
            file.close();
            return QuantumChemistryProgram::Nwchem;  // NWChem
        }
        if (loclabel(file, "CP2K|", nskip, true, false, 200))
        {
            file.close();
            return QuantumChemistryProgram::Cp2k;  // CP2K
        }
        if (loclabel(file, "vasp", nskip, true, false, 200))
        {
            file.close();
            return QuantumChemistryProgram::Vasp;  // VASP
        }
        file.close();
        return QuantumChemistryProgram::Unknown;  // Undetermined
    }

    // Output data to .otm file
    void outotmfile(SystemData& sys)
    {
        if (sys.inputfile.empty())
        {
            throw std::runtime_error("Error: inputfile is empty");
        }
        size_t itmp = sys.inputfile.rfind('.');
        if (itmp == std::string::npos)
        {
            throw std::runtime_error("Error: inputfile has no extension");
        }
        std::string otmpath = sys.inputfile.substr(0, itmp) + ".otm";
        std::cout << "Outputting data to " << otmpath << "\n";

        std::ofstream file(otmpath, std::ios::out);
        if (!file.is_open())
        {
            throw std::runtime_error("Error: Could not open file " + otmpath + " for writing");
        }

        // Write electronic energy
        file << "*E  //Electronic energy (a.u.)\n";
        file << std::fixed << std::setprecision(10) << std::setw(20) << sys.E << "\n";

        // Write wavenumbers
        file << "*wavenum  //Wavenumbers (cm-1).\n";
        if (sys.nfreq != static_cast<int>(sys.wavenum.size()))
        {
            file.close();
            throw std::runtime_error("Error: nfreq does not match wavenum size");
        }
        for (int ifreq = 0; ifreq < sys.nfreq; ++ifreq)
        {
            file << std::fixed << std::setprecision(4) << std::setw(10) << sys.wavenum[ifreq] << "\n";
        }

        // Write atoms
        file << "*atoms  //System infor: Name, mass (amu), X, Y, Z (Angstrom)\n";
        if (sys.ncenter != static_cast<int>(sys.a.size()))
        {
            file.close();
            throw std::runtime_error("Error: ncenter does not match atom array size");
        }
        for (int iatm = 0; iatm < sys.ncenter; ++iatm)
        {
            if (sys.a[iatm].index < 1 || sys.a[iatm].index > nelesupp)
            {
                file.close();
                throw std::runtime_error("Error: Invalid element index for atom " + std::to_string(iatm + 1));
            }
            file << std::left << std::setw(4) << ind2name[sys.a[iatm].index] << std::fixed << std::setprecision(6)
                 << std::setw(12) << sys.a[iatm].mass << std::setw(12) << sys.a[iatm].x << std::setw(12)
                 << sys.a[iatm].y << std::setw(12) << sys.a[iatm].z << "\n";
        }

        // Write energy levels
        file << "*elevel  //Energy (eV) and degeneracy of electronic energy levels\n";
        if (sys.nelevel != static_cast<int>(sys.elevel.size()) || sys.nelevel != static_cast<int>(sys.edegen.size()))
        {
            file.close();
            throw std::runtime_error("Error: nelevel does not match elevel or edegen size");
        }
        for (int ie = 0; ie < sys.nelevel; ++ie)
        {
            file << std::fixed << std::setprecision(6) << std::setw(12) << sys.elevel[ie] << std::setw(4)
                 << static_cast<int>(sys.edegen[ie]) << "\n";
        }

        file.close();
        std::cout << " " << otmpath << " has been successfully generated!\n";
    }

    // Create default settings.ini file with all default parameters
    void create_default_settings_file()
    {
        std::string   filename = "settings.ini";
        std::ofstream file(filename);

        if (!file.is_open())
        {
            throw std::runtime_error("Error: Could not create settings.ini file");
        }

        // Write header
        file << "# OpenThermo Settings File" << "\n";
        file << "# This file contains all available parameters with their default values" << "\n";
        file << "# Lines starting with # are comments" << "\n";
        file << "# parameter = value" << "\n";
        file << "# Values can be quoted: parameter = \"value with spaces\"" << "\n";
        file << "\n";

        // Electronic energy
        file << "# Electronic energy (a.u.) - overrides electronic energy from input file" << "\n";
        file << "# E = -76.384729" << "\n";
        file << "\n";

        // Temperature settings
        file << "# Temperature settings (K)" << "\n";
        file << "# Single temperature: T = 298.15" << "\n";
        file << "# Temperature scan: T = 200.0 400.0 25.0  (start, end, step)" << "\n";
        file << "T = 298.15" << "\n";
        file << "\n";

        // Pressure settings
        file << "# Pressure settings (atm)" << "\n";
        file << "# Single pressure: P = 1.0" << "\n";
        file << "# Pressure scan: P = 0.5 2.0 0.2  (start, end, step)" << "\n";
        file << "P = 1.0" << "\n";
        file << "\n";

        // Frequency scaling factors
        file << "# Vibrational frequency scaling factors" << "\n";
        file << "sclZPE = 1.0" << "\n";
        file << "sclheat = 1.0" << "\n";
        file << "sclS = 1.0" << "\n";
        file << "sclCV = 1.0" << "\n";
        file << "\n";

        // Low frequency treatment
        file << "# Low frequency treatment method" << "\n";
        file << "# 0 = Standard RRHO (harmonic approximation)" << "\n";
        file << "# 1 = Truhlar's QRRHO (frequency raising)" << "\n";
        file << "# 2 = Grimme's entropy interpolation (default)" << "\n";
        file << "# 3 = Minenkov's entropy + energy interpolation" << "\n";
        file << "lowvibmeth = 2" << "\n";
        file << "\n";

        // Low frequency parameters
        file << "# Parameters for low frequency treatments" << "\n";
        file << "ravib = 100.0" << "\n";
        file << "intpvib = 100.0" << "\n";
        file << "\n";

        // Calculation mode
        file << "# Calculation mode" << "\n";
        file << "# 0 = Gas phase (include translational/rotational)" << "\n";
        file << "# 1 = Condensed phase (remove translational/rotational)" << "\n";
        file << "ipmode = 0" << "\n";
        file << "\n";

        // Imaginary frequency treatment
        file << "# Imaginary frequency treatment" << "\n";
        file << "# Treat imaginary frequencies with |ν| < threshold as real" << "\n";
        file << "imagreal = 0.0" << "\n";
        file << "\n";

        // Concentration
        file << "# Concentration for phase correction" << "\n";
        file << "# conc = 1.0" << "\n";
        file << "\n";

        // Mass assignment
        file << "# Default atomic mass assignment" << "\n";
        file << "# 1 = Element average mass" << "\n";
        file << "# 2 = Most abundant isotope mass" << "\n";
        file << "# 3 = Masses from input file (default)" << "\n";
        file << "massmod = 3" << "\n";
        file << "\n";

        // Point group
        file << "# Force specific point group (auto-detect if not set)" << "\n";
        file << "# PGname = C2v" << "\n";
        file << "\n";

        // Output options
        file << "# Output options" << "\n";
        file << "# Print vibration contributions: 0=no, 1=yes, -1=to file" << "\n";
        file << "prtvib = 0" << "\n";
        file << "# Output .otm file: 0=no, 1=yes" << "\n";
        file << "outotm = 0" << "\n";
        file << "\n";

        // VASP energy selection
        file << "# VASP energy selection" << "\n";
        file << "# Select which energy to use from OUTCAR 'energy without entropy' line" << "\n";
        file << "# false/no/0 = energy  without entropy (default), true/yes/1 = energy(sigma->0)" << "\n";
        file << "extrape = false" << "\n";
        file << "\n";

        // Mass modifications section
        file << "# Mass modifications (optional)" << "\n";
        file << "# Uncomment and modify the following section to change specific atomic masses" << "\n";
        file << "# modmass" << "\n";
        file << "# 1 H 1.007825  # Atom 1: Hydrogen with specific mass" << "\n";
        file << "# 2 C 12.0      # Atom 2: Carbon-12 isotope" << "\n";
        file << "# 3 O 15.994915 # Atom 3: Oxygen-16 isotope" << "\n";
        file << "\n";

        file.close();
        std::cout << "Default settings.ini file created successfully!" << "\n";
        std::cout << "File location: " << filename << "\n";
        std::cout << "You can now edit this file to customize your default parameters." << "\n";
    }

}  // namespace util