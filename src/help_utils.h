/**
 * @file help_utils.h
 * @brief Help utility functions for OpenThermo program
 * @author Le Nhan Pham
 * @date 2025
 *
 * This file contains declarations for help-related utility functions
 * used throughout the OpenThermo molecular thermochemistry program.
 */

#ifndef HELP_UTILS_H
#define HELP_UTILS_H

#include <string>

namespace HelpUtils
{

    /**
     * @brief Print general help information
     * @param program_name Name of the program (default: "OpenThermo")
     */
    void print_help(const std::string& program_name = "OpenThermo");

    /**
     * @brief Print help for a specific command-line option
     * @param option The option to show help for (without the leading '-')
     * @param program_name Name of the program (default: "OpenThermo")
     */
    void print_option_help(const std::string& option, const std::string& program_name = "OpenThermo");

    /**
     * @brief Print input file format help
     */
    void print_input_help();

    /**
     * @brief Print output format help
     */
    void print_output_help();

    /**
     * @brief Print settings file help
     */
    void print_settings_help();

}  // namespace HelpUtils

#endif  // HELP_UTILS_H