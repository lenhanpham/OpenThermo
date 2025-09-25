/**
 * @file atommass.h
 * @brief Header for atomic mass and isotope data management
 * @author Le Nhan Pham
 * @date 2025
 *
 * This header file declares functions and data structures for handling
 * atomic mass calculations, isotope compositions, and element property
 * initialization in OpenThermo.
 */

#ifndef ATOMMASS_H
#define ATOMMASS_H

#include "chemsys.h"

/**
 * @brief Namespace containing atomic mass and isotope-related functions
 *
 * This namespace encapsulates all functions related to atomic mass data,
 * isotope compositions, and element property calculations.
 */
namespace atommass {

/**
 * @brief Initialize atomic mass data and isotope compositions
 *
 * This function initializes the global arrays containing isotope masses,
 * natural abundances, and average atomic masses for all supported elements.
 * The data is used throughout OpenThermo for mass-weighted calculations.
 *
 * @param sys Reference to the SystemData structure (currently unused but
 *            provided for future extensibility)
 *
 * @note This function must be called before any mass-related calculations
 * @note Initializes the external arrays: isomass, isowei, and elemass
 */
void initmass(SystemData& sys);

} // namespace atommass

#endif // ATOMMASS_H