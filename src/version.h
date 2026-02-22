/**
 * @file version.h
 * @brief Centralized version, release date, and author information for OpenThermo
 * @author Le Nhan Pham
 * @date 2026
 *
 * All version-related constants are defined here so that updating
 * the version number or author list only needs to happen in one place.
 */

#ifndef VERSION_H
#define VERSION_H

#include <string>

namespace OpenThermo
{
    /// Program version string (e.g. "0.001.3")
    inline constexpr const char* VERSION       = "0.001.6";

    /// Release year shown in the banner
    inline constexpr const char* RELEASE_DATE  = "2026";

    /// Short developer list for the banner
    inline constexpr const char* AUTHORS       = "Le Nhan Pham, Romain Despoullains";

    /// Project URL
    inline constexpr const char* PROJECT_URL   = "https://github.com/lenhanpham/openthermo";

}  // namespace OpenThermo

#endif  // VERSION_H
