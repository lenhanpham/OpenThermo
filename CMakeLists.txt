cmake_minimum_required(VERSION 3.10)
project(OpenThermo VERSION 0.000.1 LANGUAGES CXX)

# Prevent in-source builds
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.")
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set compiler flags
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Od /Zi")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
endif()

# Platform-specific definitions
if(WIN32)
    add_definitions(-D_WIN32)
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Source files
set(SOURCES
    src/main.cpp
    src/atommass.cpp
    src/calc.cpp
    src/loadfile.cpp
    src/symmetry.cpp
    src/util.cpp
    src/help_utils.cpp
)

# Header files
set(HEADERS
    src/defvar.h
    src/atommass.h
    src/calc.h
    src/loadfile.h
    src/symmetry.h
    src/util.h
    src/help_utils.h
)

# Add Windows resource file if building on Windows
if(WIN32)
    # Enable RC language BEFORE checking files
    enable_language(RC)
    
    set(RESOURCE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/resources/icon.rc")
    set(ICON_FILE "${CMAKE_CURRENT_SOURCE_DIR}/resources/icon.ico")
    
    if(NOT EXISTS "${RESOURCE_FILE}")
        message(FATAL_ERROR "Resource file ${RESOURCE_FILE} not found")
    endif()
    if(NOT EXISTS "${ICON_FILE}")
        message(FATAL_ERROR "Icon file ${ICON_FILE} not found")
    endif()
    
    message(STATUS "Adding Windows resource file: ${RESOURCE_FILE}")
    
    # Set resource compiler for MinGW/Clang
    if(CMAKE_RC_COMPILER STREQUAL "windres" OR MINGW)
        set(CMAKE_RC_COMPILE_OBJECT "<CMAKE_RC_COMPILER> <FLAGS> -O coff <DEFINES> -i <SOURCE> -o <OBJECT>")
    endif()
    
    # Add resource file to sources
    list(APPEND SOURCES "${RESOURCE_FILE}")
endif()

# Create executable
add_executable(OpenThermo ${SOURCES} ${HEADERS})

# Set output directory - ensures binary is created in build directory, not source directory
set_target_properties(OpenThermo PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
message(STATUS "OpenThermo binary will be created in: ${CMAKE_BINARY_DIR}")

# Installation rules
install(TARGETS OpenThermo
    RUNTIME DESTINATION bin
)

# Optional: Create a test target
enable_testing()
add_test(NAME BasicTest COMMAND OpenThermo --help)