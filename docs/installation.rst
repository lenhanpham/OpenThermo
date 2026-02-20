Installation Guide
==================

This guide covers the installation of OpenThermo on various platforms including Linux, macOS, and Windows.

System Requirements
-------------------

**Minimum Requirements:**

- **Operating System**: Linux, macOS 10.14+, or Windows 10+
- **Processor**: x86_64 architecture
- **Memory**: 512 MB RAM (2 GB recommended)
- **Storage**: 50 MB free space

**Recommended Requirements:**

- **Operating System**: Linux (Ubuntu 18.04+, CentOS 7+, Fedora 28+)
- **Compiler**: Clang 6.0+, GCC 7.0+, or Intel C++ Compiler 18.0+
- **Memory**: 2 GB+ RAM for large systems
- **Storage**: 100 MB+ for processing large output files

Supported Compilers
-------------------

OpenThermo supports multiple C++17 compilers:

- **Clang**: 6.0+ (recommended)
- **GCC**: 7.0+ (alternative)
- **Intel C++ Compiler**: 18.0+ (alternative)
- **MSVC**: 2017+ (Windows)

.. note::
   Clang is the primary compiler used for development and testing.

Build Dependencies
-------------------

**Required:**

- **GNU Make**: 3.8.1+ (build system, tested)
- **CMake**: 3.26.5+ (alternative build system, tested)
- **Standard Library**: C++17 standard library

**Optional:**

- **Doxygen**: 1.8+ (documentation generation)
- **Graphviz**: For documentation diagrams
- **Python**: 3.6+ (for test automation scripts)
- **OpenMP**: For parallel computation support

Installation Methods
====================

Quick Start Build
-----------------

**1. Clone or download the repository:**

.. code-block:: bash

   git clone https://github.com/lenhanpham/OpenThermo.git
   cd OpenThermo

**2. Build with Make (recommended):**

.. code-block:: bash

   make clean && make

**3. Alternative CMake build:**

.. code-block:: bash

   mkdir build && cd build
   cmake ..
   make

GNU Make Build System
---------------------

**Standard Build:**

.. code-block:: bash

   # Standard optimized build
   make

   # Build with OpenMP parallelization
   make OPENMP=1

   # Debug build with AddressSanitizer
   make debug

   # Release build with maximum optimization
   make release

   # Force specific compiler
   make CXX=g++      # GCC
   make CXX=clang++  # Clang
   make CXX=icpc     # Intel

   # Clean build artifacts
   make clean

**Build Targets:**

- ``all`` (default): Standard optimized build
- ``debug``: Includes debug symbols and AddressSanitizer
- ``release``: Maximum optimization for production use
- ``tsan``: ThreadSanitizer build for detecting data races (use with ``OPENMP=1``)
- ``clean``: Remove all build artifacts
- ``test``: Run the regression test suite
- ``test-generate``: Regenerate reference values from current build

CMake Build System
------------------

**Standard CMake Build:**

.. code-block:: bash

   # Create build directory
   mkdir build && cd build

   # Configure and build
   cmake ..
   make

   # Build with OpenMP parallelization
   cmake .. -DENABLE_OPENMP=ON
   make

**Advanced CMake Options:**

.. code-block:: bash

   # Build types
   cmake .. -DCMAKE_BUILD_TYPE=Release
   cmake .. -DCMAKE_BUILD_TYPE=Debug

   # Specify compiler
   cmake .. -DCMAKE_CXX_COMPILER=clang++

Windows Build (MSYS2/MinGW)
----------------------------

.. note::
   **For Windows Users, there is a binary package in the**
   `Releases <https://github.com/lenhanpham/OpenThermo/releases>`_.
   **Download this package and unzip it for further installation instructions.**

For those who would like to compile OpenThermo from source code on Windows:

**1. Install MSYS2:**

Download and install `MSYS2 <https://www.msys2.org/>`_ and open a MINGW64 terminal.

**2. Install dependencies:**

.. code-block:: bash

   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-clang mingw-w64-x86_64-cmake mingw-w64-x86_64-make

**3. Build:**

.. code-block:: bash

   make

The Makefile auto-detects Windows and adds ``-static`` linking so the binary runs without MinGW DLLs.

Post-Build Setup
================

Add OpenThermo to PATH
-----------------------

**Linux/macOS:**

Add the following line to your shell configuration file (``~/.bashrc``, ``~/.zshrc``, or equivalent),
replacing ``/path/to/OpenThermo`` with the actual path to the directory containing the ``OpenThermo``
binary (e.g. the repo root or ``build/``):

.. code-block:: bash

   export PATH=$PATH:/path/to/OpenThermo

Then reload your shell:

.. code-block:: bash

   source ~/.bashrc   # or source ~/.zshrc

**Windows:**

1. Open *Settings → System → Advanced system settings → Environment Variables*
2. Add the directory containing ``OpenThermo.exe`` to the ``Path`` variable
3. Open a new Command Prompt and test:

   .. code-block:: batch

      OpenThermo --help

Verification
------------

**1. Verify executable:**

.. code-block:: bash

   OpenThermo --help

**2. Create settings file (optional):**

.. code-block:: bash

   OpenThermo --create-config

This creates ``settings.ini`` with all available parameters and their default values.

**3. Test installation:**

.. code-block:: bash

   # Run a test calculation (if test files are available)
   OpenThermo test_molecule.log

Successful Compilation Output
------------------------------

Expected output during compilation:

.. code-block::

   Using compiler: g++
   g++ -std=c++17 -Wall -Wextra -O2 -c main.cpp -o build/main.o
   g++ -std=c++17 -Wall -Wextra -O2 -c atommass.cpp -o build/atommass.o
   ...
   g++ build/main.o ... -o OpenThermo -lrt -lstdc++fs

Environment Setup
-----------------

.. note::
   After adding OpenThermo to PATH (see :ref:`Post-Build Setup` above), you can call ``OpenThermo``
   directly from any directory without specifying its full path.

Troubleshooting
===============

Common Build Issues
-------------------

**Compiler Not Found:**

.. code-block:: bash

   # Check available compilers
   which g++ clang++ icpc

   # Install missing compiler
   sudo apt install g++-10  # Ubuntu/Debian
   sudo yum install gcc-c++  # CentOS/RHEL

**C++17 Support Missing:**

.. code-block:: bash

   # Check compiler version
   g++ --version

   # Upgrade compiler if needed
   sudo apt install g++-10
   sudo update-alternatives --config g++

**Library Issues:**

.. code-block:: bash

   # Check for required libraries
   ldconfig -p | grep stdc++

   # Rebuild if libraries are missing
   make clean && make

**Permission Issues:**

.. code-block:: bash

   # Fix permissions
   chmod +x build/OpenThermo

**OpenMP Issues:**

If OpenMP is not found during build:

.. code-block:: bash

   # Install OpenMP development libraries
   sudo apt install libomp-dev  # Ubuntu/Debian
   sudo yum install libgomp    # CentOS/RHEL

Runtime Issues
--------------

**File Not Found:**

.. code-block:: bash

   # Check file path and spelling
   ls -la molecule.log

   # Ensure file is in current directory or provide full path
   OpenThermo /full/path/to/molecule.log

**Unrecognized Program Format:**

Ensure the input file contains frequency analysis output from a supported quantum chemistry program.

**Invalid Arguments:**

.. code-block:: bash

   # Check argument syntax
   OpenThermo --help-T

   # Use proper format for temperature scan
   OpenThermo molecule.log -T 200 400 25

Performance Optimization
========================

**Compiler Selection:**

- **Clang**: Recommended for best compatibility
- **GCC**: Good general performance
- **Intel**: Best performance on Intel systems

**Build Optimization:**

.. code-block:: bash

   # Release build with full optimization
   make release -j $(nproc)

   # Use all available cores for compilation
   make -j $(nproc)

**OpenMP Parallelization:**

.. code-block:: bash

   # Build with OpenMP support
   make OPENMP=1

   # Or with CMake
   cmake .. -DENABLE_OPENMP=ON
   make

**Runtime Optimization:**

OpenThermo automatically uses half of detected physical CPU cores by default. For dedicated compute nodes:

.. code-block:: bash

   # Use explicit thread count
   OpenThermo molecule.log -omp-threads 8

Testing
=======

**Run Test Suite:**

.. code-block:: bash

   # Run the full test suite
   make test

   # Run with verbose output
   python3 tests/run_tests.py --verbose

   # Run a single test by ID
   python3 tests/run_tests.py --test gaussian_h2co_freq

**Regenerate Reference Values:**

After modifying calculation logic:

.. code-block:: bash

   make test-generate

This runs OpenThermo on all test inputs and writes the results to ``tests/reference_values.json``.

Uninstallation
==============

**Source Installation:**

.. code-block:: bash

   # Remove build artifacts
   make clean

   # Remove binary
   rm build/OpenThermo

   # Remove settings file (if created)
   rm settings.ini

Getting Help
============

If you encounter issues:

1. Check the :doc:`usage` guide for proper usage
2. Use ``OpenThermo --help`` for command-line help
3. Check system requirements and compiler compatibility
4. Report issues on GitHub with system information

.. code-block:: bash

   # Get version information
   OpenThermo --help

   # Check compiler information
   make compiler-info
