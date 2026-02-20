Developer Guide
===============

This guide is for developers who want to contribute to the OpenThermo project. It covers the codebase structure, development setup, coding standards, and contribution guidelines.

Project Overview
================

Architecture
------------

OpenThermo is a C++17 application designed for high-performance molecular thermochemistry calculations. The codebase follows a modular architecture with clear separation of concerns:

.. code-block::
   ├── src/
   │   ├── main.cpp                           # Application entry point
   │   ├── atommass.cpp                       # Atomic mass handling
   │   ├── atommass.h
   │   ├── calc.cpp                           # Thermochemistry calculations
   │   ├── calc.h
   │   ├── chemsys.h                          # Chemical system data structures
   │   ├── help_utils.cpp                     # Help and utility functions
   │   ├── help_utils.h
   │   ├── loadfile.cpp                       # File loading and parsing
   │   ├── loadfile.h
   │   ├── omp_config.h                       # OpenMP configuration
   │   ├── symmetry.cpp                       # Symmetry analysis
   │   ├── symmetry.h
   │   ├── util.cpp                           # Utility functions
   │   ├── util.h
   │   └── version.h                          # Version information
   ├── tests/                                 # Regression test suite
   ├── docs/                                  # Documentation
   ├── resources/                             # Resources (logos, etc.)
   ├── CMakeLists.txt                         # CMake build configuration
   ├── Makefile                               # Make build system
   ├── Doxyfile                               # Doxygen configuration (if present)
   └── README.md                              # User documentation

Key Components
--------------

**Main Entry Point (main.cpp)**
   - Orchestrates the entire calculation workflow
   - Processes command-line arguments
   - Coordinates file loading, calculations, and output

**File Loading (loadfile.cpp/h)**
   - Parses quantum chemistry output files
   - Supports multiple formats (Gaussian, ORCA, GAMESS, NWChem, CP2K, VASP)
   - Extracts geometry, frequencies, and electronic energies

**Thermochemistry Calculations (calc.cpp/h)**
   - Implements statistical mechanics methods
   - Handles partition function calculations
   - Applies low-frequency treatment methods

**Symmetry Analysis (symmetry.cpp/h)**
   - Detects molecular point groups
   - Calculates rotational symmetry numbers
   - Handles linear molecule detection

**Atomic Mass Handling (atommass.cpp/h)**
   - Manages atomic mass assignments
   - Supports isotopic substitution
   - Handles custom mass modifications

**Utility Functions (util.cpp/h)**
   - Mathematical utilities
   - String processing
   - File I/O helpers

**OpenMP Configuration (omp_config.h)**
   - Thread management
   - HPC scheduler detection
   - Parallelization strategy selection

Key Design Principles
---------------------

**Modularity**
   - Each module has a single responsibility
   - Clear interfaces between components
   - Easy to test and maintain

**Performance**
   - Optimized C++17 implementation
   - Optional OpenMP parallelization
   - Efficient algorithms for large systems

**Accuracy**
   - Rigorous statistical mechanics implementation
   - Multiple low-frequency treatment methods
   - Comprehensive validation

**Usability**
   - Intuitive command-line interface
   - Extensive help system
   - Configuration file support

Development Setup
=================

Prerequisites
-------------

**Required Tools:**

- **C++ Compiler**: Clang 6.0+, GCC 7.0+, or Intel C++ Compiler 18.0+
- **Build System**: GNU Make 3.8.1+ or CMake 3.26.5+
- **Standard Library**: C++17 standard library
- **Git**: Version control system

**Optional Tools:**

- **Doxygen**: 1.8+ (API documentation generation)
- **Python**: 3.6+ (for test automation scripts)
- **OpenMP**: For parallel computation support
- **Valgrind**: Memory debugging
- **Clang-Tidy**: Code analysis

Getting the Source Code
-----------------------

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/lenhanpham/OpenThermo.git
   cd OpenThermo

   # Create a development branch
   git checkout -b feature/your-feature-name

Building for Development
------------------------

**Debug Build:**

.. code-block:: bash

   # Build with debug symbols and AddressSanitizer
   make debug

   # Or with CMake
   mkdir build && cd build
   cmake -DCMAKE_BUILD_TYPE=Debug ..
   make

**Release Build:**

.. code-block:: bash

   # Optimized release build
   make release

   # Or with CMake
   mkdir build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release ..
   make

**Development Build with OpenMP:**

.. code-block:: bash

   # Build with OpenMP support
   make OPENMP=1

   # Or with CMake
   mkdir build && cd build
   cmake -DENABLE_OPENMP=ON ..
   make

Testing
-------

**Running Tests:**

.. code-block:: bash

   # Build and run tests
   make test

   # Run with verbose output
   python3 tests/run_tests.py --verbose

   # Run specific test by ID
   python3 tests/run_tests.py --test gaussian_h2co_freq

**Test Coverage:**

The suite includes 27 tests organized by category:

- **Format coverage** (10 tests) — one per supported QC program at default settings
- **Option coverage** (9 tests) — different temperatures, pressures, low-frequency methods, scaling, concentration, imaginary frequency handling, and condensed phase mode
- **Output file tests** (3 tests) — verify ``.otm``, ``.vibcon``, and temperature scan (``.UHG``/``.SCq``) file generation
- **Error handling** (2 tests) — non-existent and invalid input files
- **Batch/ensemble** (1 test) — conformer ensemble with Boltzmann averaging

**Regenerating Reference Values:**

After modifying calculation logic:

.. code-block:: bash

   make test-generate

This runs OpenThermo on all test inputs and writes the results to ``tests/reference_values.json``.

Code Quality Tools
------------------

**Static Analysis:**

.. code-block:: bash

   # Run clang-tidy
   clang-tidy src/*.cpp -- -std=c++17 -Isrc

   # Run cppcheck
   cppcheck --enable=all --std=c++17 src/

**Code Formatting:**

.. code-block:: bash

   # Format code with clang-format
   find src/ -name "*.cpp" -o -name "*.h" | xargs clang-format -i

   # Check formatting
   find src/ -name "*.cpp" -o -name "*.h" | xargs clang-format --dry-run -Werror

Documentation
-------------

**Building Documentation:**

.. code-block:: bash

   # Install Sphinx
   pip install sphinx sphinx-rtd-theme furo

   # Build HTML documentation
   cd docs
   make html

   # View documentation
   firefox _build/html/index.html

**API Documentation:**

OpenThermo uses Doxygen-style comments for API documentation. To generate API docs:

.. code-block:: bash

   # Generate Doxygen documentation (if Doxyfile exists)
   doxygen Doxyfile

   # View API docs
   firefox doxygen/html/index.html

Coding Standards
================

Code Style
----------

**Naming Conventions:**

.. code-block:: cpp

   // Classes and structs
   class ChemicalSystem;
   struct CalculationParams;

   // Functions and methods
   void calculate_thermodynamics(const ChemicalSystem& system);
   double compute_partition_function(double temperature);

   // Variables
   int atom_count;
   std::string input_file;

   // Constants
   const double R = 8.3144648;           // Gas constant (J/mol/K)
   const double kb = 1.3806503e-23;      // Boltzmann constant (J/K)

   // Member variables (with m_ prefix)
   class MyClass {
   private:
       int m_atom_count;
       std::string m_input_file;
   };

**File Organization:**

- **Headers (.h)**: Class declarations, function prototypes, constants
- **Implementations (.cpp)**: Function definitions, implementation details
- **One class per file** when possible
- **Related functionality grouped** in modules

Documentation Standards
------------------------

**Doxygen Comments:**

OpenThermo uses Doxygen-style comments for all public APIs:

.. code-block:: cpp

   /**
    * @file calc.cpp
    * @brief Thermochemistry calculation functions
    * @author Le Nhan Pham
    * @date 2025
    *
    * This file contains functions for calculating thermodynamic properties
    * using statistical mechanics methods.
    */

   /**
    * @brief Calculate Gibbs free energy
    *
    * Computes the Gibbs free energy using partition functions and
    * statistical mechanics relationships.
    *
    * @param temperature Temperature in Kelvin
    * @param enthalpy Enthalpy in atomic units
    * @param entropy Entropy in J/mol/K
    * @return Gibbs free energy in atomic units
    *
    * @note This function assumes ideal gas behavior
    * @see calculate_enthalpy()
    * @see calculate_entropy()
    */
   double calculate_gibbs_free_energy(double temperature, double enthalpy, double entropy);

**Inline Comments:**

.. code-block:: cpp

   // Use comments for complex logic
   if (frequency < threshold) {
       // Apply low-frequency treatment method
       apply_low_frequency_treatment(frequency, method);
   }

   // Use TODO comments for future improvements
   // TODO: Optimize partition function calculation for large systems

Error Handling
-------------

**Exception Safety:**

.. code-block:: cpp

   try {
       // Operation that might fail
       load_molecular_data(filename);
   } catch (const std::invalid_argument& e) {
       // Handle invalid arguments
       std::cerr << "Invalid argument: " << e.what() << std::endl;
       return 1;
   } catch (const std::runtime_error& e) {
       // Handle runtime errors
       std::cerr << "Runtime error: " << e.what() << std::endl;
       return 2;
   } catch (const std::exception& e) {
       // Handle all other exceptions
       std::cerr << "Unexpected error: " << e.what() << std::endl;
       return 3;
   }

**Return Codes:**

.. code-block:: cpp

   /**
    * @return 0 on success
    * @return 1 on general error
    * @return 2 on invalid arguments
    * @return 3 on resource unavailable
    * @return 4 on operation interrupted
    */
   int process_molecular_data(const std::string& input_file);

Memory Management
-----------------

**RAII Pattern:**

.. code-block:: cpp

   class FileProcessor {
   public:
       FileProcessor(const std::string& filename)
           : m_file(filename) {
           if (!m_file.is_open()) {
               throw std::runtime_error("Failed to open file");
           }
       }

       ~FileProcessor() {
           // Automatic cleanup
           if (m_file.is_open()) {
               m_file.close();
           }
       }

   private:
       std::ifstream m_file;
   };

**Smart Pointers:**

.. code-block:: cpp

   // Use unique_ptr for exclusive ownership
   std::unique_ptr<ChemicalSystem> system = std::make_unique<ChemicalSystem>();

   // Use shared_ptr for shared ownership
   std::shared_ptr<CalculationParams> params = std::make_shared<CalculationParams>();

Thread Safety
-------------

**OpenMP Parallelization:**

OpenThermo uses OpenMP for parallel computation. When adding parallel code:

.. code-block:: cpp

   #ifdef _OPENMP
   #pragma omp parallel for
   for (size_t i = 0; i < frequencies.size(); ++i) {
       // Thread-safe operations
       results[i] = calculate_contribution(frequencies[i], temperature);
   }
   #endif

**Threading Guidelines:**

- Document thread safety guarantees
- Use appropriate synchronization primitives
- Avoid global mutable state
- Test concurrent access patterns

Contributing
============

Development Workflow
--------------------

**1. Choose an Issue:**

Check available issues on `GitHub <https://github.com/lenhanpham/OpenThermo/issues>`_.

**2. Create a Branch:**

.. code-block:: bash

   # Create and switch to feature branch
   git checkout -b feature/descriptive-name

   # Or for bug fixes
   git checkout -b bugfix/issue-number-description

**3. Make Changes:**

- Follow coding standards
- Add tests for new functionality
- Update documentation as needed
- Ensure all tests pass

**4. Test Your Changes:**

.. code-block:: bash

   # Build and test
   make debug
   make test

   # Run code quality checks
   clang-tidy src/*.cpp -- -std=c++17 -Isrc

**5. Commit Your Changes:**

.. code-block:: bash

   # Stage your changes
   git add .

   # Commit with descriptive message
   git commit -m "feat: add new feature description

   - What was changed
   - Why it was changed
   - How it was tested"

**6. Push and Create Pull Request:**

.. code-block:: bash

   # Push your branch
   git push origin feature/your-feature-name

   # Create pull request on GitHub

Pull Request Guidelines
-----------------------

**PR Title Format:**

.. code-block::

   type(scope): description

   Types: feat, fix, docs, style, refactor, test, chore

**PR Description Template:**

.. code-block::

   ## Description
   Brief description of the changes

   ## Type of Change
   - [ ] Bug fix
   - [ ] New feature
   - [ ] Breaking change
   - [ ] Documentation update

   ## Testing
   - [ ] Unit tests added/updated
   - [ ] Integration tests added/updated
   - [ ] Manual testing performed

   ## Checklist
   - [ ] Code follows style guidelines
   - [ ] Documentation updated
   - [ ] Tests pass
   - [ ] No breaking changes

Code Review Process
-------------------

**Review Checklist:**

- [ ] Code follows established patterns
- [ ] Appropriate error handling
- [ ] Thread safety considerations (if applicable)
- [ ] Performance implications
- [ ] Documentation updated
- [ ] Tests included
- [ ] No security vulnerabilities

**Review Comments:**

- Be constructive and specific
- Suggest improvements, don't just point out problems
- Reference coding standards when applicable
- Acknowledge good practices

Testing Guidelines
===================

Unit Testing
------------

**Test Structure:**

OpenThermo uses a regression test suite. When adding new features:

.. code-block:: python

   # tests/run_tests.py
   def test_new_feature():
       """Test new feature functionality."""
       result = run_openthermo("test_input.log", ["-new-option", "value"])
       assert result["gibbs_free_energy"] == expected_value

**Running Tests:**

.. code-block:: bash

   # Run all tests
   make test

   # Run specific test
   python3 tests/run_tests.py --test test_name

Integration Testing
-------------------

**End-to-End Tests:**

The regression test suite includes integration tests that verify:

- Complete calculation workflows
- Multiple input formats
- Various calculation methods
- Output file generation

Performance Testing
-------------------

**Benchmarking:**

.. code-block:: bash

   # Profile application
   valgrind --tool=callgrind OpenThermo molecule.log

   # Memory profiling
   valgrind --tool=massif OpenThermo molecule.log

Continuous Integration
======================

CI/CD Pipeline
--------------

**Automated Testing:**

OpenThermo uses GitHub Actions for CI/CD:

- **Build**: Compile on multiple platforms (Linux, macOS, Windows)
- **Test**: Run regression test suite
- **Sanitizers**: AddressSanitizer and ThreadSanitizer builds
- **OpenMP**: Test with and without OpenMP support

**GitHub Actions Workflow:**

Tests run automatically on every push via GitHub Actions:

- **Linux and macOS**: Build + regression tests (serial and OpenMP at 1 and 4 threads)
- **Windows**: MSYS2/MinGW build + regression tests (serial and OpenMP)
- **AddressSanitizer**: Separate build with ASan on Linux to detect memory errors
- **ThreadSanitizer**: OpenMP build with TSan on Linux to detect data races

Release Process
===============

Version Numbering
-----------------

**Semantic Versioning:**

.. code-block::

   MAJOR.MINOR.PATCH

   - MAJOR: Breaking changes
   - MINOR: New features (backward compatible)
   - PATCH: Bug fixes (backward compatible)

**Release Checklist:**

- [ ] Update version in ``src/version.h``
- [ ] Update CHANGELOG.md (if present)
- [ ] Update documentation
- [ ] Create release branch
- [ ] Run full test suite
- [ ] Create GitHub release
- [ ] Tag release

**Release Commands:**

.. code-block:: bash

   # Create release branch
   git checkout -b release/v0.001.5

   # Update version in version.h
   # Commit and tag
   git add src/version.h
   git commit -m "Release v0.001.5"
   git tag -a v0.001.5 -m "Release v0.001.5"

   # Push release
   git push origin release/v0.001.5
   git push origin v0.001.5

Support and Communication
=========================

**Communication Channels:**

- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: General questions and discussions
- **Pull Request Comments**: Code review discussions

**Getting Help:**

- Check existing issues and documentation first
- Use descriptive titles for issues
- Provide minimal reproducible examples
- Include system information and versions

**Community Guidelines:**

- Be respectful and constructive
- Help newcomers learn and contribute
- Follow the code of conduct
- Acknowledge contributions from others

This developer guide provides comprehensive information for contributing to the OpenThermo project. Following these guidelines ensures high-quality, maintainable code that benefits the entire community.
