<p align="center">
  <a href="https://github.com/lenhanpham/OpenThermo">
    <picture>
      <img src="resources/ot-logo.svg" alt="OpenThermo">
    </picture>
  </a>
</p>
<p align="center">
  <a href="https://en.cppreference.com/w/cpp/17">
    <img src="https://img.shields.io/badge/C%2B%2B-17-blue.svg" alt="C++17">
  </a>
  <a href="./LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
  </a>
  <a href="https://github.com/lenhanpham/OpenThermo/actions/workflows/ot.yml">
    <img src="https://github.com/lenhanpham/OpenThermo/actions/workflows/ot.yml/badge.svg" alt="Build Status">
  </a>
</p>

A comprehensive C++ program for calculating molecular thermochemistry properties from quantum chemistry output files (Gaussian, Orca, GAMESS, NWCHEM, CP2K, VASP). OpenThermo implements state-of-the-art methods for computing thermodynamic quantities including Gibbs free energy, enthalpy, entropy, and heat capacities using statistical mechanics.

# Important Note

**The project is still in the early stage and not fully tested. Therefore, errors and inaccuracy may happen. Users are suggested to check calculated data against original data from outputs of corresponding quantum chemical programs (with default temperature, conentration, and pressure)**

A graphical user interface (GUI) version is under development and can be found here: 
[https://github.com/lenhanpham/OpenThermoGUI](https://github.com/lenhanpham/OpenThermoGUI) 

## Table of Contents

1. [Project Overview](#project-overview)
2. [Key Features](#key-features)
3. [System Requirements](#system-requirements)
4. [Installation](#installation)
5. [Configuration](#configuration)
6. [Usage](#usage)
7. [Input Formats](#input-formats)
8. [Output Formats](#output-formats)
9. [Calculation Methods](#calculation-methods)
10. [Advanced Features](#advanced-features)
11. [Examples](#examples)
12. [Troubleshooting](#troubleshooting)
13. [Contributing](#contributing)
14. [License](#license)

## Project Overview

**OpenThermo** is a high-performance C++17 program designed for molecular thermochemistry calculations. It processes quantum chemistry output files to compute thermodynamic properties using rigorous statistical mechanics methods. The program supports multiple quantum chemistry packages and implements both standard and advanced treatment methods for accurate thermochemical analysis.

### Purpose

OpenThermo addresses the need for accurate thermochemical data in computational chemistry without re-running quantum chemical programs by providing:

- **Accurate thermodynamic calculations** using statistical mechanics
- **Multiple treatment methods** for low-frequency vibrational modes
- **Comprehensive analysis** of molecular systems
- **Batch processing capabilities** for high-throughput workflows
- **Extensible architecture** for future enhancements

### Architecture

The codebase is organized into modular components located in the `src/` directory:

```
src/
├── atommass.cpp
├── atommass.h
├── calc.cpp
├── calc.h
├── chemsys.h
├── help_utils.cpp
├── help_utils.h
├── loadfile.cpp
├── loadfile.h
├── main.cpp
├── omp_config.h
├── symmetry.cpp
├── symmetry.h
├── util.cpp
├── util.h
└── version.h
```

### Development Status

- **Version**: 0.001.4
- **Language**: C++17 with clang++ compiler
- **Build System**: GNU Make with CMake support
- **Testing**: 27-test regression suite with CI on Linux, macOS, and Windows (+ AddressSanitizer, ThreadSanitizer)
- **Documentation**: Doxygen-generated API documentation (not ready)

## Key Features

### Core Functionality

- **Multi-format Support**: Gaussian, ORCA, GAMESS-US, NWChem, CP2K, VASP
- **Advanced Thermochemistry**: Standard RRHO and quasi-RRHO treatments for low-frequency modes
- **Statistical Mechanics**: Rigorous implementation of partition functions and thermodynamic properties
- **Symmetry Analysis**: Automatic point group detection and rotational symmetry number calculation

### Calculation Methods

- **Standard RRHO**: Rigid-rotor harmonic oscillator approximation
- **Quasi-RRHO Treatments**:
  - Truhlar's frequency raising method
  - Grimme's entropy interpolation
  - Minenkov's energy + entropy interpolation
  - Head-Gordon's energy interpolation (with optional entropy) [Experimental]
- **Electronic Contributions**: Multi-level electronic structure support
- **Concentration Corrections**: Solution-phase Gibbs energy adjustments

### Flexibility & Usability

- **Temperature/Pressure Scanning**: Automated parameter sweeps
- **Batch Processing**: Multi-file analysis with ensemble averaging
- **Mass Modifications**: Isotopic substitution and custom atomic masses
- **Configurable Scaling**: Frequency scaling factors for different properties
- **Multiple Output Formats**: Console, text files, and native OpenThermo format

### Technical Excellence

- **High Performance**: Optimized C++17 implementation
- **Memory Efficient**: Minimal memory footprint for large systems
- **Error Handling**: Comprehensive validation and graceful error recovery
- **Extensible Design**: Modular architecture for easy enhancement
- **Cross-Platform**: Linux, macOS, and Windows support

## System Requirements

### Minimum Requirements

- **Operating System**: Linux, macOS 10.14+, or Windows 10+
- **Processor**: x86_64 architecture
- **Memory**: 512 MB RAM (2 GB recommended)
- **Storage**: 50 MB free space

### Compiler Requirements

- **Primary**: Clang 6.0+ or GCC 7.0+
- **Alternative**: Intel C++ Compiler 18.0+
- **Standard**: C++17 compliant compiler

### Build Dependencies

- **GNU Make**: 3.8.1+ (build system, tested)
- **CMake**: 3.26.5+ (alternative build system, tested)
- **Standard Library**: C++17 standard library

### Optional Dependencies

- **Doxygen**: 1.8+ (documentation generation)
- **Graphviz**: For documentation diagrams
- **Python**: 3.6+ (for test automation scripts)

### Supported Platforms

- **Linux**: Ubuntu 18.04+, CentOS 7+, Fedora 28+
- **macOS**: 10.14+ with Xcode Command Line Tools
- **Windows**: 10+ with MSYS2/MinGW or Visual Studio 2017+

## Installation

### Quick Start Build

1. **Clone or download the repository:**
   
   ```bash
   git clone https://github.com/lenhanpham/OpenThermo.git
   cd OpenThermo
   ```

2. **Build with Make (recommended):**
   
   ```bash
   make clean && make
   ```

3. **Alternative CMake build:**
   
   ```bash
   mkdir build && cd build
   cmake ..
   make
   ```

### Build Options

#### GNU Make Build System

```bash
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
```

#### CMake Build System

```bash
# Create build directory
mkdir build && cd build

# Configure and build
cmake ..
make

# Build with OpenMP parallelization
cmake .. -DENABLE_OPENMP=ON
make

# Build types
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake .. -DCMAKE_BUILD_TYPE=Debug

# Specify compiler
cmake .. -DCMAKE_CXX_COMPILER=clang++
```

#### Windows Build (MSYS2/MinGW)

1. Install [MSYS2](https://www.msys2.org/) and open a MINGW64 terminal
2. Install dependencies:
   ```bash
   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-clang mingw-w64-x86_64-cmake mingw-w64-x86_64-make
   ```
3. Build:
   ```bash
   make
   ```
   The Makefile auto-detects Windows and adds `-static` linking so the binary runs without MinGW DLLs.

### Post-Build Setup

1. **Verify executable:**
   
   ```bash
   ./build/OpenThermo --help
   ```

2. **Create settings file (optional):**
   
   ```bash
   ./build/OpenThermo --create-config
   ```

3. **Compiler options:**
   
   ```bash
   make CXX=g++     # Force GCC compilation
   make CXX=icpc    # Force Intel compiler
   ```

4. **Clean build:**
   
   ```bash
   make clean
   ```

### Build Targets

- `all` (default): Standard optimized build
- `debug`: Includes debug symbols and AddressSanitizer
- `release`: Maximum optimization for production use
- `tsan`: ThreadSanitizer build for detecting data races (use with `OPENMP=1`)
- `clean`: Remove all build artifacts
- `test`: Run the regression test suite
- `test-generate`: Regenerate reference values from current build

### Successful Compilation Output

```
Using compiler: g++
g++ -std=c++17 -Wall -Wextra -O2 -c main.cpp -o build/main.o
g++ -std=c++17 -Wall -Wextra -O2 -c atommass.cpp -o build/atommass.o
...
g++ build/main.o ... -o OpenThermo -lrt -lstdc++fs
```

## Configuration

### Configuration Files

OpenThermo supports configuration through:

1. **Command-line arguments** (highest priority)
2. **Local settings file**: `./settings.ini`
3. **Environment settings**: `$OPENTHERMOPATH/settings.ini`
4. **Program defaults** (lowest priority)

### Creating Settings File

Generate a default settings file:

```bash
./build/OpenThermo --create-config
```

This creates `settings.ini` with all available parameters and their default values.

### Settings File Format

```ini
# OpenThermo Settings File
# Lines starting with # are comments

# Electronic energy override
E = -76.384729

# Temperature settings (K)
T = 298.15
# For temperature scan: T = 200.0 400.0 25.0

# Pressure settings (atm)
P = 1.0
# For pressure scan: P = 0.5 2.0 0.2

# Frequency scaling factors
sclZPE = 1.0      # Zero-point energy scaling
sclheat = 1.0     # Thermal energy scaling
sclS = 1.0        # Entropy scaling
sclCV = 1.0       # Heat capacity scaling

# Low frequency treatment
lowvibmeth = Grimme  # 0/Harmonic=RRHO, 1/Truhlar, 2/Grimme, 3/Minenkov, 4/HeadGordon
ravib = 100.0     # Raising threshold for Truhlar method
intpvib = 100.0   # Interpolation threshold for Grimme/Minenkov/HeadGordon
hg_entropy = true # Enable entropy interpolation for Head-Gordon method
# bav = qchem    # Bav for HeadGordon only: grimme or qchem (default: qchem). Ignored for Grimme/Minenkov.

# Calculation options
ipmode = 0         # 0=gas phase, 1=condensed phase
imagreal = 0.0    # Imaginary frequency threshold
massmod = 3       # Mass assignment: 1=average, 2=abundant, 3=file
PGname = "?"     # Point group (auto-detect if "?")

# Output options
prtlevel = 1      # Verbosity: 0=minimal, 1=default, 2=verbose, 3=full
prtvib = 0        # Vibration contributions: 0=no, 1=screen, -1=file
outotm = 0        # Output .otm file: 0=no, 1=yes

# Concentration (for solution phase)
conc = 1.0

# VASP energy selection
# false/no/0 = energy  without entropy (default), true/yes/1 = energy(sigma->0)
extrape = false

# OpenMP threading (command-line or settings.ini)
# Default: half of physical CPU cores or scheduler-allocated CPUs
# On HPC: SLURM_CPUS_PER_TASK / PBS_NP / NSLOTS are auto-detected
# omp-threads = 4
# Override with CLI: -omp-threads N (takes precedence)
# Strategy (outer T/P vs inner vibrational) is auto-selected

# Mass modifications (optional section)
# modmass
# 1 H 1.007825  # Atom 1: Hydrogen with specific mass
# 2 C 12.0      # Atom 2: Carbon-12 isotope
```

### Parameters

| Parameter    | Meaning                                                          | Default Value |
| ------------ | ---------------------------------------------------------------- | ------------- |
| `conc`       | Concentration string for solution phase Gibbs energy corrections | `0`           |
| `prtlevel`   | Output verbosity level (0=minimal, 1=default, 2=verbose, 3=full) | `1`           |
| `prtvib`     | Print vibration contributions (0=no, 1=screen, -1=file)          | `0`           |
| `lowvibmeth` | Low frequency treatment method                                   | `2`           |
| `massmod`    | Mass assignment mode (1=average, 2=abundant, 3=file)             | `3`           |
| `outotm`     | Output .otm file flag (0=no, 1=yes)                              | `0`           |
| `ipmode`     | Calculation mode (0=gas phase, 1=condensed phase)                | `0`           |
| `T`          | Temperature in Kelvin                                            | `298.15`      |
| `P`          | Pressure in atmospheres                                          | `1.0`         |
| `sclZPE`     | ZPE scaling factor                                               | `1.0`         |
| `sclheat`    | Thermal energy scaling factor                                    | `1.0`         |
| `sclS`       | Entropy scaling factor                                           | `1.0`         |
| `sclCV`      | Heat capacity scaling factor                                     | `1.0`         |
| `ravib`      | Raising threshold for Truhlar method (cm⁻¹)                      | `100.0`       |
| `intpvib`    | Interpolation threshold for Grimme/Minenkov/HeadGordon (cm⁻¹)    | `100.0`       |
| `hg_entropy` | Entropy interpolation for Head-Gordon method                     | `true`        |
| `bav`        | Bav preset for HeadGordon free-rotor entropy (grimme/qchem)          | `qchem`       |
| `imagreal`   | Imaginary frequency threshold (cm⁻¹)                             | `0.0`         |
| `Eexter`     | External electronic energy override (a.u.)                       | `0.0`         |
| `extrape`    | VASP electronic energy selection                                 | `false`       |
| `PGname`     | Point group name ("?" for auto-detect)                           | `"?"`         |
| `omp-threads`| OpenMP thread count (clamped to effective cores)                 | Half of effective cores |

**Details for `extrape`:**

- `false/no/0` → energy without entropy  
- `true/yes/1` → energy (σ→0)

**Details for `lowvibmeth`:**

- 0/Harmonic = Rigid-Rotor Harmonic Oscillator approximation (RRHO)
- 1/Truhlar =  Truhlar's QRRHO method
- 2/Grimme = Grimme's method
- 3/Minenkov = Minenkov's method
- 4/HeadGordon = Head-Gordon's energy interpolation (with optional entropy) [Experimental]

### Parameter Precedence

1. **Command-line arguments** (highest priority)
2. **Current directory** `settings.ini`
3. **Environment path** `settings.ini`
4. **Program defaults** (lowest priority) 

### Environment Variables

- **`OPENTHERMOPATH`**: Directory containing `settings.ini`
- **`OPENTHERMO_OUTPUT`**: Default output directory

### Configuration Priority

Settings are applied in this order:

1. Command-line arguments override everything
2. Current directory `settings.ini`
3. Environment path `settings.ini`
4. Built-in program defaults

### Command Line Options

### Input and Output Options

#### `-E <value>`

- **Description**: Override electronic energy
- **Units**: Atomic units (a.u.)
- **Default**: Use value from input file
- **Example**: `-E -76.384729`

#### `-T <T>` or `-T <T1 T2 step>`

- **Description**: Temperature specification
- **Units**: Kelvin (K)
- **Default**: 298.15 K
- **Examples**:
  - `-T 300` (single temperature)
  - `-T 200 400 25` (scan from 200K to 400K in 25K steps)

#### `-P <P>` or `-P <P1 P2 step>`

- **Description**: Pressure specification
- **Units**: Atmospheres (atm)
- **Default**: 1.0 atm
- **Examples**:
  - `-P 2.0` (single pressure)
  - `-P 0.5 2.0 0.2` (scan from 0.5 to 2.0 atm in 0.2 atm steps)

### Thermochemistry Options

#### `-sclZPE <factor>`

- **Description**: Scale factor for ZPE frequencies
- **Default**: 1.0
- **Example**: `-sclZPE 0.98`

#### `-sclheat <factor>`

- **Description**: Scale factor for thermal energy frequencies
- **Default**: 1.0
- **Example**: `-sclheat 0.99`

#### `-sclS <factor>`

- **Description**: Scale factor for entropy frequencies
- **Default**: 1.0
- **Example**: `-sclS 0.99`

#### `-sclCV <factor>`

- **Description**: Scale factor for heat capacity frequencies
- **Default**: 1.0
- **Example**: `-sclCV 0.99`

### Low Frequency Treatment Options

#### `-lowvibmeth <mode>`

- **Description**: Low frequency treatment method
- **Values**:
  - `0` or `Harmonic`: Standard RRHO (harmonic approximation)
  - `1` or `Truhlar`: Truhlar's QRRHO (frequency raising)
  - `2` or `Grimme`: Grimme's entropy interpolation
  - `3` or `Minenkov`: Minenkov's entropy + energy interpolation
  - `4` or `HeadGordon`: Head-Gordon's energy interpolation (+ optional entropy) [Experimental]
- **Default**: Grimme
- **Example**: `-lowvibmeth 1` or `-lowvibmeth Truhlar`

#### `-ravib <value>`

- **Description**: Raising value for low frequencies (Truhlar's method)
- **Units**: cm⁻¹
- **Default**: 100.0
- **Example**: `-ravib 50.0`

#### `-intpvib <value>`

- **Description**: Interpolation frequency threshold for Grimme, Minenkov, and Head-Gordon methods
- **Units**: cm⁻¹
- **Default**: 100.0
- **Example**: `-intpvib 50.0`

#### `-hg_entropy <bool>`

- **Description**: Enable or disable entropy interpolation for the Head-Gordon method
- **Values**: `true`/`false` or `1`/`0`
- **Default**: true
- **Example**: `-hg_entropy false`

#### `-bav <preset>`

- **Description**: Average moment of inertia (Bav) used in the free-rotor entropy term. **Only applicable to the HeadGordon method** (`lowvibmeth=4`). Grimme and Minenkov methods always use `grimme` (1×10⁻⁴⁴ kg m²); this option is ignored for them.
- **Values**:
  - `grimme`: I_av = 1×10⁻⁴⁴ kg m² (μ'_av from Grimme 2012, used by ORCA/xtb/GoodVibes/Shermo)
  - `qchem`: I_av = 2.79928×10⁻⁴⁶ kg m² (B_av = 1 cm⁻¹, as specified in the Q-Chem manual)
- **Default**: `qchem` (for HeadGordon)
- **Example**: `-bav grimme`

#### `-ipmode <mode>`

- **Description**: Calculation mode
- **Values**:
  - `0`: Gas phase (include translational/rotational)
  - `1`: Condensed phase (remove translational/rotational)
- **Default**: 0
- **Example**: `-ipmode 1`

#### `-imagreal <value>`

- **Description**: Treat imaginary frequencies as real if |ν| < value
- **Units**: cm⁻¹
- **Default**: 0.0 (no treatment)
- **Example**: `-imagreal 50.0`

### Mass and Symmetry Options

#### `-massmod <type>`

- **Description**: Default atomic mass assignment
- **Values**:
  - `1`: Element average mass
  - `2`: Most abundant isotope mass
  - `3`: Masses from input file (default)
- **Default**: 3
- **Example**: `-massmod 2`

#### `-PGname <name>`

- **Description**: Force specific point group
- **Default**: Auto-detect
- **Example**: `-PGname C2v`

#### `-conc <string>`

- **Description**: Concentration specification
- **Default**: None
- **Example**: `-conc "2.5"`

### Output Control Options

#### `-prtlevel <level>`

- **Description**: Controls the amount of information printed to the console
- **Values**:
  - `0`: Minimal — banner + final thermodynamic data only
  - `1`: Default — parameter summary, compact system info (atom count, point group, frequency range), and final data
  - `2`: Verbose — full system data (per-atom masses, full frequency list, moments of inertia, rotational constants) and component breakdown (Translation/Rotation/Vibration/Electronic sections)
  - `3`: Full — everything in level 2 plus per-mode vibrational detail tables (automatically enables `-prtvib 1` unless explicitly overridden)
- **Default**: 1
- **Examples**:
  - `-prtlevel 0` (minimal output for scripting)
  - `-prtlevel 2` (full output matching pre-0.001.1 behavior)
- **Settings file**: `prtlevel = 1`

#### `-prtvib <mode>`

- **Description**: Print vibration contributions
- **Values**:
  - `0`: No (default)
  - `1`: Yes, to screen
  - `-1`: Yes, to `basename.vibcon` file
- **Default**: 0
- **Example**: `-prtvib 1`

#### `-outotm <mode>`

- **Description**: Output OpenThermo format file
- **Values**:
  - `0`: No (default)
  - `1`: Yes
- **Default**: 0
- **Example**: `-outotm 1`

### Special Options

#### `-noset`

- **Description**: Skip loading settings from `settings.ini`
- **Example**: `-noset`

### OpenMP / Performance Options

#### `-omp-threads <N>`

- **Description**: Set the number of OpenMP threads for parallel computation
- **Default**: Half of the detected physical CPU cores (minimum 1), or half of scheduler-allocated CPUs on HPC
- **HPC scheduler support**: Automatically detects allocated CPUs from `SLURM_CPUS_PER_TASK`, `PBS_NUM_PPN`, `PBS_NP`, `NSLOTS`, `LSB_DJOB_NUMPROC`
- **Behavior**:
  - If not specified: uses half of effective cores (scheduler allocation or physical)
  - If `N` ≤ effective cores: uses `N` threads
  - If `N` > effective cores: clamps to default (half) with a warning
- **Strategy auto-selection**: OpenThermo automatically selects the best parallelization strategy:
  - **Outer** (T/P scan parallel): when many T/P scan points are available (≥ nthreads)
  - **Inner** (vibrational parallel): when few T/P points but many frequencies (>50)
- **Precedence**: CLI `-omp-threads` > `settings.ini` `omp-threads` > auto-detect
- **Examples**:
  - `-omp-threads 4` (use 4 threads, if ≤ effective cores)
  - `-omp-threads 2` (use 2 threads)
- **Settings file**: `omp-threads = 4` (in settings.ini)

### Help Options

#### `--help`

- **Description**: Show general help
- **Example**: `./OpenThermo --help`

#### `--help-<option>`

- **Description**: Show help for specific option
- **Examples**:
  - `--help-T` (temperature help)
  - `--help-lowvibmeth` (low frequency help)
  - `--help-input` (input formats help)

## Input File Formats

### Supported Formats

#### 1. OpenThermo Format (.otm)

**Native format containing all molecular data:**

```text
*E  //Electronic energy (a.u.)
-76.384729

*wavenum  //Wavenumbers (cm^-1)
1234.5
2345.6
...

*atoms  //Name, mass (amu), X, Y, Z (Angstrom)
C   12.000000   0.000000   0.000000   0.000000
H    1.007825   0.000000   0.000000   1.089000
...

*elevel  //Energy (eV) and degeneracy
0.0 1
1.5 3
```

#### 2. Quantum Chemistry Output Files

##### Gaussian (.log, .out)

- **Requirements**: Frequency analysis output
- **Extracts**: Geometry, frequencies, electronic energy
- **Features**: Automatic format detection

##### ORCA (.out)

- **Requirements**: Frequency calculation results
- **Features**: Supports various ORCA output formats

##### GAMESS-US (.log)

- **Requirements**: Hessian/frequency analysis
- **Features**: Standard GAMESS output parsing

##### NWChem (.out)

- **Requirements**: Frequency analysis output
- **Features**: Comprehensive NWChem support

##### CP2K (.out)

- **Requirements**: Vibrational analysis output
- **Features**: Supports molecular and periodic systems
- **Note**: For condensed phase systems (ipmode=1): contributions of translation and rotation are ignored

##### VASP (OUTCAR)

- **Requirements**: Vibrational analysis output in OUTCAR, sysmtem information in CONTCAR
- **Features**: Supports molecular and periodic systems
- **Note**: For condensed phase systems (ipmode=1): contributions of translation and rotation are ignored

#### 3. List Files (.txt)

**Batch processing of multiple files:**

```text
molecule-a.log
molecule-b.out
/path/to/molecule-c.otm
```

### File Detection

- **Automatic**: Program detects format from file content
- **No manual specification** required
- **Priority**: OpenThermo → specific program signatures

## Output Files and Formats

### Console Output

```
OpenMP threads: 24 (default: half of 48 physical cores). Use -omp-threads N to override.

                   --- Summary of Current Parameters ---

Running parameters:
 Print level: 1 (0=minimal, 1=default, 2=verbose, 3=full)
Printing individual contribution of vibration modes: No
 Temperature:          298.150 K
 Pressure:             1.000 atm
 Scaling factor of vibrational frequencies for ZPE:         1.0000
 Scaling factor of vibrational frequencies for U(T)-U(0):   1.0000
 Scaling factor of vibrational frequencies for S(T):        1.0000
 Scaling factor of vibrational frequencies for CV:          1.0000
 Low frequencies treatment: Grimme's interpolation for entropy
 Vibrational frequency threshold used in the interpolation is 100.00 cm^-1
                      -------- End of Summary --------

OpenThermo started to process BIH-conformers-1.log at Tue Feb 17 14:10:27 2026

                      -------- Chemical System Data -------
                      -------------------------------------
 Electronic energy:      -690.56452500 a.u.
 Electronic energy level 1     E =     0.000000 eV     Degeneracy =   1
 Atoms: 33 (16 H, 15 C, 2 N)  Total mass: 224.131420 amu
 Point group: Cs     Rotational symmetry number:   1
 Frequencies: 93 (range: 25.1 -- 3233.0 cm^-1)


                       ----------------------------
                       -------- Final data --------
                       ----------------------------
 Total q(V=0):        8.578987e+42
 Total q(bot):        8.806492e-87
 Total q(V=0)/NA:     1.424574e+19
 Total q(bot)/NA:    1.462352e-110
 Total CV:     237.069 J/mol/K      56.661 cal/mol/K
 Total CP:     245.383 J/mol/K      58.648 cal/mol/K
 Total S:      484.246 J/mol/K     115.738 cal/mol/K    -TS:   -34.507 kcal/mol
 Zero point energy (ZPE):    736.269 kJ/mol    175.973 kcal/mol   0.280430 a.u.
 Thermal correction to U:    773.673 kJ/mol    184.912 kcal/mol   0.294677 a.u.
 Thermal correction to H: 776.152210 kJ/mol 185.504830 kcal/mol   0.295621 a.u.
 Thermal correction to G: 631.774188 kJ/mol 150.997655 kcal/mol   0.240630 a.u.
 Electronic energy:       -690.5645250 a.u.
 Sum of electronic energy and ZPE, namely U/H/G at 0 K:       -690.2840949 a.u.
 Sum of electronic energy and thermal correction to U:        -690.2698485 a.u.
 Sum of electronic energy and thermal correction to H:        -690.2689043 a.u.
 Sum of electronic energy and thermal correction to G:        -690.3238950 a.u.

Calculation completed at: Tue Feb 17 14:10:27 2026

                    ---------- Happy calculation ----------
                    ---- OpenThermo normally terminated ---
```

### Generated Files

#### basename.UHG (Temperature/Pressure Scan)

**Thermal corrections and total energies:**

```text
Ucorr, Hcorr and Gcorr are in kcal/mol; U, H and G are in a.u.

     T(K)     P(atm)    Ucorr     Hcorr     Gcorr            U                H                G
   298.15     1.000     4.567     4.789     4.123   -76.380162   -76.379940   -76.384063
   323.15     1.000     5.123     5.345     4.567   -76.379506   -76.379284   -76.383851
```

#### basename.SCq (Temperature/Pressure Scan)

**Entropy and partition functions:**

```text
S, CV and CP are in cal/mol/K; q(V=0)/NA and q(bot)/NA are unitless

    T(K)     P(atm)      S         CV        CP       q(V=0)/NA      q(bot)/NA
   298.15     1.000    45.67     12.34     13.56   1.234567e-05   2.345678e-03
   323.15     1.000    48.90     13.45     14.67   1.456789e-05   2.567890e-03
```

#### basename.vibcon (Vibration Contributions)

**Individual mode contributions (when prtvib = 1 or -1):**

```text
Vibrational mode contributions at 298.15 K:

Mode     Freq(cm-1)    ZPE(kJ/mol)    U-T(kJ/mol)    S(J/mol/K)    CV(J/mol/K)
   1        456.7         0.005         0.008         2.34         1.23
   2        789.0         0.008         0.012         3.45         2.34
   ...
```

#### *.otm (OpenThermo Format)

**Native format file (when outotm = 1):**

```text
*E  //Electronic energy (a.u.)
-76.384729

*wavenum  //Wavenumbers (cm^-1)
456.7
789.0
...

*atoms  //Name, mass (amu), X, Y, Z (Angstrom)
C   12.000000   0.000000   0.000000   0.000000
H    1.007825   0.000000   0.000000   1.089000
...
```

## Thermochemistry Calculation Methods

### Standard RRHO Model

**Rigid-Rotor Harmonic Oscillator approximation**

#### Partition Functions

- **Translational**: `q_trans = (2πmkT/h²)^(3/2) * V`
- **Rotational**: `q_rot = √(πI_a I_b I_c (kT/h²)^3) / σ` (non-linear)
- **Vibrational**: `q_vib = ∏(1 - exp(-hν_i/kT))⁻¹`
- **Electronic**: `q_elec = ∑ g_i exp(-ε_i/kT)`

#### Thermodynamic Properties

- **ZPE**: `∑ (1/2) h ν_i`
- **Thermal Energy**: `U - U(0) = RT² ∂(ln Q)/∂T`
- **Entropy**: `S = R ln Q + RT ∂(ln Q)/∂T`
- **Heat Capacity**: `CV = T ∂S/∂T`

### Quasi-RRHO Treatments

#### 1. Truhlar's QRRHO Method (`lowvibmeth = 1`)

**Frequency raising approach:**

- Frequencies below threshold (`ravib`) are raised to `ravib` value
- Applied to: ZPE, thermal energy, entropy, heat capacity
- **Mathematical**: Replace ν < ν_threshold with ν_threshold
- **Advantage**: Simple and computationally efficient

#### 2. Grimme's Interpolation (`lowvibmeth = Grimme`)

**Entropy interpolation between RRHO and free rotor:**

- **Weighting**: `w = 1 / (1 + (ν_threshold/ν)^4)`
- **Entropy**: `S = w × S_RRHO + (1-w) × S_free`
- **Free rotor entropy**: `S_free = R [ 1/2 + ln(√(8π³I kT / h²)) ]`
- **Moment of inertia I**: Uses Bav = 1×10⁻⁴⁴ kg m² (Grimme's original μ'_av, fixed for this method)
- **Threshold**: Configurable via `intpvib` parameter

#### 3. Minenkov's Interpolation (`lowvibmeth = 3`)

**Extended Grimme's method with energy interpolation:**

- **Same entropy treatment** as Grimme's
- **Additional energy interpolation**: `U = (1/α) U_RRHO + (1-1/α) U_free`
- **Where**: `α = 1 + (ν_threshold/ν)^4`
- **Free rotor energy**: `U_free = RT/2`
- **Note**: ZPE = 0 for free rotor contribution

#### 4. Head-Gordon's Interpolation (`lowvibmeth = 4`)

**Energy interpolation with optional entropy interpolation (Li et al., 2015):**

- **Energy interpolation**: `U = w × U_RRHO + (1-w) × U_free`
- **Weighting**: `w = 1 / (1 + (ν_threshold/ν)^4)` (same damping function as Grimme/Minenkov)
- **Free rotor energy**: `U_free = RT/2`
- **ZPE handling**: ZPE is folded into the interpolated total vibrational energy (damped for low-frequency modes, as per eq. 4 of the paper). Separate ZPE is not displayed.
- **Cv interpolation**: `Cv = w × Cv_HO + (1-w) × R/2`
- **Optional entropy interpolation** (`hg_entropy = true`, default): uses the same Grimme free-rotor entropy formula
- **Entropy off** (`hg_entropy = false`): entropy uses the standard harmonic oscillator model (paper's original behavior)
- **Bav**: Configurable via `-bav` option (only for this method):
  - `qchem` (default): I_av = 2.79928×10⁻⁴⁶ kg m² (B_av = 1 cm⁻¹, Q-Chem manual)
  - `grimme`: I_av = 1×10⁻⁴⁴ kg m² (Grimme 2012)
- **Threshold**: Configurable via `intpvib` parameter
- **Reference**: Li, Guo, Head-Gordon, Bell, *J. Phys. Chem. C*, 2015

### Low Frequency Detection

- **Automatic threshold**: Based on `ravib` or `intpvib` parameters
- **Default threshold**: 100.0 cm⁻¹
- **Treatment scope**: Applied to all vibrational modes below threshold

## Advanced Features

### Temperature and Pressure Scanning

```bash
# Temperature range
./OpenThermo molecule.log -T 200 400 25

# Pressure range
./OpenThermo molecule.log -P 0.5 2.0 0.2

# Combined range
./OpenThermo molecule.log -T 273 373 50 -P 0.5 2.0 0.5
```

### Mass Modifications

**Custom atomic masses in settings.ini:**

```ini
modmass
1 H 1.007825  # Atom 1: Hydrogen with specific mass
2 C 12.0      # Atom 2: Carbon-12
3 O 15.994915 # Atom 3: Oxygen-16
```

### Symmetry Analysis

- **Automatic detection** of point groups
- **Manual override** with `-PGname` option
- **Rotational symmetry** number calculation
- **Linear molecule** detection

### Imaginary Frequency Handling

- **Automatic detection** of imaginary frequencies
- **Optional treatment** as real frequencies
- **Configurable threshold** with `-imagreal` option

### VASP Energy Selection

- **Energy line format**: In VASP OUTCAR files, the energy line contains two values:
  
  ```
  energy  without entropy=      -27.39346935  energy(sigma->0) =      -27.39346935
  ```
- **Selection options**:
  - `extrape = false` (default): Use first energy value (4th token)
  - `extrape = true`: Use last energy value (final token)
- **Configuration**: Set in `settings.ini` or use default behavior

### Batch Processing

**Process multiple files:**

```bash
# Create file list
echo "molecule-1.log" > batch.txt
echo "molecule-2.out" >> batch.txt
echo "molecule-3.otm" >> batch.txt

# Process batch
./build/OpenThermo batch.txt
```

## Usage Examples

### Basic Calculations

```bash
# Standard calculation at 298.15 K, 1 atm
./build/OpenThermo water.log

# Custom conditions
./build/OpenThermo methane.out -T 300 -P 2.0

# High precision calculation
./build/OpenThermo benzene.otm -T 298.15 -lowvibmeth 2 -sclZPE 0.98
```

### Advanced Calculations

```bash
# Temperature scan with Grimme's method
./build/OpenThermo molecule.log -T 200 400 25 -lowvibmeth 2

# Pressure scan with custom scaling
./build/OpenThermo molecule.log -P 0.1 10 0.5 -sclS 0.99 -sclCV 0.99

# Transition state calculation
./build/OpenThermo ts.out -imagreal 100 -lowvibmeth 1 -ravib 50

# Condensed phase calculation
./OpenThermo crystal.out -ipmode 1 -conc "1.0 M"

# Minimal output (final data only)
./build/OpenThermo molecule.log -prtlevel 0

# Verbose output with per-mode vibration detail
./build/OpenThermo molecule.log -prtlevel 3
```

### Batch Processing Examples

```bash
# Process all .log files in directory
ls *.log > files.txt
./OpenThermo files.txt

# Custom analysis for multiple molecules
./OpenThermo molecules.txt -T 298 -lowvibmeth 3 -prtvib 1 -outotm 1
```

### Settings File Examples

```ini
# Standard thermochemistry
T = 298.15
P = 1.0
lowvibmeth = Grimme
sclZPE = 1.0

# High-precision calculation
T = 298.15
P = 1.0
lowvibmeth = Minenkov
sclZPE = 0.98
sclheat = 0.99
sclS = 0.99
sclCV = 0.99

# Head-Gordon quasi-RRHO (energy + entropy)
T = 298.15
P = 1.0
lowvibmeth = HeadGordon
intpvib = 100.0
hg_entropy = true

# Low temperature analysis
T = 100
P = 1.0
lowvibmeth = 1
ravib = 50.0
```

## Troubleshooting

### Common Issues

#### 1. File Not Found

**Error**: `Error: Unable to find molecule.log`
**Solution**:

- Check file path and spelling
- Ensure file is in current directory or provide full path
- For .txt batch files, verify all listed files exist

#### 2. Unrecognized Program Format

**Error**: `Error: Unable to identify the program that generated this file`
**Solution**:

- Verify the input file contains frequency analysis output
- Check if the quantum chemistry program is supported
- Ensure the calculation completed successfully

#### 3. Invalid Arguments

**Error**: `Error: Invalid value for -T`
**Solution**:

- Check argument syntax: `-T 298.15` (single) or `-T 200 400 25` (range)
- Ensure numeric values are valid
- Use `--help-T` for detailed syntax

#### 4. Compilation Errors

**Error**: Undefined references or compilation failures
**Solution**:

- Ensure all source files are present
- Check compiler compatibility (GCC 7+, Intel, Clang)
- Run `make clean && make` to rebuild
- Verify Makefile is in the correct directory

#### 5. Settings File Issues

**Warning**: `Warning: settings.ini cannot be found`
**Solution**:

- Create `settings.ini` in current directory
- Set `OPENTHERMOPATH` environment variable
- Use `-noset` to skip settings loading

### Performance Issues

#### Large Files

- Enable batch processing with smaller chunks
- Build with OpenMP (`make OPENMP=1`) for temperature/pressure scans
- Use `-omp-threads N` to control thread count (default: half of physical cores)
- For dedicated compute nodes, use `-omp-threads` with the full core count for maximum performance
- On shared HPC headnodes, the conservative default (half cores) prevents overloading

### Validation

#### Compare Results

- Use different low frequency treatments (`lowvibmeth 0,1,2,3`)
- Compare with literature values
- Validate against experimental data

#### Debug Output

- Use `-prtlevel 2` or `-prtlevel 3` for more detailed output
- Use `-prtvib 1` for detailed vibration analysis (auto-enabled at level 3)
- Check `basename.vibcon` for individual contributions
- Review scan files for temperature/pressure dependence

## Technical Details

### Physical Constants

```cpp
const double R = 8.3144648;           // Gas constant (J/mol/K)
const double kb = 1.3806503e-23;      // Boltzmann constant (J/K)
const double NA = 6.02214179e23;      // Avogadro number
const double h = 6.62606896e-34;      // Planck constant (J·s)
const double au2kJ_mol = 2625.5;      // Hartree to kJ/mol
const double au2cm_1 = 219474.6363;   // Hartree to cm⁻¹
```

### Algorithm Overview

#### 1. Input Processing

- File format detection
- Molecular geometry extraction
- Frequency and energy parsing
- Symmetry analysis

#### 2. Thermochemistry Calculation

- Mass assignment and moment of inertia calculation
- Partition function evaluation
- Thermodynamic property computation
- Low frequency treatment application

#### 3. Output Generation

- Console results display
- File output creation
- Error handling and validation

### Memory Management

- **Dynamic arrays** for molecular data
- **Vector containers** for flexible sizing
- **Efficient parsing** of large output files
- **Minimal memory footprint** for typical calculations

### Error Handling

- **Comprehensive validation** of input parameters
- **Graceful degradation** for invalid data
- **Detailed error messages** with suggestions
- **File existence checks** before processing

### Performance Optimizations

- **OpenMP parallelization** (optional) for temperature/pressure scans, vibrational loops, and symmetry calculations
- **Auto-strategy selection**: OpenThermo chooses between parallelizing the T/P scan (outer) or the vibrational loop (inner) based on workload characteristics
- **HPC-friendly defaults**: Uses half of physical CPU cores to avoid headnode overload
- **Efficient algorithms** for partition function calculation
- **Optimized mathematical operations**

#### OpenMP Parallelization

When built with OpenMP enabled (`make OPENMP=1` or `cmake -DENABLE_OPENMP=ON`), OpenThermo automatically parallelizes computation with conservative, HPC-friendly defaults.

**Thread management:**

- **Default**: Half of detected physical CPU cores (minimum 1)
- **HPC scheduler**: Automatically detects allocated CPUs from `SLURM_CPUS_PER_TASK`, `PBS_NUM_PPN`, `PBS_NP`, `NSLOTS`, `LSB_DJOB_NUMPROC` — when detected, uses half of the scheduler allocation instead of hardware cores
- **Override**: Use `-omp-threads N` (CLI) or `omp-threads = N` (settings.ini); CLI takes precedence
- **Safety clamp**: If `N` exceeds effective cores (scheduler or hardware), the default is used with a warning
- **Rationale**: Prevents overloading shared HPC headnodes while leaving resources for other users and system processes

**Strategy auto-selection:**

- **Outer strategy**: When many T/P scan points are available (≥ nthreads), the T/P loop is parallelized and vibrational calculations run serially inside each iteration
- **Inner strategy**: When few T/P points but many vibrational frequencies (>50), the vibrational loop is parallelized instead
- The strategy is selected automatically — users don't need to understand the internal parallelism patterns

**Examples:**

```bash
# Default: uses half of physical cores, auto-selects strategy
./build/OpenThermo molecule.log -T 100 5000 1

# Explicit thread count for a dedicated compute node
./build/OpenThermo molecule.log -T 100 5000 1 -omp-threads 8

# Verbose output shows thread count and strategy selection
./build/OpenThermo molecule.log -T 100 500 10 -prtlevel 2
```

**Example output notifications:**

```
# Default (4 physical cores detected, no scheduler)
OpenMP threads: 2 (default: half of 4 physical cores). Use -omp-threads N to override.
Strategy: outer (parallel T/P scan, 32 points, 30 frequencies)

# On HPC with SLURM_CPUS_PER_TASK=8
OpenMP threads: 4 (default: half of 8 scheduler-allocated CPUs). Use -omp-threads N to override.
Strategy: outer (parallel T/P scan, 32 points, 30 frequencies)

# Explicit valid request
OpenMP threads: 4 (user-specified)
Strategy: inner (serial T/P scan, parallel vibrational, 200 frequencies)

# Excessive request (clamped)
Warning: Requested 100 threads exceeds 4 physical cores. Using default (2).
Use -omp-threads <= 4 to use more threads.
```

The code compiles and runs correctly both with and without OpenMP.

## Testing

OpenThermo includes a regression test suite that validates thermodynamic calculations against reference values across all supported formats.

### Running Tests

```bash
# Run the full test suite
make test

# Run with verbose output (shows details on failure)
python3 tests/run_tests.py --verbose

# Run a single test by ID
python3 tests/run_tests.py --test gaussian_h2co_freq
```

### Test Coverage

The suite includes 25 tests organized by category:

**Format coverage** (10 tests) — one per supported QC program at default settings:

| ID | Format | Molecule |
|----|--------|----------|
| `gaussian_h2co_freq` | Gaussian | H2CO |
| `gaussian_large` | Gaussian | Large molecule |
| `gaussian_c2h5` | Gaussian | C2H5 |
| `orca_optfreq` | ORCA | Opt+freq |
| `orca_ethanol` | ORCA | Ethanol |
| `gamess_ch4` | GAMESS | CH4 |
| `gamess_h2co` | GAMESS | H2CO |
| `nwchem_nh2coh` | NWChem | NH2COH |
| `cp2k_co2_freq` | CP2K | CO2 |
| `vasp_co_freq` | VASP | CO |

**Option coverage** (9 tests) — different temperatures, pressures, low-frequency methods, scaling, concentration, imaginary frequency handling, and condensed phase mode.

**Output file tests** (3 tests) — verify `.otm`, `.vibcon`, and temperature scan (`.UHG`/`.SCq`) file generation.

**Error handling** (2 tests) — non-existent and invalid input files.

**Batch/ensemble** (1 test) — conformer ensemble with Boltzmann averaging.

Each regression test compares thermodynamic values against reference values with defined tolerances.

### Regenerating Reference Values

After modifying calculation logic, regenerate the reference file:

```bash
make test-generate
```

This runs OpenThermo on all test inputs and writes the results to `tests/reference_values.json`.

### CI

Tests run automatically on every push via GitHub Actions:
- **Linux and macOS**: Build + regression tests (serial and OpenMP at 1 and 4 threads)
- **Windows**: MSYS2/MinGW build + regression tests (serial and OpenMP)
- **AddressSanitizer**: Separate build with ASan on Linux to detect memory errors
- **ThreadSanitizer**: OpenMP build with TSan on Linux to detect data races

## Contributing

We welcome contributions to OpenThermo! Here's how you can help:

### Development Setup

1. **Fork the repository** on GitHub

2. **Clone your fork:**
   
   ```bash
   git clone https://github.com/your-username/openthermo.git
   cd openthermo
   ```

3. **Set up development environment:**
   
   ```bash
   make clean && make debug  # Build with debug symbols
   ./build/OpenThermo --help # Verify installation
   ```

4. **Create a feature branch:**
   
   ```bash
   git checkout -b feature/your-feature-name
   ```

### Code Standards

- **Language**: C++17 with clang-format styling
- **Documentation**: Doxygen comments for all public APIs
- **Testing**: Add tests for new features
- **Style**: Follow LLVM coding standards
- **Commits**: Use descriptive commit messages

### Contribution Guidelines

1. **Code Quality**:
   
   - Ensure all code compiles without warnings
   - Add comprehensive documentation
   - Follow existing code patterns and conventions
   - Test your changes thoroughly

2. **Documentation**:
   
   - Update README.md for new features
   - Add Doxygen comments to new functions
   - Update examples if needed

3. **Testing**:
   
   - Test with various quantum chemistry output files
   - Verify results against known reference values
   - Check edge cases and error conditions

4. **Pull Requests**:
   
   - Provide clear description of changes
   - Reference any related issues
   - Ensure CI/CD checks pass
   - Request review from maintainers

### Reporting Issues

- **Important Note**: OpenThermo is not fully tested agaist several ouputs
- **Bug Reports**: Use the issue template with sample input files
- **Feature Requests**: Describe the use case and expected behavior
- **Performance Issues**: Include system specs and benchmark data

### Development Commands

```bash
# Build and test
make clean && make
make test

# Generate documentation
doxygen Doxyfile

# Format code
clang-format -i src/*.cpp src/*.h

# Run static analysis
clang-tidy src/*.cpp -- -std=c++17
```

## Cite

Cite this project if you use OpenThermo for your research:

```bibtex
@article{pham_openthermo_2025,
    title = {{OpenThermo} {A} {Comprehensive} {C}++ {Program} for {Calculation} of {Thermochemical} {Properties}},
    url = {https://www.researchgate.net/doi/10.13140/RG.2.2.22380.63363},
    doi = {10.13140/RG.2.2.22380.63363},
    urldate = {2025},
    author = {Pham, Le Nhan},
    year = {2025},
    note = {Publisher: Preprint},
}
```

## License

OpenThermo is licensed under the **MIT License**.

### Acknowledgments

- **Le Nhan Pham**: Developer and maintainer
- **Romain Despoullains**: developer
- **Open-source community**: For contributions and feedback
- **Quantum chemistry community**: For validation and testing

---

**OpenThermo v0.001.4** - High-performance molecular thermochemistry calculations
Developed by Le Nhan Pham | [GitHub](https://github.com/lenhanpham/OpenThermo)

For more information, visit the project documentation or use `OpenThermo --help`
