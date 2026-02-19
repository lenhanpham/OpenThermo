.. OpenThermo documentation master file, created by
    sphinx-quickstart on Thu Sep 05 12:00:00 2025.
    You can adapt this file completely to your liking, but it should at least
    contain the root `toctree` directive.

OpenThermo User Manual
=======================

**OpenThermo** is a comprehensive C++17 program for calculating molecular thermochemistry properties from quantum chemistry output files (Gaussian, Orca, GAMESS, NWCHEM, CP2K, VASP). OpenThermo implements state-of-the-art methods for computing thermodynamic quantities including Gibbs free energy, enthalpy, entropy, and heat capacities using statistical mechanics.

.. image:: _static/ot-logo.svg
    :alt: OpenThermo Logo
    :align: center
    :scale: 75%

Overview
--------

OpenThermo provides a comprehensive suite of tools for molecular thermochemistry calculations:

* **Multi-format Support**: Gaussian, ORCA, GAMESS-US, NWChem, CP2K, VASP
* **Advanced Thermochemistry**: Standard RRHO and quasi-RRHO treatments for low-frequency modes
* **Statistical Mechanics**: Rigorous implementation of partition functions and thermodynamic properties
* **Symmetry Analysis**: Automatic point group detection and rotational symmetry number calculation
* **High Performance**: Optimized C++17 implementation with optional OpenMP parallelization
* **Batch Processing**: Multi-file analysis with ensemble averaging

Key Features
------------

**Scientific Capabilities**
   - Calculate thermodynamic properties (Gibbs free energy, enthalpy, entropy, heat capacity)
   - Support for multiple low-frequency treatment methods (RRHO, Truhlar, Grimme, Minenkov, Head-Gordon)
   - Automatic symmetry detection and rotational symmetry number calculation
   - Temperature and pressure scanning capabilities

**Performance & Safety**
   - OpenMP parallelization for temperature/pressure scans and vibrational calculations
   - Auto-strategy selection for optimal parallelization
   - HPC-friendly defaults (uses half of physical CPU cores)
   - Memory-efficient algorithms for large systems

**Workflow Integration**
   - Command-line interface with extensive options
   - Configuration file support (settings.ini)
   - Batch processing of multiple files
   - Multiple output formats (console, text files, native .otm format)

Getting Started
---------------

New to OpenThermo? Start here:

1. :doc:`installation` - Install the software on your system
2. :doc:`usage` - Learn how to use all features with examples
3. :doc:`developer` - Developer guide for contributors

Quick Start
-----------

.. code-block:: bash

   # Basic usage - calculate thermochemistry from Gaussian output
   ./build/OpenThermo molecule.log

   # Custom temperature and pressure
   ./build/OpenThermo molecule.log -T 300 -P 2.0

   # Temperature scan with Grimme's method
   ./build/OpenThermo molecule.log -T 200 400 25 -lowvibmeth Grimme

   # Get help
   ./build/OpenThermo --help

Important Note
--------------

**The project is still in the early stage and not fully tested. Therefore, errors and inaccuracy may happen. Users are suggested to check calculated data against original data from outputs of corresponding quantum chemical programs (with default temperature, concentration, and pressure)**

A graphical user interface (GUI) version is under development and can be found here: 
`OpenThermoGUI <https://github.com/lenhanpham/OpenThermoGUI>`_

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   installation
   usage

.. toctree::
   :maxdepth: 2
   :caption: Developer Guide:

   developer

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Support & Contributing
======================

For bug reports, feature requests, or questions:

- **Issues**: `GitHub Issues <https://github.com/lenhanpham/OpenThermo/issues>`_
- **Documentation**: This user manual
- **Help**: Use ``OpenThermo --help`` for command-line help

License
=======

OpenThermo is released under the MIT License. See the LICENSE file for details.

Version Information
====================

Current Version: **v0.001.4**

- **v0.001.4**: Current release featuring multi-format support, advanced thermochemistry methods, OpenMP parallelization, and comprehensive testing

.. note::
    This documentation is for OpenThermo v0.001.4. For older versions, please refer to the archived documentation.
