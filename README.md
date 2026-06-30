# AIMPRO Analysis Tools

A collection of utilities for analysing AIMPRO input and output files. This package is intended to be imported as a Python module, providing an abstraction layer for low-level file parsing and data handling in AIMPRO workflows.

## Features

- Custom data structures representing individual atoms and atomic systems, storing properties such as atom index, species, and position.
- Methods for rotating atomic systems about arbitrary axes, angles, and centres of rotation.
- Parsing of lattice vectors and lattice constants from user-specified optimisation iterations.
- Extraction of atomic positions in both Cartesian and fractional coordinates of the lattice vectors for selected optimisation steps.
- Retrieval of key calculation parameters and results, including:
  - Final total energy
  - Final maximum atomic force on any atom in the system
  - Plane-wave cut-off energy
  - Additional simulation metadata.

## Purpose

The module is designed to simplify the development of higher-level analysis scripts by providing reusable tools for interacting with AIMPRO calculations in a consistent and automated manner.

Please credit the code author, James Ramsey, in all works arising from the use of this software.

Copyright (c) J Ramsey 2026 All Rights Reserved.

