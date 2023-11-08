# MechanochemicalPatterns

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://SteffenPL.github.io/MechanochemicalPatterns.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://SteffenPL.github.io/MechanochemicalPatterns.jl/dev/)
[![Build Status](https://github.com/SteffenPL/MechanochemicalPatterns.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SteffenPL/MechanochemicalPatterns.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/SteffenPL/MechanochemicalPatterns.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/SteffenPL/MechanochemicalPatterns.jl)


## Setup

Use the usual approach to `instantiate` the Julia environment.
```julia
using Pkg
Pkg.instantiate()  # or in pkg> mode
```

In addition, there is the file `scripts\config.toml` which is not shipped with the package
for local configuration (such as paths). Copy the file to `scripts\config_template.toml` to `scripts\config.toml` and edit the file to your needs.

The `data_folder` (see config file) refers to a shared folder which collects all 
experimental datasets. As these files are large we use other cloud drives to share these datasets (Google Drive or OneDrive).