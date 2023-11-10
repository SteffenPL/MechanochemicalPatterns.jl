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

### Local configuration file

To avoid that changes in the repository depend on the computer one is using, **we use a local configuration file which is not syncronized with the repository.**


To initialize, copy `scripts\config_template.toml` to `scripts\config.toml` and edit the file according to your needs.

Options:
- `data_dir` a shared folder which collects all experimental datasets (shared via Google Drive or OneDrive).
- `input_dir` folder to read parameter files from.
- `output_dir` target folder for generated data like plots, videos and CSV tables. (Should be outside of this git repository.)
- `makie_inline` determines if plots appear in an extra window or inside VS Code.

Please use the following line to obtain the configuration in your scripts: 
```julia
using MechanochemicalPatterns
config = load_config()
```
or 
```julia
using TOML
config = TOML.parsefile("scripts/config.toml")
```

## Guidelines 

- Quantities should always be in units with respect to `μm`, `hours` (and `pg`). 
  - For nondimensional code, please define global scaling factors explicitly in the code, such that one can easily see convert the model to proper units.