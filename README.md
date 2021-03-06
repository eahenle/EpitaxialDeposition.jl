# EpitaxialDeposition.jl

## Chemical Reactors Term Project

### Ian Harreschou and Adrian Henle

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![.github/workflows/ci_build.yml](https://github.com/eahenle/EpitaxialDeposition.jl/actions/workflows/ci_build.yml/badge.svg)](https://github.com/eahenle/EpitaxialDeposition.jl/actions/workflows/ci_build.yml)

#### Purpose

`EpitaxialDeposition.jl` is a pure-Julia package for simulating the epitaxial growth of single-crystalline silicon thin films on silicon wafers by a chemical vapor deposition process.  The heavy lifting is done by the [`Catalyst.jl`](https://github.com/SciML/Catalyst.jl) and [`JFVM.jl`](https://github.com/simulkade/JFVM.jl) packages for chemical kinetics and mass transport, respectively.

#### Setup

`EpitaxialDeposition.jl` is not a package on the Julia general registry.  To obtain it, clone the repo:

```
$ git clone https://github.com/eahenle/EpitaxialDeposition.jl
```

The [`JFVM.jl`](https://github.com/simulkade/JFVM.jl) dependency is also not registered, and must be obtained from its GitHub repo.  To easily set up this and all other dependencies, run:

```
$ julia setup.jl
```

#### Testing

To check that the package is working, do the following in the REPL `Pkg` mode:

```julia
Pkg> activate .
Pkg> test
```

#### Disclaimer

**This package is a demonstration only.**  Attempting to use this package as-is for actual chemistry could result in harm to equipment or personnel.
