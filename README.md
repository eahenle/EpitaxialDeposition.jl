# EpitaxialDeposition.jl

## Chemical Reactors Term Project

### Ian Harreschou and Adrian Henle

#### Purpose

`EpitaxialDeposition.jl` is a pure-Julia package for simulating the epitaxial growth of single-crystalline silicon thin films on silicon wafers by a chemical vapor deposition process.  The heavy lifting is done by the `Catalyst.jl` and `JFVM.jl` packages for chemical kinetics and mass transport, respectively.

#### Setup

`EpitaxialDeposition.jl` is not a package on the Julia general registry.  To obtain it, clone the repo:

```
$ git clone https://github.com/eahenle/EpitaxialDeposition.jl
```

The `JFVM.jl` dependency is also not registered, and must be obtained from its GitHub repo.  To set up this and all other dependencies, run:

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

**This package is a demonstration only.**  The assumptions are many, and the accuracy is untested.  Attempting to use this package for actual chemistry could result in harm to equipment or personnel.
