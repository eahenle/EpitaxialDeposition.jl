import Pkg

Pkg.activate(".")
Pkg.add(url="https://github.com/simulkade/JFVM.jl")
Pkg.instantiate()
Pkg.build()

using EpitaxialDeposition
