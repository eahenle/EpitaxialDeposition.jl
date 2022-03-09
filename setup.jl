import Pkg

Pkg.activate(".")
Pkg.add("https://github.com/simulkade/JFVM.jl")
Pkg.instantiate()
Pkg.build()

using EpitaxialDeposition
