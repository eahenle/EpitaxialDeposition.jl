using EpitaxialDeposition, Test
import Aqua

@testset "notebooks run" begin
    @info "Running kinetics notebook"
    @time include("../kinetics.jl")
    @test true

    @info "Running transport notebook"
    @time include("../transport.jl")
    @test true

    @info "Running presentation notebook"
    @time include("../Presentation.jl")
    @test true

    @info "Running coupled notebook"
    @time include("../coupled.jl")
    @test true
end

Aqua.test_ambiguities(EpitaxialDeposition, recursive=false) ##! `recursive=false` b/c Aqua finds an ambiguity IN ITSELF!
Aqua.test_unbound_args(EpitaxialDeposition)
#Aqua.test_undefined_exports(EpitaxialDeposition) ##! name clashes in deps
Aqua.test_project_extras(EpitaxialDeposition)
Aqua.test_stale_deps(EpitaxialDeposition, ignore=[:LsqFit])
Aqua.test_deps_compat(EpitaxialDeposition)
Aqua.test_project_toml_formatting(EpitaxialDeposition)
