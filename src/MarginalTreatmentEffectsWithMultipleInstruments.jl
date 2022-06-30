module MarginalTreatmentEffectsWithMultipleInstruments

    using JuMP
    using Clp
    using DataFrames
    using Parameters
    using LinearAlgebra
    using CSV
    using Suppressor
    using NamedArrays
    using DataStructures

    include("basis.jl")
    include("mtr.jl")
    include("dgp.jl")
    include("target_parameters.jl")
    include("ivlike.jl")
    include("weights.jl")
    include("mutual_consistency.jl")
    include("bounds.jl")
    include("paper_simulations.jl")
    include("util_save.jl")
    include("toplevel.jl")

end
