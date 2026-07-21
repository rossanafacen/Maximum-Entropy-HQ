module MaximumEntropy

    const fmGeV= 1/0.1973261 
    const invfmGeV= 1/fmGeV
    
    using TensorCast
    using Tullio
    using Symbolics
    using BenchmarkTools
    using NonlinearSolve
    using Interpolations
    using CSV
    using HCubature
    using StaticArrays
    using Fluidum
    
    include("functions_HQ.jl")
    include("2D_ME_HQ_all_hadrons.jl")
    include("2D_ME_HQ_single_particle.jl")
    include("spectra_and_multiplicity.jl")
    include("spectra_and_multiplicity_old.jl")
    
    export charm_diff_current, lambdar, charm_density, lambdan_lambdar
end 
