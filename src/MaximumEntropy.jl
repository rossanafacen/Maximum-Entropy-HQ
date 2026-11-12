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
    include("2D_HQ_all_hadrons.jl")
    include("2D_HQ_single_particle.jl")
    include("part_dic_spectra_and_multiplicity.jl")
    
    include("LF\\2D_LF_all_hadrons.jl")
    include("LF\\2D_LF_single_particle.jl")
    include("LF\\3D_LF_all_hadrons.jl")
    
    include("LF\\1D_LF_all_hadrons_LRF.jl")
    include("LF\\2D_LF_all_hadrons_LRF.jl")
    include("LF\\3D_LF_all_hadrons_LRF.jl")
    export charm_diff_current, lambdar, charm_density, lambdan_lambdar, lagrange_multiplier_2D
end 
