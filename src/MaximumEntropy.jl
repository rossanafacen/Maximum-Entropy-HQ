module MaximumEntropy

    const fmGeV= 1/0.1973261 
    const invfmGeV= 1/fmGeV
    
    using TensorCast
    using Tullio
    using Symbolics
    using BenchmarkTools
    using NonlinearSolve
    using CSV
    using HCubature
    using StaticArrays
    using Fluidum
    
    #include("functions_HQ.jl")
    include("functions_HQ_equilibrium.jl")
    # include("functions_2d_HQ.jl")
    # include("functions_LF.jl")
    # include("functions_2d_LF.jl")
    # include("functions_LF_no_symmetries.jl")
    # include("lambda_evolution.jl")
    # include("lambda_evolution_2d.jl")
    
    export charm_diff_current, lambdar, charm_density, lambdan_lambdar
end 
