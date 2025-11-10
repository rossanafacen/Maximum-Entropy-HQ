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
    using Fluidum
    

    include("maximum_entropy_functions_LF.jl")
    include("maximum_entropy_functions_HQ.jl")
    include("maximum_entropy_functions_2d_HQ.jl")
    include("lambda_evolution.jl")
    include("lambda_evolution_2d.jl")
    
    export charm_diff_current, lambdar
end 
