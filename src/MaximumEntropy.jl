module MaximumEntropy

    using TensorCast

    using Tullio
    using Symbolics
    using BenchmarkTools
    using NonlinearSolve
    using CSV
    using Fluidum

    #include("maximum_entropy_functions.jl")
    #include("maximum_entropy_functions_2d.jl")
    include("maximum_entropy_functions_LF.jl")
    include("maximum_entropy_functions_HQ.jl")
    include("maximum_entropy_functions_2d_HQ.jl")
    include("maximum_entropy_HQ.jl")
end 
