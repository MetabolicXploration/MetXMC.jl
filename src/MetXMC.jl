module MetXMC

    using LinearAlgebra
    using SparseArrays, JuMP, GLPK
    using MetXBase
    using MetXBase: _dense
    using MetXOptim
    using Random

    import Distributions

    #! include Types
    include("Types/HRModel.jl")
    include("Types/MC0Model.jl")
    
    #! include HRModelUtils
    include("HRModelUtils/interfaces.jl")
    include("HRModelUtils/sample.jl")
    
    #! include MC0ModelUtils
    include("MC0ModelUtils/interfaces.jl")
    include("MC0ModelUtils/sample.jl")

    #! include .

end