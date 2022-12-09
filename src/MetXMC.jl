module MetXMC

    using LinearAlgebra
    using SparseArrays, JuMP, GLPK
    using MetXBase
    using MetXBase: _dense
    using MetXOptim
    using Random

    import Distributions

    #! include Types
    include("Types/MetHRModel.jl")
    
    #! include MetHRModelUtils
    include("MetHRModelUtils/sample.jl")

    #! include .

end