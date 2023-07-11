module MetXMC

    using LinearAlgebra
    using SparseArrays, JuMP, GLPK
    using MetXBase
    using MetXBase: _dense, _histogram!
    using MetXOptim
    using MetXEP
    using Random

    import Distributions
    import Distributions: pdf

    #! include Types
    include("Types/0_AbstractHitOrDropSampler.jl")
    include("Types/EPBoxedMvNormal.jl")
    include("Types/EPJoinMvNormal.jl")
    include("Types/EPTiltedMvNormal.jl")
    include("Types/HRModel.jl")
    include("Types/MC0Model.jl")
    
    #! include HRModelBase
    include("HRModelBase/dist_interface.jl")
    include("HRModelBase/hd_interface.jl")
    include("HRModelBase/walk.jl")
    include("HRModelBase/warnup.jl")
    
    #! include MC0ModelBase
    include("MC0ModelBase/dist_interface.jl")
    include("MC0ModelBase/hd_interface.jl")
    include("MC0ModelBase/sample.jl")

    #! include AbstractHitOrDropSamplerBase
    include("AbstractHitOrDropSamplerBase/dist_interface.jl")
    include("AbstractHitOrDropSamplerBase/extras_interface.jl")
    include("AbstractHitOrDropSamplerBase/hd_interface.jl")
    include("AbstractHitOrDropSamplerBase/lep_interface.jl")
    
    #! include EPJoinMvNormalBase
    include("EPJoinMvNormalBase/dist_interface.jl")
    include("EPJoinMvNormalBase/hd_interface.jl")
    
    #! include EPBoxedMvNormalBase
    include("EPBoxedMvNormalBase/dist_interface.jl")
    include("EPBoxedMvNormalBase/hd_interface.jl")
    
    #! include EPTiltedMvNormalBase
    include("EPTiltedMvNormalBase/dist_interface.jl")
    include("EPTiltedMvNormalBase/hd_interface.jl")

    #! include .

    # exports
    @_exportall_non_underscore()

end