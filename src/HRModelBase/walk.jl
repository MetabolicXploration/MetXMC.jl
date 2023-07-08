## ------------------------------------------------------------------
function _walk!(hrm::HRModel{T}, rng) where {T}
    
    v::Vector{T} = hrm.v            # position
    dv::Vector{T} = hrm.dv          # step
    r::Vector{T} = hrm.r            # rand direction
    base::Matrix{T} = hrm.base      # nullspace(lep.S)
    lb::Vector{T} = hrm.lep.lb
    ub::Vector{T} = hrm.lep.ub
    damp::Float64 = config(hrm, :damp, 1.0)::Float64
    
    n = length(v)
    k = size(base,2)

    # pick a random direction
    r .= base * randn(rng, k)
    dv .= r / norm(r)
    # compute intersection
    @inbounds l::T = maximum(
        min(
            (lb[i]-v[i])/dv[i], 
            (ub[i]-v[i])/dv[i]
        ) 
        for i=1:n if !iszero(dv[i])
    )
    @inbounds u::T = minimum(
        max(
            (lb[i]-v[i])/dv[i], 
            (ub[i]-v[i])/dv[i]
        ) 
        for i=1:n if !iszero(dv[i])
    )

    # find a random point in the intersection
    λ = l + (u-l) * rand(rng)
    v .+= λ * dv * damp
    return v
end

# ## ------------------------------------------------------------------
# export sample!
# import Distributions.sample!
# function sample!(onhit::Function, hrm::HRModel, niters::Int;
#         rng = hrm.rng
#     )
#     for _ in 1:niters
#         _sample!(hrm, rng)
#         onhit(hrm) === true && break
#     end
#     return nothing
# end

# function sample!(hrm::HRModel, nsamples::Int; 
#         dt::Int = 1, 
#         rng = hrm.rng,
#         rxns = nothing
#     ) 
#     rxns = isnothing(rxns) ? rxns : rxnindex(hrm.lep, rxns)
#     return _sample!(hrm, rxns, nsamples, dt, rng)
# end

# function sample!(hrm::HRModel; rxns = nothing, rng = hrm.rng)
#     rxns = isnothing(rxns) ? rxns : rxnindex(hrm.lep, rxns)
#     return _sample!(hrm, rxns, rng)
# end
