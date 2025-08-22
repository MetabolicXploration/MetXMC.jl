## ------------------------------------------------------------------
function _walk!(hrm::HRModel{T}, rng) where {T}
    
    v::Vector{T} = hrm.v            # position
    dv::Vector{T} = hrm.dv          # direction
    r::Vector{T} = hrm.r            # rand direction
    base::Matrix{T} = hrm.base      # nullspace(lep.S)
    lb::Vector{T} = hrm.lep.lb
    ub::Vector{T} = hrm.lep.ub
    # damp::Float64 = config(hrm, :damp, 1.0)::Float64

    n = length(v)
    k = size(base,2)

    while true # keep trying till success
        # step
        # The direction is uniform in the free variable space
        # 
        r .= base * randn(rng, k)      
        dv .= r / norm(r)         # |r| = 1

        # compute feasible region of x1 = x0 .+ λ * dv line
        # that is, λ ∈ [λ0, λ1]
        # for each dimension, you have a λ0 and λ1 which do not
        # break it constraint
        # We need to find the smaller [λ0, λ1] interval, so,
        # no constraint is broken across the sampling line.
        λ0 = -Inf
        λ1 = Inf
        @inbounds for i in 1:n
            dvi = dv[i]
            dvi == 0 && continue
            vi = v[i]
            _λ0 = (lb[i] - vi)/dvi
            _λ1 = (ub[i] - vi)/dvi
            if dvi > 0 
                λ0 = max(_λ0, λ0) # _λ0 <= 0
                λ1 = min(_λ1, λ1) # _λ1 >= 0
            else
                # dvi < 0
                # reverse direction
                λ0 = max(_λ1, λ0) # _λ1 <= 0
                λ1 = min(_λ0, λ1) # _λ0 >= 0
            end
        end
        # fail to resample => repeat
        isinf(λ0) && continue
        isinf(λ1) && continue

        # find a random point in the intersection
        λ = λ0 + (λ1-λ0) * rand(rng)
        # v .+= λ * dv * damp # here damp prevent explaring the close to boundary areas
        v .+= λ * dv
        return v
    end
end


# ## ------------------------------------------------------------------
# function _walk!(hrm::HRModel{T}, rng) where {T}
    
#     v::Vector{T} = hrm.v            # position
#     dv::Vector{T} = hrm.dv          # step
#     r::Vector{T} = hrm.r            # rand direction
#     base::Matrix{T} = hrm.base      # nullspace(lep.S)
#     lb::Vector{T} = hrm.lep.lb
#     ub::Vector{T} = hrm.lep.ub
#     damp::Float64 = config(hrm, :damp, 1.0)::Float64
    
#     n = length(v)
#     k = size(base,2)

#     # As described at https://doi.org/10.1371/journal.pone.0122670
#     # 1.  Choose a uniformly distributed direction, that is, a point extracted from the uniform distribution on the D-dimensional unit sphere. This can be done with the Marsaglia method, i.e. by generating D independent gaussian random variables it with zero mean and unit variance, and then normalizing the vector to unit length;
#     r .= base * randn(rng, k) # sample solution
#     dv .= r / norm(r)         # |r| = 1
#     # compute intersection
#     @inbounds l::T = maximum(
#         min(
#             (lb[i]-v[i])/dv[i], 
#             (ub[i]-v[i])/dv[i]
#         ) 
#         for i=1:n if !iszero(dv[i])
#     )
#     @inbounds u::T = minimum(
#         max(
#             (lb[i]-v[i])/dv[i], 
#             (ub[i]-v[i])/dv[i]
#         ) 
#         for i=1:n if !iszero(dv[i])
#     )

#     # find a random point in the intersection
#     λ = l + (u-l) * rand(rng)
#     v .+= λ * dv * damp
#     return v
# end

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
