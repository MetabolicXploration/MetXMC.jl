## ------------------------------------------------------------------
function _sample!(hrm::MetHRModel{T}, rng) where {T}
    
    x::Vector{T} = hrm.x
    dx::Vector{T} = hrm.dx
    v::Vector{T} = hrm.v
    base::Matrix{T} = hrm.base
    lb::Vector{T} = hrm.net.lb
    ub::Vector{T} = hrm.net.ub
    
    # preprocessing: find base
    n = length(v)
    k = size(base,2)

    # pick a random direction
    v .= base * randn(rng, k)
    dx .= v / norm(v)
    # compute intersection
    @inbounds l::T = maximum(
        min(
            (lb[i]-x[i])/dx[i], 
            (ub[i]-x[i])/dx[i]
        ) 
        for i=1:n if !iszero(dx[i])
    )
    @inbounds u::T = minimum(
        max(
            (lb[i]-x[i])/dx[i], 
            (ub[i]-x[i])/dx[i]
        ) 
        for i=1:n if !iszero(dx[i])
    )

    # l, u = -Inf, Inf
    # @inbounds for i in 1:n
    #     dx[i] == 0 && continue
    #     _l = (lb[i]-x[i])/dx[i]
    #     _u = (ub[i]-x[i])/dx[i]
    #     l = max(l, min(_l, _u))
    #     u = min(u, max(_l, _u))
    # end

    # find a random point in the intersection
    t = l + (u-l) * rand(rng)
    x .+= t * dx * 1e-3
    return x
end

function _sample!(hrm::MetHRModel, idx::Int, sn::Int, dt::Int, rng)
    h = zeros(sn)
    for si in 1:sn
        for t in 1:dt; _sample!(hrm, rng); end
        h[si] = _sample!(hrm, rng)[idx]
    end
    return h
end

_sample!(hrm::MetHRModel, rxn::Int, rng) = _sample!(hrm, rng)[rxn]
_sample!(hrm::MetHRModel, ::Nothing, rng) = _sample!(hrm, rng)

function _sample!(hrm::MetHRModel, idxs, sn::Int, dt::Int, rng)
    N = length(idxs)
    h = zeros(sn, N)
    for si in 1:sn
        for t in 1:dt; _sample!(hrm, rng); end
        h[si, :] .= _sample!(hrm, rng)[idxs]
    end
    return h
end

function _sample!(hrm::MetHRModel, ::Nothing, sn::Int, dt::Int, rng)
    _, N = size(hrm.net.S)
    h = zeros(sn, N)
    for si in 1:sn
        for t in 1:dt; _sample!(hrm, rng); end
        h[si, :] .= _sample!(hrm, rng)
    end
    return h
end

## ------------------------------------------------------------------
export sample!
import Distributions.sample!
function sample!(hrm::MetHRModel, nsamples::Int; 
        dt::Int = 3, 
        rng = hrm.rng,
        rxns = nothing
    ) 
    rxns = isnothing(rxns) ? rxns : rxnindex(hrm.net, rxns)
    return _sample!(hrm, rxns, nsamples, dt, rng)
end

function sample!(hrm::MetHRModel; rxns = nothing, rng = hrm.rng)
    rxns = isnothing(rxns) ? rxns : rxnindex(hrm.net, rxns)
    return _sample!(hrm, rxns, rng)
end

## ------------------------------------------------------------------
function warmup(S, b, lb, ub, jump_args...)
    
    T = eltype(S)
    M, N = size(S)
    c = zeros(T, N)
    opm = FBAFluxOpModel(S, b, lb, ub, c, jump_args...)
    x0 = zeros(T, N)
    c0, c1 = zero(T), one(T)
    
    for i=1:N
        
        set_linear_obj!(opm, i, c1)
        fba!(opm)
        x0 .+= solution(opm) / 2N
        
        set_linear_obj!(opm, i, -c1)
        fba!(opm)
        x0 .+= solution(opm) / 2N

        set_linear_obj!(opm, i, c0)

    end
    return x0
end
