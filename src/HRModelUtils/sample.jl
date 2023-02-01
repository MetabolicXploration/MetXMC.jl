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
        optimize!(opm)
        x0 .+= solution(opm) / 2N
        
        set_linear_obj!(opm, i, -c1)
        optimize!(opm)
        x0 .+= solution(opm) / 2N

        set_linear_obj!(opm, i, c0)

    end
    return x0
end

## ------------------------------------------------------------------
function _sample!(hrm::HRModel{T}, rng) where {T}
    
    x::Vector{T} = hrm.x
    dx::Vector{T} = hrm.dx
    v::Vector{T} = hrm.v
    base::Matrix{T} = hrm.base
    lb::Vector{T} = hrm.net.lb
    ub::Vector{T} = hrm.net.ub
    damp::Float64 = config(hrm, :damp, 1.0)::Float64
    
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

    # find a random point in the intersection
    t = l + (u-l) * rand(rng)
    x .+= t * dx * damp
    return x
end

function _sample!(hrm::HRModel, idx::Int, sn::Int, dt::Int, rng)
    h = zeros(sn)
    for si in 1:sn
        for t in 1:dt; _sample!(hrm, rng); end
        h[si] = _sample!(hrm, rng)[idx]
    end
    return h
end

_sample!(hrm::HRModel, rxn::Int, rng) = _sample!(hrm, rng)[rxn]
_sample!(hrm::HRModel, ::Nothing, rng) = _sample!(hrm, rng)

function _sample!(hrm::HRModel, idxs, sn::Int, dt::Int, rng)
    N = length(idxs)
    h = zeros(sn, N)
    for si in 1:sn
        for t in 1:dt; _sample!(hrm, rng); end
        h[si, :] .= _sample!(hrm, rng)[idxs]
    end
    return h
end

function _sample!(hrm::HRModel, ::Nothing, sn::Int, dt::Int, rng)
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
function sample!(onhit::Function, hrm::HRModel, niters::Int;
        rng = hrm.rng
    )
    for _ in 1:niters
        _sample!(hrm, rng)
        onhit(hrm) === true && break
    end
    return nothing
end

function sample!(hrm::HRModel, nsamples::Int; 
        dt::Int = 1, 
        rng = hrm.rng,
        rxns = nothing
    ) 
    rxns = isnothing(rxns) ? rxns : rxnindex(hrm.net, rxns)
    return _sample!(hrm, rxns, nsamples, dt, rng)
end

function sample!(hrm::HRModel; rxns = nothing, rng = hrm.rng)
    rxns = isnothing(rxns) ? rxns : rxnindex(hrm.net, rxns)
    return _sample!(hrm, rxns, rng)
end
