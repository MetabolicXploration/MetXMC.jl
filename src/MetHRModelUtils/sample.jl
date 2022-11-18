## ------------------------------------------------------------------
function _sample!(hrm::MetHRModel)
    
    x = hrm.x
    v = hrm.v
    base = hrm.base
    lb = hrm.net.lb
    ub = hrm.net.ub
    
    # preprocessing: find base
    n = length(v)
    k = size(base,2)

    # pick a random direction
    v .= base * randn(hrm.rng, k)
    dx = v/norm(v)
    # compute intersection
    l = maximum(min((lb[i]-x[i])/dx[i], (ub[i]-x[i])/dx[i]) for i=1:n if dx[i] != 0)
    u = minimum(max((lb[i]-x[i])/dx[i], (ub[i]-x[i])/dx[i]) for i=1:n if dx[i] != 0)
    # find a random point in the intersection
    t = l + (u-l) * rand(hrm.rng)
    x .+= t * dx
    return x
end

function _sample!(hrm::MetHRModel, idx::Int, sn::Int, dt::Int)
    h = zeros(sn)
    for si in 1:sn
        for t in 1:dt; _sample!(hrm); end
        h[si] = _sample!(hrm)[idx]
    end
    return h
end

_sample!(hrm::MetHRModel, rxn::Int) = _sample!(hrm)[rxn]
_sample!(hrm::MetHRModel, ::Nothing) = _sample!(hrm)

function _sample!(hrm::MetHRModel, idxs, sn::Int, dt::Int)
    N = length(idxs)
    h = zeros(sn, N)
    for si in 1:sn
        for t in 1:dt; _sample!(hrm); end
        h[si, :] .= _sample!(hrm)[idxs]
    end
    return h
end

function _sample!(hrm::MetHRModel, ::Nothing, sn::Int, dt::Int)
    _, N = size(hrm.net.S)
    h = zeros(sn, N)
    for si in 1:sn
        for t in 1:dt; _sample!(hrm); end
        h[si, :] .= _sample!(hrm)
    end
    return h
end

## ------------------------------------------------------------------
export sample!
function sample!(hrm::MetHRModel; 
        dt::Int = 3, 
        nsamples = 1,
        rxns = nothing
    ) 
    rxns = isnothing(rxns) ? rxns : rxnindex(hrm.net, rxns)
    return nsamples == 1 ? 
        _sample!(hrm, rxns) : 
        _sample!(hrm, rxns, nsamples, dt)
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
