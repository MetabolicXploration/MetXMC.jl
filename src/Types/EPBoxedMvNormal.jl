export EPBoxedMvNormal
struct EPBoxedMvNormal{T} <: AbstractHitOrDropSampler where {T<:AbstractFloat}

    Qi::MultivariateNormal
    be::Vector{T}
    G::Matrix{T}
    
    idxi::Vector{Int}
    idxd::Vector{Int}

    lb::Vector{T}
    ub::Vector{T}

    v::Vector{T}
    vi::Vector{T}
    rng::Random.AbstractRNG

    extras::Dict

end

function EPBoxedMvNormal(epm::FluxEPModelT0{T}; 
        rng = Random.GLOBAL_RNG
    ) where T

    Qi = MultivariateNormal(
        _dense(epm.vi) * epm.scalefact,
        _dense(epm.Σi) * epm.scalefact^2
    )
    be = _dense(epm.be) .* epm.scalefact
    idxi, idxd = epm.idxi, epm.idxd
    G = _dense(epm.G)
    Nd, Ni = size(epm.G)
    vi = zeros(Ni)
    v = zeros(Nd + Ni)

    _lb, _ub = lb(epm), ub(epm) # This is already rescaled
    
    # extras
    # TODO: Add some net extras
    extras = Dict()


    EPBoxedMvNormal{T}(Qi, be, G, idxi, idxd, _lb, _ub, v, vi, rng, extras)
end