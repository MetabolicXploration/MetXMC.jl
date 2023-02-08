export EPJoinMvNormal
struct EPJoinMvNormal{T} <: AbstractHitOrDropSampler where {T<:AbstractFloat}
    
    # epm::FluxEPModelT0{T}

    Qi::MultivariateNormal
    be::Vector{T}
    G::Matrix{T}
    
    idxi::Vector{Int}
    idxd::Vector{Int}

    v::Vector{T}
    vi::Vector{T}
    rng::Random.AbstractRNG

    extras::Dict

end

function EPJoinMvNormal(epm::FluxEPModelT0{T}; 
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

    # extras
    # TODO: Add some net extras
    extras = Dict()

    return EPJoinMvNormal{T}(Qi, be, G, idxi, idxd, v, vi, rng, extras)
end