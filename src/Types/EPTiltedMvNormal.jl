struct EPTiltedMvNormal{T} <: AbstractHitOrDropSampler where {T<:AbstractFloat}

    Qi::MultivariateNormal
    be::Vector{T}
    G::Matrix{T}
    
    idxi::Vector{Int}
    idxd::Vector{Int}
    
    rxni::Int
    lb::T
    ub::T

    v::Vector{T}
    vi::Vector{T}
    rng::Random.AbstractRNG

    extras::Dict

end

function EPTiltedMvNormal(epm::FluxEPModelT0{T}, rxn; 
        rng = Random.GLOBAL_RNG
    ) where T

    rxni = rxnindex(epm, rxn)
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

    _lb, _ub = lb(epm, rxni), ub(epm, rxni) # This is already rescaled
    
    # extras
    # TODO: Add some lep extras
    extras = Dict()

    EPTiltedMvNormal{T}(Qi, be, G, idxi, idxd, rxni, _lb, _ub, v, vi, rng, extras)
end