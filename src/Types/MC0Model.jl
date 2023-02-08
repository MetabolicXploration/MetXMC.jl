# ------------------------------------------------------------------
struct MultivariateUniform{T} 
    v0::Vector{T}
    Δv::Vector{T}

    MultivariateUniform(lb::Vector{T}, ub::Vector{T}) where T = new{T}(lb, ub - lb)
end

import Base.rand
rand(rng::Random.AbstractRNG, d::MultivariateUniform) = d.v0 .+ rand(rng, length(d.Δv)) .* d.Δv

import Distributions.pdf
pdf(::MultivariateUniform{T}, x) where T = one(T)

## ------------------------------------------------------------------
export MC0Model
struct MC0Model{T} <: AbstractHitOrDropSampler where {T<:AbstractFloat}
    # net
    enet::EchelonMetNet

    # sample! stuf
    v::Vector{T}
    vi::Vector{T}

    rng::Random.AbstractRNG

    # Default sampler
    U::MultivariateUniform{T}

    # extras
    extras::Dict

end

function MC0Model(enet::EchelonMetNet; 
        rng = Random.GLOBAL_RNG
    )
    
    Nf, Nd = length(enet.idxf), length(enet.idxd)
    v, vi = zeros(Nf + Nd), zeros(Nf)

    T = eltype(enet.net.S)
    U = MultivariateUniform(
        Vector{T}(enet.net.lb[enet.idxf]), 
        Vector{T}(enet.net.ub[enet.idxf])
    )
    
    return MC0Model{T}(enet, v, vi, rng, U, Dict())
end

function MC0Model(net::MetNet; 
        rng = Random.GLOBAL_RNG, 
        echtol = 1e-10
    )
    enet = EchelonMetNet(net; tol = echtol)
    return MC0Model(enet; rng)
end
