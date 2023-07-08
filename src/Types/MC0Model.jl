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
    # lep
    elep::EchelonLEPModel

    # sample! stuf
    v::Vector{T}
    vi::Vector{T}

    rng::Random.AbstractRNG

    # Default sampler
    U::MultivariateUniform{T}

    # extras
    extras::Dict

end

function MC0Model(elep::EchelonLEPModel; 
        rng = Random.GLOBAL_RNG
    )
    
    Nf, Nd = length(elep.idxf), length(elep.idxd)
    v, vi = zeros(Nf + Nd), zeros(Nf)

    T = eltype(elep.lep.S)
    U = MultivariateUniform(
        Vector{T}(elep.lep.lb[elep.idxf]), 
        Vector{T}(elep.lep.ub[elep.idxf])
    )
    
    return MC0Model{T}(elep, v, vi, rng, U, Dict())
end

function MC0Model(lep::LEPModel; 
        rng = Random.GLOBAL_RNG, 
        echtol = 1e-10
    )
    elep = EchelonLEPModel(lep; tol = echtol)
    return MC0Model(elep; rng)
end
