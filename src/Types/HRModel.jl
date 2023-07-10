# TODO: add burnin samples kwarg option see [[@wilhelmTmvtnormPackageTruncated2010]]

struct HRModel{T} <: AbstractHitOrDropSampler where {T<:AbstractFloat}
    # lep
    lep::LEPModel

    # sample! stuf
    v::Vector{T}                # position
    dv::Vector{T}               # step
    r::Vector{T}                # rand direction
    base::Matrix{T}             # nullspace(lep.S)

    rng::Random.AbstractRNG

    # extras
    extras::Dict
end

function HRModel(lep::LEPModel, jump_args...; rng = Random.GLOBAL_RNG)

    v = _warmup(lep.S, lep.b, lep.lb, lep.ub, jump_args...)
    dv = zeros(length(v))
    r = zeros(length(v))
    base = nullspace(_dense(lep.S))
    
    return HRModel{eltype(lep.S)}(lep, v, dv, r, base, rng, Dict())
end