
export HRModel
struct HRModel{T} <: AbstractHitOrDropSampler where {T<:AbstractFloat}
    # net
    net::MetNet

    # sample! stuf
    v::Vector{T}                # position
    dv::Vector{T}               # step
    r::Vector{T}                # rand direction
    base::Matrix{T}             # nullspace(net.S)

    rng::Random.AbstractRNG

    # extras
    extras::Dict
end

function HRModel(net::MetNet, jump_args...; rng = Random.GLOBAL_RNG)

    v = _warmup(net.S, net.b, net.lb, net.ub, jump_args...)
    dv = zeros(length(v))
    r = zeros(length(v))
    base = nullspace(_dense(net.S))
    
    return HRModel{eltype(net.S)}(net, v, dv, r, base, rng, Dict())
end