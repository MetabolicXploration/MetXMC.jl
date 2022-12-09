export MetHRModel
struct MetHRModel{T<: AbstractFloat}
    # net
    net::MetNet

    # sample! stuf
    x::Vector{T}
    dx::Vector{T}
    v::Vector{T}
    base::Matrix{T}
    rng::Random.AbstractRNG

    # extras
    extras::Dict
end

function MetHRModel(net::MetNet, jump_args...; rng = Random.GLOBAL_RNG)

    x = warmup(net.S, net.b, net.lb, net.ub, jump_args...)
    dx = zeros(length(x))
    v = zeros(length(x))
    base = nullspace(_dense(net.S))
    
    return MetHRModel{eltype(net.S)}(net, x, dx, v, base, rng, Dict())
end