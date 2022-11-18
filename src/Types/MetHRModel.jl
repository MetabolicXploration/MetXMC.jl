export MetHRModel
struct MetHRModel{T<: AbstractFloat}
    # net
    net::MetNet

    # sample! stuf
    x::Vector{T}
    v::Vector{T}
    base::Matrix{T}
    rng::MersenneTwister


    # extras
    extras::Dict
end

function MetHRModel(net::MetNet, jump_args...; rng = Random.GLOBAL_RNG)

    x = warmup(net.S, net.b, net.lb, net.ub, jump_args...)
    v = zeros(length(x))
    base = nullspace(_dense(net.S))
    
    return MetHRModel{eltype(net.S)}(net, x, v, base, rng, Dict())
end