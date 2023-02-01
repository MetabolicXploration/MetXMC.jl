export MC0Model
struct MC0Model{T}
    # net
    enet::EchelonMetNet

    # sample! stuf
    x::Vector{T}
    vf::Vector{T}
    rng::Random.AbstractRNG

    # extras
    extras::Dict

end

function MC0Model(enet::EchelonMetNet; 
        rng = Random.GLOBAL_RNG
    )
    
    Nf, Nd = length(enet.idxf), length(enet.idxd)
    x, vf = zeros(Nf + Nd), zeros(Nf)
    
    T = eltype(enet.net.S)
    return MC0Model{T}(enet, x, vf, rng, Dict())
end

function MC0Model(net::MetNet; 
        rng = Random.GLOBAL_RNG, 
        echtol = 1e-10
    )
    enet = EchelonMetNet(net; tol = echtol)
    return MC0Model(enet; rng)
end
