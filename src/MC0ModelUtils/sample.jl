## ------------------------------------------------------------------
function _sample!(onhit::Function, onmiss::Function, mcm::MC0Model{T}, rng, niters) where {T}
    
    x::Vector{T} = mcm.x
    vf::Vector{T} = mcm.vf
    enet::EchelonMetNet = mcm.enet
    lb::Vector{T} = enet.net.lb
    ub::Vector{T} = enet.net.ub

    vf0::Vector{T} = lb[enet.idxf]
    Δvf::Vector{T} = ub[enet.idxf] - vf0

    it = 0
    while true

        it += 1
        
        for i in eachindex(vf)
            vf[i] = vf0[i] + rand(rng) * Δvf[i]
        end
        
        isfea = isfeasible_vf!(x, enet, vf; testfree = false)
        flag = isfea ? onhit(mcm) : onmiss(mcm)
        flag === true && break

        it >= niters &&  break
    end

    return nothing
end

function _sample!(mcm::MC0Model, nsamples::Int, idxs::AbstractVector;
        niters = Inf,
        rng = mcm.rng,
    )
    si = 0
    samples = zeros(nsamples, length(idxs))
    _sample!(_do_nothing, mcm, rng, niters) do _mcm
        si += 1
        samples[si, idxs] .= _mcm.x[idxs]
        return nsamples == si
    end
    return samples
end

_sample!(mcm::MC0Model, nsamples::Int, idxs::Nothing; kwargs...) = 
    _sample!(mcm::MC0Model, nsamples::Int, eachindex(mcm.enet.net.rxns); kwargs...)

function _sample!(mcm::MC0Model, nsamples::Int, idx::Int;
        niters = Inf,
        rng = mcm.rng,
    )
    si = 0
    samples = zeros(nsamples)
    _sample!(_do_nothing, mcm, rng, niters) do _mcm
        si += 1
        samples[si] = _mcm.x[idx]
        return nsamples == si
    end
    return samples
end


_do_nothing(x...) = nothing

## ------------------------------------------------------------------
export sample!
import Distributions.sample!
function sample!(onhit::Function, mcm::MC0Model, niters; 
        onmiss::Function = _do_nothing,
        rng = mcm.rng
    ) 
    _sample!(onhit, onmiss, mcm, rng, niters)
end

function sample!(mcm::MC0Model, nsamples::Int; 
        niters = Inf,
        rng = mcm.rng,
        rxns = nothing
    ) 
    idxs = isnothing(rxns) ? nothing : rxnindex(mcm.enet.net, rxns)
    return _sample!(mcm, nsamples, idxs; rng, niters)
end