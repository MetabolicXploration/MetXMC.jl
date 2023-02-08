# ## ------------------------------------------------------------------
# function _sample!(
#         onhit::Function, onmiss::Function, 
#         mcm::MC0Model{T}, G, rng, niters;
#         testfree = true
#     ) where {T}
    
#     x::Vector{T} = mcm.x
#     vf::Vector{T} = mcm.vf
#     enet::EchelonMetNet = mcm.enet

#     it = 0
#     while true
#         it += 1
#         vf .= rand(rng, G)
        
#         isfea = isfeasible_vf!(x, enet, vf; testfree)
#         flag = isfea ? onhit(mcm) : onmiss(mcm)
#         flag === true && break

#         it >= niters &&  break
#     end

#     return nothing
# end

# function _sample!(
#         mcm::MC0Model, G, nsamples::Int, idxs::AbstractVector;
#         niters = Inf,
#         rng = mcm.rng,
#         testfree = true
#     )
#     si = 0
#     ws = zeros(nsamples)
#     samples = zeros(nsamples, length(idxs))
#     _sample!(_do_nothing, mcm, G, rng, niters; testfree) do _mcm
#         si += 1
#         samples[si, idxs] .= _mcm.x[idxs]
#         ws[si] = pdf(G, _mcm.vf)
#         return nsamples == si
#     end
#     return ws, samples
# end

# _sample!(mcm::MC0Model, G, nsamples::Int, idxs::Nothing; kwargs...) = 
#     _sample!(mcm, G, nsamples, eachindex(mcm.enet.net.rxns); kwargs...)

# function _sample!(mcm::MC0Model, G, nsamples::Int, idx::Int;
#         niters = Inf,
#         rng = mcm.rng,
#         testfree = true
#     )
#     si = 0
#     ws, samples = zeros(nsamples), zeros(nsamples)
#     _sample!(_do_nothing, mcm, G, rng, niters; testfree) do _mcm
#         si += 1
#         samples[si] = _mcm.x[idx]
#         ws[si] = pdf(G, _mcm.vf)
#         return nsamples == si
#     end
#     return ws, samples
# end

# _do_nothing(x...) = nothing

# ## ------------------------------------------------------------------
# export sample!
# import Distributions.sample!
# function sample!(onhit::Function, mcm::MC0Model, niters; 
#         G = mcm.U,
#         rng = mcm.rng,
#         testfree = true,
#         onmiss::Function = _do_nothing
#     ) 
#     _sample!(onhit, onmiss, mcm, G, rng, niters; testfree)
# end

# function sample!(mcm::MC0Model, nsamples::Int;
#         niters = Inf,
#         G = mcm.U,
#         rng = mcm.rng,
#         rxns = nothing, 
#         testfree = true,
#     ) 
#     idxs = isnothing(rxns) ? nothing : rxnindex(mcm.enet.net, rxns)
#     return _sample!(mcm, G, nsamples, idxs; rng, niters, testfree)
# end