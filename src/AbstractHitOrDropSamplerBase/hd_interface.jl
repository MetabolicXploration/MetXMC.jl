# _hit_or_drop!(s::AbstractHitOrDropSampler, span::Bool)
# This is the main interface method of AbstractHitOrDropSampler
# sample! the Sampler and return a weight.
# If the weight == 0.0 the sample shoulf be droped
_hit_or_drop!(m::AbstractHitOrDropSampler, ::Bool) = error("_hit_or_drop!($(typeof(m)), ::Bool) not implemented")

## ------------------------------------------------------------------
_do_nothing(x...) = nothing

## ------------------------------------------------------------------
# sample interface
function _sample!(onhit::Function, ondrop::Function, rw::Function,
        mcm::AbstractHitOrDropSampler, niters, dospan::Bool
    )

    it = 0
    while true
        it += 1
        w0 = _hit_or_drop!(mcm, dospan)

        flag = iszero(w0) ? ondrop(rw(w0)) : onhit(rw(w0))
        flag === true && break

        it >= niters &&  break
    end

    return nothing
end

function _sample!(
        mcm::AbstractHitOrDropSampler, nsamples::Int, idxs::AbstractVector;
        niters = Inf,
        dospan = true, 
        rw = identity
    )
    si = 0
    ws = zeros(nsamples)
    samples = zeros(nsamples, length(idxs))
    _sample!(_do_nothing, rw, mcm, niters, dospan) do w
        si += 1
        samples[si, idxs] .= vpoint(mcm)[idxs]
        ws[si] = w
        return nsamples == si
    end
    return ws, samples
end

_sample!(mcm::AbstractHitOrDropSampler, nsamples::Int, idxs::Nothing; kwargs...) = 
    _sample!(mcm, nsamples, eachindex(vpoint(mcm)); kwargs...)

function _sample!(mcm::AbstractHitOrDropSampler, nsamples::Int, idx::Int;
        niters = Inf,
        dospan = true, 
        rw = identity
    )
    si = 0
    ws, samples = zeros(nsamples), zeros(nsamples)
    _sample!(_do_nothing, rw, mcm, niters, dospan) do w
        si += 1
        samples[si] = vpoint(mcm)[idx]
        ws[si] = w
        return nsamples == si
    end
    return ws, samples
end

## ------------------------------------------------------------------
# TODO: Think about saparating 'wsample -> ws, samples' and 'sample -> samples'
import Distributions.sample!
function sample!(onhit::Function, mcm::AbstractHitOrDropSampler, niters; 
        ondrop::Function = _do_nothing, 
        dospan = true, 
        rw = identity
    ) 
    _sample!(onhit, ondrop, rw, mcm, niters, dospan)
end

function sample!(mcm::AbstractHitOrDropSampler, nsamples::Int;
        niters = Inf,
        rxns = nothing, 
        dospan = true, 
        rw = identity
    ) 
    idxs = isnothing(rxns) ? nothing : rxnindex(mcm, rxns)
    return _sample!(mcm, nsamples, idxs; rw, niters, dospan)
end

## ------------------------------------------------------------------
function sample_histogram!(
        mvn::AbstractHitOrDropSampler, rxni::Int, 
        bins::AbstractVector, hist::AbstractVector; 
        nsamples,
        sample_kwargs...
    )

    sample!(mvn, nsamples; sample_kwargs...) do w
        vi = vpoint(mvn)[rxni]
        _histogram!(bins, hist, vi, w)
    end
    return (bins, hist)
end

function sample_histogram!(
        mvn::AbstractHitOrDropSampler, rxni::Int, 
        bins::AbstractVector;
        kwargs...
    )

    hist = zeros(length(bins))
    return sample_histogram!(mvn, rxni, bins, hist; kwargs...)
end

function sample_histogram!(
        mvn::AbstractHitOrDropSampler, rxni::Int, 
        v0::Real, v1::Real, nbins::Int; 
        kwargs...
    )

    bins = range(v0, v1; length = nbins)
    hist = zeros(length(bins))
    return sample_histogram!(mvn, rxni, bins, hist; kwargs...)
end
