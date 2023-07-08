import MetXBase.span!
export span!

function span!(mcm::AbstractHitOrDropSampler) 
    mcm.v[mcm.idxd] .= mcm.be - mcm.G * mcm.vi
    mcm.v[mcm.idxi] .= mcm.vi
    return mcm.v
end
