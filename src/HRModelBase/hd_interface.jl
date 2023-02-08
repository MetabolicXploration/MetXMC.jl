function _hit_or_drop!(mvn::HRModel, ::Bool)
    _walk!(mvn, mvn.rng)
    return viwdf(mvn)
end