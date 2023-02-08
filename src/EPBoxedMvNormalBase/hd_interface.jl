function _hit_or_drop!(mvn::EPBoxedMvNormal, ::Bool)
    mvn.vi .= rand(mvn.rng, mvn.Qi)
    span!(mvn)
    return viwdf(mvn)
end