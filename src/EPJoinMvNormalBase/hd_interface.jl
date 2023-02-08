function _hit_or_drop!(mvn::EPJoinMvNormal, ::Bool)
    mvn.vi .= rand(mvn.rng, mvn.Qi)
    span!(mvn)
    return viwdf(mvn)
end