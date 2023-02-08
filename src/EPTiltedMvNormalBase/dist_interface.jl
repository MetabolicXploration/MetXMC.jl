# weight of a independent point
function viwdf(mvn::EPTiltedMvNormal{T}) where T
    v = mvn.v[mvn.rxni]
    v < mvn.lb && return zero(T)
    v > mvn.ub && return zero(T)
    return pdf(mvn.Qi, mvn.vi)
end

# weight of a full point
wdf(mvn::EPTiltedMvNormal) = viwdf(mvn)