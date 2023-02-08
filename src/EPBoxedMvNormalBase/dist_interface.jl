# weight of a independent point
function viwdf(mvn::EPBoxedMvNormal{T}) where T
    for i in eachindex(mvn.v)
        v = mvn.v[i]
        v < mvn.lb[i] && return zero(T)
        v > mvn.ub[i] && return zero(T)
    end
    return pdf(mvn.Qi, mvn.vi)
end

# weight of a full point
wdf(mvn::EPBoxedMvNormal) = viwdf(mvn)