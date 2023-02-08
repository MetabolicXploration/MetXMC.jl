# weight of a independent point
viwdf(mvn::EPJoinMvNormal, vi::AbstractVector) = pdf(mvn.Qi, vi)
viwdf(mvn::EPJoinMvNormal) = viwdf(mvn, mvn.vi)

# weight of a full point
wdf(mvn::EPJoinMvNormal, v) = viwdf(mvn, v[mvn.epm.idxi])
wdf(mvn::EPJoinMvNormal) = viwdf(mvn, mvn.v)