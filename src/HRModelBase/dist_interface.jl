# weight of a independent point
viwdf(::HRModel{T}) where T = one(T)

# weight of a full point
wdf(mvn::HRModel) = viwdf(mvn)