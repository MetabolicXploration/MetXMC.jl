# weight of a independent point
viwdf(::MC0Model{T}) where T = one(T)

# weight of a full point
wdf(mvn::MC0Model) = viwdf(mvn)