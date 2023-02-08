export viwdf
export wdf


export vipoint, vipoint!
vipoint(m::AbstractHitOrDropSampler) = m.vi
vipoint!(m::AbstractHitOrDropSampler, vi) = (m.vi .= vi)
export vpoint, vpoint!
vpoint(m::AbstractHitOrDropSampler) = m.v
vpoint!(m::AbstractHitOrDropSampler, v) = (m.v .= v)