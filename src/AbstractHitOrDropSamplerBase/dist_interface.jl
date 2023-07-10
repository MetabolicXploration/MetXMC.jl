export viwdf
export wdf

vipoint(m::AbstractHitOrDropSampler) = m.vi
vipoint!(m::AbstractHitOrDropSampler, vi) = (m.vi .= vi)

vpoint(m::AbstractHitOrDropSampler) = m.v
vpoint!(m::AbstractHitOrDropSampler, v) = (m.v .= v)