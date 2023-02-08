function _hit_or_drop!(hrm::MC0Model{T}, ::Bool) where {T}
    hrm.vf .= rand(rng, hrm.U)
    isfea = isfeasible_vf!(hrm.v, hrm.enet, hrm.vf; testfree = false)
    return isfea ? viwdf(hrm) : zero(T)
end