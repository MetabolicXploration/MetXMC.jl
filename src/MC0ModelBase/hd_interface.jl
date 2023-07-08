function _hit_or_drop!(hrm::MC0Model{T}, ::Bool) where {T}
    hrm.vi .= rand(hrm.rng, hrm.U)
    isfea = isfeasible_vf!(hrm.v, hrm.elep, hrm.vi; testfree = false)
    return isfea ? viwdf(hrm) : zero(T)
end