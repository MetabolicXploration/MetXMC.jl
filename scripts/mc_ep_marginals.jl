@time begin
    using MetXMC
    using MetXNetHub
    using MetXOptim
    using MetXOptim: GLPK
    using MetXBase
    using MetXEP
    using Plots
    using Distributions
end

# TODO: Positive define net (struct)

## ------------------------------------------------------------------
let
    fn = "/Users/Pereiro/Downloads/HEK293.mat"
    global net = MetXBase.load_net(fn)
    opm = fba(net, GLPK.Optimizer)
end

## ------------------------------------------------------------------
let
    global net = MetXNetHub.pull_net("iJO1366")
    global net = box(net, GLPK.Optimizer)
    lb0, ub0 = fva(net, GLPK.Optimizer; verbose = true)
    # sort(ub0 - lb0)
    # met_rxns(net, 1)
    
    # nothing
end

## ------------------------------------------------------------------
# TODO: Support maxent sampling exp(βx)
## ------------------------------------------------------------------
# Utils
# TODO: This is wrong!!! Fix it on MetXPlots
function _plot_marginal!(p::AbstractPlot, epm::FluxEPModelT0, rxn; 
        pkwargs...
    )

    rxn = rxnindex(epm, rxn)

    μ = [epm.μd; epm.μi][epm.idxmap_inv[rxn]]
    σ = [epm.sd; epm.si][epm.idxmap_inv[rxn]]
    lb = [epm.lbd; epm.lbi][epm.idxmap_inv[rxn]]
    ub = [epm.ubd; epm.ubi][epm.idxmap_inv[rxn]]

    xs, ys = MetXEP.sample_tnorm(
        μ * epm.scalefact, 
        σ * epm.scalefact^2, 
        lb * epm.scalefact, 
        ub * epm.scalefact; 
        xbins = 1000, digits = 15
    )
    plot!(p, xs, ys; pkwargs...)
end

function _marginal(epm::FluxEPModelT0, rxn)

    rxn = rxnindex(epm, rxn)

    μ = [epm.μd; epm.μi][epm.idxmap_inv[rxn]] * epm.scalefact
    σ = [epm.sd; epm.si][epm.idxmap_inv[rxn]] * epm.scalefact^2
    lb = [epm.lbd; epm.lbi][epm.idxmap_inv[rxn]] * epm.scalefact
    ub = [epm.ubd; epm.ubi][epm.idxmap_inv[rxn]] * epm.scalefact

    return TruncatedNormal(μ, sqrt(σ), lb, ub)
end

## ------------------------------------------------------------------
let
    # global net = MetXNetHub.pull_net("ecoli_core")
    global net = MetXNetHub.pull_net("toy_net")
    global net = box(net, GLPK.Optimizer)
    global mcm = MetHRModel(net, GLPK.Optimizer)
    global epm = FluxEPModelT0(net)
    config!(epm, :epsconv, 1e-7)
    converge!(epm)
end

## ------------------------------------------------------------------
# TODO: sample the epm distribution
let
    rxn = rand(reactions(net))
    nsamples = Int(3e5)
    
    p = plot(; title = rxn)

    # HR
    # p = plot()
    hr_samples = MetXMC.sample!(mcm; nsamples, rxns = rxn)
    histogram!(p, hr_samples;
        label = "MC",
        alpha = 0.6, 
        normalize = :pdf,
        c = :black
    )

    # EP
    _plot_marginal!(p, epm, rxn; 
        title = rxn, 
        label = "EP",
        alpha = 0.6, 
        lw = 5
    )
    
    p
end

## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
let
    # rxns = rand(reactions(net))
    # rxns = "EX_o2_e"
    rxns = "BIOMASS_Ecoli_core_w_GAM"
    nsamples = Int(8e3)
    # @time samples = sample!(mcm; rxns, nsamples)

    # MaxEnt
    β = 1e3
    l, u = bounds(net, rxns)
    # TODO: make numerically stable
    
    samples = zeros(nsamples)
    c = 1
    while c <= nsamples
        x = sample!(mcm; rxns)
        # P(x, β) < rand() && continue
        samples[c] = x
        c += 1
    end

    

    # Plots
    # histogram(samples; 
    #     label = reactions(net, rxns),
    #     xlabel = "flx",
    #     ylabel = "pdf",
    #     normalize = :pdf,
    #     c = :black
    # )
end

# Solve this

## ------------------------------------------------------------------
## ------------------------------------------------------------------

let
    net = MetXNetHub.pull_net("ecoli_core")
    S, b, lb, ub = net.S, net.b, net.lb, net.ub
    x0 = hrsample(S, b, lb, ub; nsamples = 2000)
    histogram(x0[:, rand(eachindex(lb))]; label = "", c = :black)
end



## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
let
    net = MetXNetHub.pull_net("toy_net")
    S, b, lb, ub = net.S, net.b, net.lb, net.ub
    x0 = warmup(S, b, lb, ub)
    println.(lb, "\t-\t", x0, "\t-\t", ub)
    return nothing
end


## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
function P(x, β, l, u) 
    iszero(β) && return inv(u - l)
    return exp(β * x) / ((exp(β * u) - exp(β * l)) / β )
end

let
    β = 1e1
    l, u = 0.0, 1.0
    xs = range(l, u; length = 1000)
    plot(xs, P.(xs, β, l, u))
end