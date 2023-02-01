using MetXMC
using MetXOptim
using MetXBase
using MetXBase: _isapprox, _histogram
using MetXNetHub
using GLPK
using Test

@testset "MetXMC.jl" begin
    
    # MC comparizon
    let 
        model_id = "toy_net4D"
        net0 = pull_net(model_id)
        net = box(net0, GLPK.Optimizer)
        enet = EchelonMetNet(net; tol = 1e-10)

        global data_pool = Dict()
        
        for (mcm_id, mcm) in [
                ("MC0Model", MC0Model(enet)), 
                ("HRModel", HRModel(net, GLPK.Optimizer)), 
            ]

            println("="^60)
            @show mcm_id

            # samples
            niters = 500000
            @time samples = sample!(mcm, niters)
            @show size(samples)
            
            for id0 in net.rxns

                id0i = rxnindex(net, id0)

                dists = get!(data_pool, id0, [])

                # marginals
                bins, hist = _histogram(samples[:, id0i], 100)
                hist ./= sum(hist)
                
                push!(dists, hist)

            end # for i in 1:N
        
        end # for (mcm_id, mcm)

        for id0 in net.rxns
            dists = data_pool[id0]
            @test _isapprox(dists...; rtol = 5e-2)
        end

    end
end
