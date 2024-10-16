using MetXMC
using MetXOptim
using MetXBase
using MetXBase: _isapprox, _histogram
using MetXNetHub
using GLPK
using Test
using Random

@testset "MetXMC.jl" begin

    Random.seed!(123)
    
    # MC comparizon
    let 
        model_id = "toy_net4D"
        net0 = pull_net(model_id)
        lep0 = lepmodel(net0)
        lep = fva_strip(lep0, GLPK.Optimizer)
        elep = EchelonLEPModel(lep; tol = 1e-10)

        data_pool = Dict()
        
        for (mcm_id, mcm) in [
                ("MC0Model", MC0Model(elep)), 
                ("HRModel", HRModel(lep, GLPK.Optimizer)), 
            ]

            println("="^60)
            @show mcm_id

            # samples
            niters = 500_000
            @time ws, samples = sample!(mcm, niters)
            @show size(samples)
            
            for id0 in lep.colids

                id0i = colindex(lep, id0)

                dists = get!(data_pool, id0, [])

                # marginals
                bins, hist = _histogram(samples[:, id0i], 100)
                hist ./= sum(hist)
                
                push!(dists, hist)

            end # for i in 1:N
        
        end # for (mcm_id, mcm)

        for (_, dists) in data_pool
            Ss = [sum(p .+ log.(p)) for p in dists]
            @test _isapprox(Ss...; rtol = 1e-2)
        end

    end
end
