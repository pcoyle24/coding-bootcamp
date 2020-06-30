# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 06/25/2020
# lec8.jl
# ------------------------------------------------------------------------------

## Birthday Problem
using Random, Plots, Distributions, Statistics, Parameters

function birthday(n::Int64, sims::Int64)
    results = zeros(sims)

    for i = 1:sims
        days = rand(1:365, n)
        results[i] = length(unique(days)) #returns unique elements passed to it
    end
    results
end

res_20, res_50, res_70 = birthday(20, 10000), birthday(50, 10000), birthday(70, 10000)

h20 = histogram(res_20, label = "n = 20", normed=true);
h50 = histogram(res_50, label = "n = 50", normed=true);
h70 = histogram(res_70, label = "n = 70", normed=true);

plot(h20, h50, h70, layout=(1,3))

## Distance of points in a cube
function Point_Distance(sim::Int64)
    results = zeros(sim)
    for i = 1:sim
        p1, p2 = rand(3), rand(3)
        results[i] = sqrt(sum((p1.-p2).^2))
    end
    return mean(results)
end

Point_Distance(100000)

## Labor Supply Problem
@with_kw struct Primitives
    β_0::Float64 = 2.7
    β_1::Float64 = 0.47
    σ::Float64 = 0.597
    α::Float64 = 1.0
    B::Float64 = 5.0
end

mutable struct Results
    emax::Array{Float64,1}
    lfp::Array{Float64,1}
end

function Solve_model(sims::Int64)
    prim = Primitives()
    res = Results(zeros(2), zeros(2))

    compute_emax(prim, res, sims)

    return prim, res
end

function compute_emax(prim::Primitives, res:: Results, sims::Int64)
    @unpack β_0, β_1, σ, α, B = prim
    dist = Normal(0,σ)
    val_nwork = α + log(B)
    utils, lfps = zeros(2, sims), zeros(2, sims)

    for s = 0:1
        for i = 1:sims
            ε = rand(dist)
            wage = exp(β_0 + β_1*s + ε)
            util = max(val_nwork, log(wage))
            utils[s+1, i] = util
            lfps[s+1, i] = (log(wage)>val_nwork)
        end
    end


    res.emax[1], res.emax[2] = mean(utils[1,:]), mean(utils[2,:])
    res.lfp[1], res.lfp[2] = mean(lfps[1,:]), mean(lfps[2,:])
end

prim, res = Solve_model(100)
