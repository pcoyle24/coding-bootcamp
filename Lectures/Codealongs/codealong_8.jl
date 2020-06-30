using Random, Plots, Distributions, Statistics, Parameters

######birthday problem
function birthday(n::Int64, sims::Int64)
    results = zeros(sims) #preallocate monte-carlo results vector
    for i = 1:sims #loop over simulations
        days = rand(1:365, n) #draw birthdays
        results[i] = length(unique(days)) #fill in vector
    end
    results #return
end
res_20, res_50, res_70 = birthday(20, 10000), birthday(50, 10000), birthday(70, 10000)
Plots.histogram(res_20)
Plots.histogram(res_50)
Plots.histogram(res_70)

#####Average distance between two random points in a cube
function Point_distance(sims::Int64)
    results = zeros(sims)
    for i = 1:sims #loop over simulations
        p1, p2 = rand(3), rand(3) #two points!
        results[i] = sqrt(sum((p1.-p2).^2))
    end
    return mean(results)
end
Point_distance(10000)

####something more involved: computing expected value of college factoring in wage offer shocks
@with_kw struct Primitives
    β_0::Float64 = 2.7 #wage constant
    β_1::Float64 = 0.47 #college premium
    σ::Float64 = 0.597 #wage offer SD
    α::Float64 = 1.0 #leisure
    B::Float64 = 5.0 #base consumption
end

mutable struct Results
    emax::Array{Float64,1} #first for no college, second for college
    lfp::Array{Float64,1} #lfp probabilities
end

#initializes model primitives and executes solution
function Solve_model(sims::Int64)
    prim = Primitives() #initialize primitives
    res = Results(zeros(2), zeros(2)) #initialize resutls
    Compute_emax(prim, res, sims) #solve model
    prim, res #return
end

function Compute_emax(prim::Primitives, res::Results, sims::Int64)
    #housekeeping
    @unpack β_0, β_1, σ, α, B = prim
    dist = Normal(0, σ)
    val_nwork = α + log(B)
    utils, lfps = zeros(2, sims), zeros(2, sims)

    for s = 0:1 #loop over schooling levels
        for i = 1:sims #loop over simulations
            ε = rand(dist) #draw shock and compute resultant wage
            wage = exp(β_0 + β_1*s + ε)
            util = max(log(wage), val_nwork) #max utility
            utils[s+1, i] = util #update
            lfps[s+1, i] = (log(wage)>val_nwork) #update working decision
        end
    end

    res.emax[1], res.emax[2] = mean(utils[1,:]), mean(utils[2,:]) #expected max utilities
    res.lfp[1], res.lfp[2] = mean(lfps[1,:]), mean(lfps[2,:]) #LFP probabilities
end

prim, res = Solve_model(100) #run the code


#
