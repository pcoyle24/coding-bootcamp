####parallelization example
using Distributed

#add processes
workers()
addprocs(5)

@everywhere using Distributions, Statistics, SharedArrays

@everywhere function draw(σ::Int64, n::Int64)
    dist = Normal(0, σ)
    draws = rand(dist, n)
    v1, v2 = mean(draws), std(draws)
    v1, v2
end

temp = zeros(10,2)
@elapsed for i = 1:10
    mean, sd = draw(i, 100000000)
    temp[i,1] = mean
    temp[i,2] = sd
end

temp = SharedArray{Float64}(10,2)
@elapsed @sync @distributed for i = 1:10
    mean, sd = draw(i, 100000000)
    temp[i,1] = mean
    temp[i,2] = sd
end




#
