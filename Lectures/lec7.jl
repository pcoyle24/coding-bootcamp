# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 06/23/2020
# lec7.jl
# ------------------------------------------------------------------------------

using Distributed #used to distribute across cores

## add processes
# How many cores currenty using?
workers()

# use more cores (use all 6 of available)
addprocs(5)

# make sure all cores use packages (everywhere macro)
@everywhere using Statistics, SharedArrays, Distributions

@everywhere function draw(σ::Int64, n::Int64)
    dist = Normal(0, σ)
    draws = rand(dist,n)
    v1, v2 = mean(draws), std(draws)
    return v1, v2
end


# Run on single core
temp = zeros(10,2)
@elapsed for i = 1:10
    mean, sd = draw(i, 100000000)
    temp[i, 1] = mean
    temp[i, 2] = sd
end

# Run on all cores
temp = SharedArray{Float64}(10,2)
@elapsed @sync @distributed for i = 1:10 #dist macro is parfor equiv, sync is
    mean, sd = draw(i, 100000000)
    temp[i, 1] = mean
    temp[i, 2] = sd
end

#Overhead for small tasks can mean that single core runs faster.
