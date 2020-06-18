#model primitives
@with_kw struct Primitives
    β::Float64 = 0.95 #discount rate
    ρ::Float64 = 2.0 #CRRA coefficient
    κ::Float64 = 0.7 #Ben-Porath HC coefficient

    #distribution of HC shocks
    μ::Float64 = -0.029
    σ::Float64 = 0.111

    #grids and dimensions
    a_grid::Vector{Float64} = collect(range(0.01, length = 10, stop = 0.6)) #ability grid
    na = length(a_grid) #number of ability points
    h_grid::Vector{Float64} = collect(range(50, length = 50, stop = 200)) #HC grid
    nh = length(h_grid)
    age_grid = collect(20:1:70)
    n_age = length(age_grid)
end

#model results
mutable struct Results
    val_func::Array{Float64,3} #value function
    pol_func::Array{Float64,3} #value function
end

#equal-mass discretization of continuous dsitribution method from Kennan (2006)
function discretize_lognormal(μ::Float64, σ::Float64, n::Int64)
    dist = LogNormal(μ, σ) #set up log normal distribution
    shocks = zeros(n) #preallocate market luck states

    for i = 1:n #begin loop to fill discretized vector
        quant = (2*i-1)/(2*n) #desired quantile
        shocks[i] = quantile(dist, quant) #fill element of state vector
    end
    shocks #return deliverable
end

###function for obtaining index of arbitary value in HC grid
function get_index(val::Float64, grid::Array{Float64,1})
    n = length(grid)
    index = 0 #preallocation
    if val<=grid[1] #LEQ smallest element
        index = 1
    elseif val>=grid[n] #GEQ biggest element
        index = n
    else
        index_upper = findfirst(z->z>val, grid)
        index_lower = index_upper - 1
        val_upper, val_lower = grid[index_upper], grid[index_lower] #values
        index = index_lower + (val - val_lower)  / (val_upper - val_lower) #weighted average
    end
    index #return
end

##CRRA utility
function utility(c::Float64, ρ::Float64)
    return (c^(1-ρ))/(1-ρ)
end
############
