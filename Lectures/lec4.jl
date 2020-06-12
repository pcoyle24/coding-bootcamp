# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 06/11/2020
# lec4.jl
# ------------------------------------------------------------------------------

using Optim, Interpolations, Plots, Parameters, Polynomials

## Polynomial Approximations
# Runge's Problem (why poly not always a good idea)
Runge(x) = 1/(1 + 25*x^2)

x_fine = collect(-1.:0.01:1.)
plot(x_fine, Runge.(x_fine))

# Fit 5th degree polynomial
x_coarse = collect(range(-1.0, length = 6, stop = 1.0))
poly_fit_5 = fit(x_coarse, Runge.(x_coarse), 5)
plot!(x_fine, poly_fit_5.(x_fine))

# Higher order Poly does not solve the problem
x_coarse_2 = collect(range(-1.0, length = 10, stop = 1.0))
poly_fit_9 = fit(x_coarse_2, Runge.(x_coarse_2), 9)
plot!(x_fine, poly_fit_9.(x_fine))

## Linear/Spine Interpolation (both vastly out perform polynomial fit)
x_range = -1.0:0.25:1.0

# linear interp
runge_linear = LinearInterpolation(x_range, Runge.(x_range))

# cubic spline
runge_cublic = CubicSplineInterpolation(x_range, Runge.(x_range))

plot(x_fine, [Runge.(x_fine), runge_linear.(x_fine), runge_cublic.(x_fine)])

## Other way to define interps
# Sine wave
x = collect(1.:1:10.)
x_fine = collect(1.:0.01:10.)

y = sin.(x_coarse)
y_fine = sin.(x_fine)

plot(x_fine, y_fine)
scatter!(x,y, markersize = 2) #(! is grid hold)

# Creating intepolation
y_linear = interpolate(y, BSpline(Linear()))
y_cubic = interpolate(y, BSpline(Cubic(Line(OnGrid())))) #Must specify boundary conditions to make cubic work
plot(x_fine,[y_fine, y_linear.(x_fine), y_cubic.(x_fine)])
scatter!(x,y, markersize = 2) #(! is grid hold)
# In curvy region of function, linear really really struggles

#Some Julia Optimization: Only used discritized *range* of function in interpolate command

grid_simple = collect(0:2:20)
interp = interpolate(grid_simple, BSpline(Linear())) #expects index as input
interp(0)
interp(20)
interp(1) #Expects an INDEX
interp(11)

interp(1.5) #half way between first and second index
# To evaluate at *value* must first figure out it's *index* position

## VFI Revisted
#struct to hold model primitives
@with_kw struct Primitives
    β::Float64 = 0.99 #discount factor
    θ::Float64 = 0.36 #production
    δ::Float64 = 0.025 #depreciation
    k_grid::Array{Float64,1} = collect(range(0.1, length = 101, stop= 45.0)) #capital grid is v coarse... under this code, VF looks like crap
    nk::Int64 = length(k_grid) #number of capital grid states
end

mutable struct Results
    val_func::Array{Float64,1} #value function
    pol_func::Array{Float64,1} #policy function
end

#function to solve the model
function Solve_model()
    #initialize primitives and results
    prim = Primitives()
    val_func, pol_func = zeros(prim.nk), zeros(prim.nk)
    res = Results(val_func, pol_func)


    error, n = 100, 0
    while error>0.0001 #loop until convergence
        n+=1
        v_next = Bellman(prim, res) #next guess of value function
        error = maximum(abs.(v_next .- res.val_func)) #check for convergence
        res.val_func = v_next #update
        # println("Current error: ", error)
    end
    println("Value function converged in ", n, " iterations")
    prim, res
end

#Bellman operator
function Bellman(prim::Primitives, res::Results)
    @unpack β, δ, θ, nk, k_grid = prim #unpack primitive structure
    v_next = zeros(nk) #preallocate next guess of value function
    k_interp = interpolate(k_grid,BSpline(Linear()))
    val_interp = interpolate(res.val_func,BSpline(Linear()))
    # val_interp = interpolate(res.val_func,(BSpline(Linear()), NoInterp)) #Only interpolating over one dimension (for stoch shock)

    for i_k = 1:nk #loop over state space
        k = k_grid[i_k] #value of capital
        budget = k^θ + (1-δ)*k #budget

        # TOTALLY DITCH GRID SEARCH
        # for i_kp = 1:nk
        #     kp = k_grid[i_kp] #value of k'
        #     c = budget - kp #consumption
        #     if c>0 #check if postiive
        #         val = log(c) + β * res.val_func[i_kp] #compute utility
        #         if val > max_util #wow new max!
        #             max_util = val #reset max value
        #             res.pol_func[i_k] = kp #update policy function
        #         end
        #     end
        # end

        # Define function (to be OPTIMIZED over)
        val(i_kp) = log(budget - k_interp(i_kp)) + β*val_interp(i_kp)
        obj(i_kp) = -val(i_kp)

        # optimize via brents method
        lower = 1.0
        upper = get_index(budget, k_grid)
        opt = optimize(obj, lower, upper)
        res.pol_func[i_k] = k_interp(opt.minimizer)
        v_next[i_k] = -opt.minimum #update value function
    end
    return v_next
end

function get_index(val::Float64, grid::Array{Float64, 1})
    n = length(grid)
    index = 0
    if val <= grid[1]
        index = 1
    elseif val >= grid[n]
        index = n
    else
        index_upper = findfirst(x->x>val, grid)
        index_lower = index_upper - 1
        val_upper, val_lower = grid[index_upper], grid[index_lower]

        index = index_lower + (val - val_lower) / (val_upper - val_lower)
    end
    return index
end





@elapsed prim, res = Solve_model() #solve the model.
plot(prim.k_grid, res.val_func)
plot(prim.k_grid, res.pol_func)
