using Optim, Interpolations, Plots, Polynomials, Parameters

###Runge's function/phenomenon##
Runge(x) = 1 /(1 + 25*x^2)
x_fine = collect(-1.0:0.01:1.0)
x_coarse = collect(range(-1.0, length = 6, stop = 1.0))
Plots.plot(x_fine, Runge.(x_fine))

#5-degree polynomial fit
poly_fit = fit(x_coarse, Runge.(x_coarse), 5)
Plots.plot!(x_fine, poly_fit.(x_fine))

#try a higher-degree -- surely that will help!
x_coarse_2 = collect(range(-1.0, length = 10, stop = 1.0))
poly_fit_2 = fit(x_coarse_2, Runge.(x_coarse_2), 9)
Plots.plot!(x_fine, poly_fit_2.(x_fine))

#OK, how about linear/spline interpolation
x_range = -1.0:0.25:1.0
runge_range = Runge.(x_range)
runge_linear = LinearInterpolation(x_range, runge_range)
runge_spline = CubicSplineInterpolation(x_range, runge_range)
Plots.plot(x_fine, [runge_grid, runge_linear.(x_fine), runge_spline.(x_fine)])

####Interpolating a sine wave
x = collect(1:1:10)
y = sin.(x)
x_fine = collect(1:0.01:10)
y_fine = sin.(x_fine)

#plot
plot(x_fine, y_fine)
scatter!(x, y, markersize = 4)

#either linear or cubic splines
interp_y = interpolate(y, BSpline(Linear()))
interp_y = interpolate(y, BSpline(Cubic(Line(OnGrid()))))
y_fine_interp = interp_y.(x_fine)
plot(x_fine, [y_fine, y_fine_interp])
scatter!(x, y, markersize = 4)

###Interpolate -- here we learn about the indexing
interp = interpolate(collect(0:2:20), BSpline(Linear()))

####Bilinear interpolation
f(x, y) = x^2 + y^2
grid_coarse = collect(0.0:0.5:5.0)
grid_fine = collect(0.0:0.01:5.0)
ncoarse, nfine = length(grid_coarse), length(grid_fine)
z_grid = zeros(nfine, nfine)

for i = 1:nfine, j = 1:nfine
    x, y = grid_fine[i], grid_fine[j]
    z_grid[i,j] = f(x, y)
end

Plots.contourf(grid_fine, grid_fine, z_grid)

#now get to bilinear inteerp
func_grid_coarse = zeros(ncoarse, ncoarse)
for i = 1:ncoarse, j = 1:ncoarse
    x, y = grid_coarse[i], grid_coarse[j]
    func_grid_coarse[i,j] = f(x, y)
end

#create the interpolation
f_grid_interp = interpolate(func_grid_coarse, BSpline(Linear()))
z_interp = zeros(nfine, nfine)
for i = 1:nfine, j = 1:nfine
    x, y = grid_fine[i], grid_fine[j]

    #note: we need to convert these values of x and y to indices in the coarsened grid
    #we can do this by noting the following mapping:
    #in the coarsened grid, zero is the first element (0->1)
    #the second element is 0.5 (0.5->2)
    #the third is 1 and so on (1->3)
    #the rule, then is index(x) = 2x+1

    x_ind, y_ind = 2x+1, 2y+1 #conversion to index
    z_interp[i,j] = f_grid_interp(x_ind, y_ind)
end

#compare contour plots
Plots.contourf(grid_fine, grid_fine, z_grid)
Plots.contourf(grid_fine, grid_fine, z_interp) #slightly more jagged, but otherwise pretty good!

#####Better version of optimal growth
######Optimal Growth with Interpolation / Optimization

using Parameters

#a struct to hold model primitives.
@with_kw struct Primitives
    β::Float64 = 0.99 #discount rate.
    θ::Float64 = 0.36 #capital share
    δ::Float64 = 0.025 #capital depreciation
    k_grid::Array{Float64,1} = collect(range(0.1, length = 50, stop = 45.0)) #capital grid. Much more coarse.
    nk::Int64 = length(k_grid) #number of capital elements
end

#mutable struct to hold model results
mutable struct Results
    val_func::Array{Float64,1}
    pol_func::Array{Float64,1}
end

#function that executes the model and returns results
function Solve_model()
    prim = Primitives() #initialize primitives
    val_func = zeros(prim.nk) #preallocate value function as a vector of zeros
    pol_func = zeros(prim.nk) #preallocate value function as a vector of zeros
    res = Results(val_func, pol_func) #initialize results
    V_iterate(prim, res) #value function iteration
    prim, res #return deliverables
end

#value function iteration program. Note that we do EVERYTHING in functions, handling as few global
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-3)
    error = 100 #starting error
    n = 0 #counter
    while error>tol #main convergence loop
        n+=1
        v_next = Bellman(prim, res) #execute Bellman operator
        error = maximum(abs.(v_next - res.val_func)) #reset error term
        res.val_func = v_next #update value function held in results vector
    end
    println("Value functions converged in ", n, " iterations.")
end

###index of any value in k_grid
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

#Bellman operator
function Bellman(prim::Primitives, res::Results)
    @unpack β, δ, θ, nk, k_grid = prim #unpack parameters from prim struct. Improves readability.
    v_next = zeros(nk)
    k_interp = interpolate(k_grid, BSpline(Linear()))
    val_interp = interpolate(res.val_func, (BSpline(Linear())))

    for i_k = 1:nk #loop over state space
        k= k_grid[i_k] #convert state indices to state values
        budget = k^θ + (1-δ)*k #budget given current state. Doesn't this look nice?

        val(kp) = log(budget- k_interp(kp)) + β * val_interp(kp)
        obj(kp) = - val(kp)

        lower = 1.0 #lower bound on k grid to search
        upper = get_index(budget, k_grid)
        opt = optimize(obj, lower, upper) #find optimal value
        res.pol_func[i_k] = k_interp(opt.minimizer[1]) #policy function
        v_next[i_k] = -opt.minimum #update next guess of value function
    end
    v_next
end

@elapsed prim, res = Solve_model() #solve for policy and value functions #
Plots.plot(prim.k_grid, res.pol_func)
Plots.plot(prim.k_grid, res.val_func)
















#########################









#####
