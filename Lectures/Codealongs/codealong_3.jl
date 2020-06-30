using Optim, Interpolations, Plots

####Brent's method
f(x, y) = (x-y)^2
opt = optimize(x->f(x, 1.0), -5.0, 5.0)
#easy!

#####Rosenbrock
function Rosenbrock(vec::Vector{Float64})
    x1, x2 = vec[1], vec[2]
    val = (1 - x1)^2 + 100 * (y - x2^2)^2
    val
end

##evaluate function at a bunch of points
x_grid = collect(-3:0.01:3)
nx = length(x_grid)
z_grid = zeros(nx, nx)

for i = 1:nx, j = 1:ny
    guess = [x_grid[i], x_grid[j]]
    z_grid[i,j] = Rosenbrock(guess)
end

#plot
Plots.surface(x_grid, x_grid, z_grid, seriescolor=:viridis, camera = (50,50))
Plots.contourf(x_grid, x_grid, log.(1 .+ z_grid), seriescolor=:inferno)

#Nelder-Mead
guess = [0.0, 0.0]
opt = optimize(Rosenbrock, guess)

###LBFGS
opt = optimize(Rosenbrock, guess, LBFGS())

#Define Gradient
function g(G, guess::Vector{Float64})
    x, y = guess[1], guess[2]
    G[1] = -2.0 * (1.0 - x) - 400.0 * (y - x^2) * x
    G[2] = 200.0 * (y - x^2)
    G #return
end

#Hessian
function h(H, guess::Vector{Float64})
    x, y = guess[1], guess[2]
    H[1,1] = 2.0 - 400.0 * y + 1200.0 * x^2
    H[1,2] = -400.0 * x
    H[2,1] = -400.0 * x
    H[2,2] = 200.0
    H #retturn
end

#Nelder_Mead
guess = [0.0, 0.0]
opt = optimize(Rosenbrock, g, h, guess_init)

####Lots of local minima
function Greiwank(x::Array{Float64,1})
    val = (1/4000)*sum(x.^2) - prod(cos.(x./sqrt(length(x)))) + 1
    val
end

##evaluate function at a bunch of points
x_grid = collect(-5:0.01:5)
nx = length(x_grid)
z_grid = zeros(nx, nx)

for i = 1:nx, j = 1:ny
    guess = [x_grid[i], x_grid[j]]
    z_grid[i,j] = Greiwank(guess)
end

##plots
Plots.surface(x_grid, x_grid, z_grid, seriescolor=:viridis, camera = (50,50))
Plots.contourf(x_grid, x_grid, z_grid, seriescolor=:inferno)

#global optimum at (0,0)
guess_init = [3.0, 3.0]
opt = optimize(Greiwank, guess_init) #this fails!

#now this works!
guess_init = [2.0, 2.0]
opt = optimize(Greiwank, guess_init) #this works!

#try multiple starts!
function Multistart()
    x_grid = collect(-5:0.5:5)
    nx = length(x_grid)
    minimum, minimizers = 100, [100, 100] #preallocate bad values for minimum and minimizes

    for i = 1:nx, j = 1:nx
        guess = [x_grid[i], x_grid[j]] #starting guess
        opt = optimize(Greiwank, guess) #nelder-mead with new starting guess
        if opt.minimum<minimum #new minimum!
            minimum = opt.minimum #update
            minimizers = opt.minimizer #update
        end
    end
    minimum, minimizers #return
end
min, minimizers = Multistart()

#####OLS example
using Distributions, Random

##run the same OLS as first class
dist = Normal(0,1)
β_0 = 1.0
β_1 = 2.0
β_2 = 3.0
n = 10000
x = rand(n).*10
x2 = x.^2
Random.seed!(1234)  ###important to remember!!
ϵ = rand(dist, n)
Y_true = β_0 .+ β_1.*x + β_2.*x2 .+ ϵ
X = hcat(ones(n), x, x2)
β_ols = inv(X' * X) * X' * Y_true

####ols-nelder function
function OLS_Nelder(β::Array{Float64,1})
    β_0, β_1, β_2 = β[1], β[2], β[3] #unpack β
    Random.seed!(1234) #DON'T USE AT FIRST
    ϵ = rand(dist, n) #draw epsilons
    Y_true = 1.0 .+ 2.0.*x + 3.0.*x2 .+ ϵ
    Y_predict = β_0 .+ β_1.*x + β_2.*x2 #true and predicted value
    error = sum((Y_true.-Y_predict).^2) #sum of squared error
    error #return
end

#do OLS with nelder-mead
guess_init = [1.0, 2.0, 3.0]
opt = optimize(OLS_Nelder, guess_init) #it works



#########################
