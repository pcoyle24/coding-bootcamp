# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 06/10/2020
# ps2.jl
# ------------------------------------------------------------------------------

using Optim, Plots

## Question 1
function Himmelblau(in::Array{Float64,1})
    x = in[1]
    y = in[2]

    out = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
    return out
end

function g(G, in::Array{Float64,1})
    x = in[1]
    y = in[2]

    G[1] = (4*x)*(x^2 + y - 11) + 2*(x + y^2 - 7)
    G[2] = 2*(x^2 + y - 11) + (4*y)*(x + y^2 - 7)
end
function h(H, in::Array{Float64,1})
    x = in[1]
    y = in[2]

    H[1,1] = 12*x^2 + 4*y - 44 + 2
    H[1,2] = 4*x + 4*y
    H[2,1] = 4*x + 4*y
    H[2,2] = 2 + 4x + 12y - 28
end

# Part A
n_grid = 1001
x_grid = range(-4, 4, length = n_grid)
y_grid = x_grid

# Preallocate space for functuon
z_grid = zeros(n_grid, n_grid)
for i in 1:n_grid
    for j = 1:n_grid
        vec_in = [x_grid[i], y_grid[j]]
        z_grid[i,j] = Himmelblau(vec_in)
    end
end

# Plot
surface(x_grid, y_grid, z_grid, camera = (50, 50))
contourf(x_grid, y_grid, log.(z_grid))

# Part B (There are 4 local minima)
for x1 = [-2.0, 2.0]
    for x2 = [-3.0, 3.0]
        guess = [x1, x2]
        opt_newton = optimize(Himmelblau, g, h,guess)
        println(opt_newton.minimizer)
    end
end


# Part C
for x1 = [-2.0, 2.0]
    for x2 = [-3.0, 3.0]
        guess = [x1, x2]
        opt_nm = optimize(Ackley, guess, NelderMead())
        println(opt_nm.minimizer)
    end
end

## Question 2
function Ackley(in::Array{Float64,1})
    x = in[1]
    y = in[2]

    out = -20*exp(-0.2*(0.5*(x^2 + y^2))^(0.5)) - exp(0.5*(cos(2*pi*x) + cos(2*pi*y))) + exp(1) + 20
    return out
end

# Part A
n_grid = 1001
x_grid = range(-4, 4, length = n_grid)
y_grid = x_grid

# Preallocate space for functuon
z_grid = zeros(n_grid, n_grid)
for i in 1:n_grid
    for j = 1:n_grid
        vec_in = [x_grid[i], y_grid[j]]
        z_grid[i,j] = Ackley(vec_in)
    end
end

# Plot
surface(x_grid, y_grid, z_grid, camera = (50, 50))
contourf(x_grid, y_grid, z_grid)

# Part B (trying different guesses)
for x1 = [-2.0, 2.0]
    for x2 = [-3.0, 3.0]
        guess = [x1, x2]
        opt_nm = optimize(Ackley, guess, NelderMead())
        opt_lbfgs = optimize(Ackley, guess, LBFGS())
        println(opt_nm.minimizer)
        println(opt_lbfgs.minimizer)
        println(" ")
    end
end

## Question 3
function Rastrigin(x::Array{Float64,1}; A = 10)
    sum = 0
    for i = 1:length(x)
        sum = sum + x[i]^2 - A*cos(2*pi*x[i])
    end

    n = length(x)
    out = A*n + sum

    return out
end

# Part A
n_grid = 1001
x_grid = range(-5.12, 5.12, length = n_grid)

# Preallocate space for functuon
z_grid_1d = zeros(n_grid,1)
for i in 1:n_grid

    z_grid_1d[i] = Rastrigin([x_grid[i]])
end

plot(x_grid,z_grid_1d)

# Part B
n_grid = 1001
x_grid = range(-5.12, 5.12, length = n_grid)

# Preallocate space for functuon
z_grid_2d = zeros(n_grid,n_grid)
for i in 1:n_grid
    for j = 1:n_grid
        x_in = [x_grid[i], x_grid[j]]
        z_grid_2d[i,j] = Rastrigin(x_in)
    end
end

surface(x_grid,x_grid,z_grid_2d, camera = (50, 50))
contourf(x_grid,x_grid,log.(z_grid_2d .+ 1))

# Part C (trying different guesses)
for x1 = [-2.0, 2.0]
    for x2 = [-3.0, 3.0]
        guess = [x1, x2]
        opt_nm = optimize(Rastrigin, guess, NelderMead())
        opt_lbfgs = optimize(Rastrigin, guess, LBFGS())
        println(opt_nm.minimizer)
        println(opt_lbfgs.minimizer)
        println(" ")

    end
end
