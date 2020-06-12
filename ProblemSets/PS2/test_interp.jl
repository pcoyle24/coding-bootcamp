using Optim, Plots, Interpolations, Statistics

function get_inx(val::Float64, grid::Array{Float64,1})
    if val >= grid[length(grid)]
        inx = length(grid)
    elseif val <= grid[1]
        inx = 1
    else
        inx_upper = findfirst(x->x>val, grid)
        inx_lower = inx_upper - 1
        val_upper, val_lower = grid[inx_upper], grid[inx_lower]
        inx = inx_lower + (val - val_lower) / (val_upper - val_lower)
    end
end

function approx_error(f::Function, grid::Array{Float64,1}, grid_fine::Array{Float64,1})
    f_grid = f.(grid)
    f_interp = interpolate(f_grid,BSpline(Linear()))
    f_apx = zeros(length(grid_fine),1)

    for (i, grid_i) = enumerate(grid_fine)
        grid_inx = get_inx(grid_i, grid)
        f_apx[i] = f_interp(grid_inx)
    end

    f_out = abs.(f_apx .- f.(grid_fine))
    f_out_avg = mean(f_out)

    return f_out, f_out_avg
end



function opt_grid(grid::Array{Float64,1})
    grid[1] = 0
    grid[end] =100

    f(x) = log(x+1)
    grid_fine = collect(0:0.1:100)

    ~, out = approx_error(f, grid, grid_fine)

    return out
end


x_grid = collect(range(0,100, length = 11))
x_grid_fine = collect(0:0.1:100)
f(x) = log(x+1)

# Part A
error, avg_error = approx_error(f, x_grid, x_grid_fine)
plot(x_grid, f.(x_grid), label = "Interpolated Function", color=:blue,lw = 2)
plot!(x_grid_fine, f.(x_grid_fine), label = "True Function",  color=:black,lw = 2)
plot!(x_grid_fine, error, label = "Approximation Error",  color=:red,lw = 2)


# Part C
opt = optimize(opt_grid, x_grid)
error, avg_error = approx_error(f, opt.minimizer, x_grid_fine)
plot(opt.minimizer, f.(opt.minimizer), label = "Interpolated Function", color=:blue,lw = 2)
plot!(x_grid_fine, f.(x_grid_fine), label = "True Function",  color=:black,lw = 2)
plot!(x_grid_fine, error, label = "Approximation Error",  color=:red,lw = 2)
