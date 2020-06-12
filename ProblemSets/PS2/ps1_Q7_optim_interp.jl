# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 05/26/2020
# ps1_Q7_optim_interp.jl
# ------------------------------------------------------------------------------

using Optim, Parameters, Plots, Interpolations
## Housekeeping
# Directory
dir = "/Users/philipcoyle/Documents/School/University_of_Wisconsin/SecondYear/Summer_2020/CodingBootcamp/ProblemSets/PS2/"

## Structures
@with_kw struct Params
    β::Float64 = 0.99
    δ::Float64 = 0.025
    θ::Float64 = 0.36

    tol::Float64 = 1e-5
    maxit::Int64 = 10000
end

@with_kw struct Shocks
    Zg::Float64 = 1.25
    Zb::Float64 = 0.2

    # Transition Probabilities
    # Pyx = Pr(Z' = Zy | Z = Zx) x∈{g,b} & y∈{g,b}
    Pgg::Float64 = 0.977
    Pbg::Float64 = 1 - Pgg
    Pbb::Float64 = 0.926
    Pgb::Float64 = 1 - Pbb
    Tmat::Array{Float64,2} = [
        Pgg Pbg
        Pgb Pbb
    ]
end

@with_kw struct Grids
    Zg::Float64 = 1.25
    Zb::Float64 = 0.2

    k_lb::Float64 = 0.1
    k_ub::Float64 = 45.0^(1/4)
    n_k::Int64 = 21
    k_grid::Array{Float64,1} = range(k_lb, stop = k_ub, length = n_k).^4


    n_Z::Int64 = 2
    Z_grid::Array{Float64,1} = [Zg, Zb]
end

mutable struct PolFuncs
    pf_c::Array{Float64,2}
    pf_k::Array{Float64,2}
    pf_v::Array{Float64,2}
end


## Functions
function solve_model()
    P = Params()
    S = Shocks()
    G = Grids()

    @unpack n_k, n_Z = G
    # Initial Guess
    pf_c = zeros(n_k, n_Z)
    pf_k = zeros(n_k, n_Z)
    pf_v = zeros(n_k, n_Z)
    PFs = PolFuncs(pf_c, pf_k, pf_v)

    converged = 0
    it = 1
    while converged == 0 && it < P.maxit
        @unpack pf_v, pf_k, pf_c = PFs

        pf_c_up, pf_k_up, pf_v_up = Bellman(P, S, G, PFs)

        diff_v = sum(abs.(pf_v_up - pf_v))
        diff_k = sum(abs.(pf_k_up - pf_k))
        diff_c = sum(abs.(pf_c_up - pf_c))

        max_diff = diff_v + diff_k + diff_c

        if mod(it, 250) == 0 || max_diff < P.tol
            println(" ")
            println("*************************************************")
            println("AT ITERATION = ", it)
            println("MAX DIFFERENCE = ", max_diff)
            println("*************************************************")

            if max_diff < P.tol
                converged = 1
            end
        end
        # Update the policy functions
        PFs = PolFuncs(pf_c_up, pf_k_up, pf_v_up)
        # η = 0.0
        # PFs = PolFuncs(η.*pf_c .+ (1-η).*pf_c_up, η.*pf_k .+ (1-η).*pf_k_up, η.*pf_v .+ (1-η).*pf_v_up)
        it = it + 1
    end

    return G, PFs
end

function Bellman(P::Params, S::Shocks, G::Grids, PFs::PolFuncs)
    @unpack β, δ, θ = P
    @unpack Zg, Zb, Pgg, Pbg, Pbb, Pgb, Tmat = S
    @unpack n_k, k_grid, n_Z, Z_grid = G
    @unpack pf_c, pf_k, pf_v = PFs

    pf_k_up = zeros(n_k, n_Z)
    pf_c_up = zeros(n_k, n_Z)
    pf_v_up = zeros(n_k, n_Z)

    k_interp = interpolate(k_grid, BSpline(Linear()))
    v_interp = interpolate(pf_v, BSpline(Linear()))

    for (i_Z, Z) in enumerate(Z_grid)
        Pr = Tmat[i_Z, :]
        for (i_k, k_today) in enumerate(k_grid)

            y_today = Z * k_today^θ
            budget = y_today + (1 - δ) * k_today

            v_tomorrow(i_kp) = Pr[1] * v_interp(i_kp, 1) + Pr[2] * v_interp(i_kp, 2)
            val_func(i_kp) = log(budget - k_interp(i_kp)) +  β*v_tomorrow(i_kp)
            obj(i_kp) = -val_func(i_kp)
            lower = 1.0
            upper = get_index(budget, k_grid)
            opt = optimize(obj, lower, upper)

            k_tomorrow = k_interp(opt.minimizer[1])
            c_today = budget - k_tomorrow
            v_today = -opt.minimum

            # Update PFs
            pf_k_up[i_k, i_Z] = k_tomorrow
            pf_c_up[i_k, i_Z] = c_today
            pf_v_up[i_k, i_Z] = v_today

        end
    end

    return pf_c_up, pf_k_up, pf_v_up
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

function plot_pfs(dir::String, G::Grids, PFs::PolFuncs)
    @unpack k_grid = G
    @unpack pf_c, pf_k, pf_v = PFs

    pf1 = plot(k_grid,pf_v[:,1],legend = false,color=:black, label = "Good State",lw = 2);
    pf_1 = plot!(pf1,k_grid,pf_v[:,2],title="Value Function",legend = true,color=:blue, label = "Bad State",lw = 2);

    pf2 = plot(k_grid,pf_k[:,1],legend = false,color=:black, lw = 2);
    pf_2 = plot!(pf2,k_grid,pf_k[:,2],title="Capital Investment",legend = false,color=:blue, lw = 2);

    pf3 = plot(k_grid,pf_c[:,1],legend = false,color=:black, lw = 2);
    pf_3 = plot!(pf3,k_grid,pf_c[:,2],title="Consumption",legend = false,color=:blue, lw = 2);

    pf = plot(pf_1,pf_2,pf_3,layout=(1,3),size = (600,400)) #Size can be adjusted so don't need to mess around with 'blank space'
    xlabel!("Initial Capital Stock")
    # savefig(dir * "Q7_PFs.pdf")
end

## Main Code
@time G, PFs = solve_model();
plot_pfs(dir, G, PFs)
