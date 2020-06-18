# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 06/18/2020
# lec6.jl
# ------------------------------------------------------------------------------

## Ben-Porath Problem (Human Capital Accumulaton)
using Parameters, Interpolations, Optim, Plots, Distributions, Random

# Master Code Code to run ALL functions

include("codealong_6_functions.jl")
include("codealong_6_model_live.jl")

@elapsed prim, res = Solve_model()


pol_func_avg = zeros(51)
for i = 1:prim.n_age
    pol_func_avg[i] = mean(res.pol_func[:,:,i])
end

plot(prim.age_grid, pol_func_avg)


#simulation function
function Simulate(prim::Primitives, res::Results, h_0_val::Float64, i_a::Int64, nsim::Int64)
    @unpack n_age, age_grid, a_grid, β, ρ, κ, μ, σ, a_grid, na, h_grid, nh = prim
    @unpack val_func, pol_func = res

    a = a_grid[i_a]
    dist = LogNormal(μ, σ)

    hc_grid = zeros(n_age, nsim)
    earnings_grid = zeros(n_age, nsim)
    inv_grid = zeros(n_age, nsim)

    for s = 1:nsim
        h_val = h_0_val
        for j = 1:n_age
            pol_interp = interpolate(res.pol_func[i_a, :, j], BSpline(Linear()))
            h_inx = get_index(h_val, h_grid)
            I = pol_interp(h_inx)

            hc_grid[j, s] = h_val
            earnings_grid[j, s] = h_val*(1-I)
            inv_grid[j, s] = I

            # Draw Shock
            hp = (h_val + a*(h_val*I)^κ)*rand(dist)

            # Update
            h_val = hp
        end
    end
    hc_avg = zeros(n_age)
    earnings_avg= zeros(n_age)
    inv_avg = zeros(n_age)

    for k = 1:n_age
        hc_avg[k] = mean(hc_grid[k,:])
        earnings_avg[k] = mean(earnings_grid[k,:])
        inv_avg[k] = mean(inv_grid[k,:])
    end

    return hc_avg, earnings_avg, inv_avg

end


#simulate the life course of a single person
hc_avg, earnings_avg, inv_avg = Simulate(prim, res, 110.0, 7, 10000)
plot(prim.age_grid, hc_avg)
plot(prim.age_grid, earnings_avg)
plot(prim.age_grid, inv_avg)
