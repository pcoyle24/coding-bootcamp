using Parameters, Interpolations, Optim, Plots, Distributions, Random
include("codealong_6_functions.jl")
include("codealong_6_model.jl")

prim, res = Solve_model()

pol_func_avgs = zeros(51)
for i = 1:prim.n_age
    pol_func_avgs[i] = mean(res.pol_func[:,:,i])
end

Plots.plot(prim.age_grid, pol_func_avgs)

#simulation function
function Simulate(prim::Primitives, res::Results, h_0_val::Float64, i_a::Int64, nsim::Int64)
    @unpack n_age, age_grid, a_grid, β, ρ, κ, μ, σ, a_grid, na, h_grid, nh = prim
    @unpack val_func, pol_func = res
    a = a_grid[i_a] #level of ability
    dist = LogNormal(μ, σ) #set up log normal distribution

    #simulate some outcomes for a representative person
    hc_grid = zeros(n_age, nsim) #tracking human capital
    earnings_grid = zeros(n_age, nsim) #tracking earnings
    inv_grid = zeros(n_age, nsim) #tracking self-investment

    for s = 1:nsim #loop over simulations
        h_val = h_0_val
        for j = 1:n_age #loop over age, starting from 20 now
            interp_pol = interpolate(pol_func[i_a, :, j], BSpline(Linear())) #inteprolated policy function
            h_ind = get_index(h_val, h_grid)
            inv = interp_pol(h_ind) #investment decision

            hc_grid[j, s] = h_val #store current HC
            earnings_grid[j, s] = h_val * (1-inv) #store earnings
            inv_grid[j, s] = inv #store investment choice

            hp = (h_val + a*(h_val*inv)^κ) * rand(dist) #next-period HC
            h_val = hp #reset h in preparation for next period
        end
    end
    hc_avg, earnings_avg, inv_avg = zeros(n_age), zeros(n_age), zeros(n_age)

    for j = 1:n_age #loop over ages to get age profiles
        hc_avg[j] = mean(hc_grid[j, :])
        earnings_avg[j] = mean(earnings_grid[j, :])
        inv_avg[j] = mean(inv_grid[j, :])
    end
    hc_avg, earnings_avg, inv_avg
end


hc_avg, earnings_avg, inv_avg= Simulate(prim, res, 110.0, 7, 10000)
Plots.plot(prim.age_grid, hc_avg, ylims=(0,140), yticks =   0:20:140)
Plots.plot(prim.age_grid, earnings_avg,  ylims=(0,140), yticks =  0:20:140)
Plots.plot(prim.age_grid, inv_avg,  ylims=(0, 0.5), yticks =  0:0.05:0.5)


#####
