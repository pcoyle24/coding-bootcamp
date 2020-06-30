#######
using Random, Distributions, Plots

dir = "C:/Users/Garrett/Documents/Grad_School/Julia_class/Problem Sets/Solutions"

#####problem 1
using Optim, NLSolversBase, Random, CSV, DataFrames
using LinearAlgebra: diag
Random.seed!(0);                            # Fix random seed generator for reproducibility
dir = "C:\\Users\\Garrett\\Box Sync\\Julia_class\\Problem Sets" #directory to read in CSV file from
temp = CSV.read("$dir\\lwage.csv"; datarow=1) #read in
data = convert(Array{Float64,2}, temp) #convert to array object for now
Y = data[:,1] #Y variables
n = length(Y) #number of lines
X = hcat(ones(n), data[:,2], data[:,3], data[:,3].^2) #X variables
nvar = 4

#LL function
function Log_Likelihood(X, Y, β, log_σ)
    σ = exp(log_σ)
    llike = -n/2*log(2π) - n/2* log(σ^2) - (sum((Y - X * β).^2) / (2σ^2)) #from online resource
    llike = -llike #flip to negative
end

guess_init = [2.0, 0.5, 0.03, 0.0, -0.693] #startring guess from OLS
func = TwiceDifferentiable(vars -> Log_Likelihood(X, Y, vars[1:nvar], vars[nvar + 1]), guess_init; autodiff=:forward); #
opt = optimize(func, guess_init; g_tol = 1e-5) #optimize
parameters = Optim.minimizer(opt)
parameters[nvar+1] = exp(parameters[nvar+1]) #estimated parameters
numerical_hessian = hessian!(func,parameters) #compute hessian
var_cov_matrix = inv(numerical_hessian) #compute var-covar matrix
stand_errs = diag(var_cov_matrix) #obtain standard errors

###problem 2 -- matching
function Match_slip(n::Int64, sims::Int64)
    results = zeros(sims)
    for i = 1:sims
        people, slips = collect(1:1:n), collect(1:1:n) #initialize people and slips
        slips = shuffle(slips) #shuffle slips
        count = 0 #number of matches
        for j = 1:n
            count += (people[j]==slips[j]) #add if a match
        end
        results[i] = count #update MC results vector
    end
    results
end

#obtain results and plot
res_10, res_20 = Match_slip(10, 10000), Match_slip(20, 10000)
Plots.histogram(res_10)
Plots.savefig("$dir/ps4_1a.png")
Plots.histogram(res_20)
Plots.savefig("$dir/ps4_1b.png")


###problem 3 -- retirement savings
function Retirement(save::Float64, sims::Int64, rng::Int64)
    dist = Normal(0.06, 0.06) #distribution of investment return
    results = zeros(sims)

    for i = 1:sims
        earnings, savings = 100, 100 #earnings/savings at 30

        for j = 1:37 #gets us through age 67
            raise, inv_return = 0.03, 0.06 #no RNG
            if rng == 1
                raise, inv_return = rand()*6/100, rand(dist)
            end
            savings *= (1 + inv_return) #adjust savings
            earnings *= (1 + raise) #adjust earnings
            savings += (earnings * save) #add to savings
        end
        results[i] = savings/earnings #divide by salary
    end
    results
end

#11.25, no RNG
res_norng = Retirement(.1125, 10000, 0)
mean(res_norng)

#with rng
res_11 = Retirement(.1125, 10000, 1) #11.25
Plots.histogram(res_11)
Plots.savefig("$dir/ps4_2a.png")

res_15 = Retirement(.15, 10000, 1) #15
Plots.histogram(res_15)
Plots.savefig("$dir/ps4_2b.png")










########
