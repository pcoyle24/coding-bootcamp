

########Question 6 bootstrap
using Distributed

workers()
addprocs(3) #add new workers!

##make sure that all our cores have access to all the functions we need
@everywhere using Optim, NLSolversBase, Random, SharedArrays, CSV, DataFrames, Statistics
@everywhere using LinearAlgebra: diag
@everywhere dir = "C:\\Users\\Garrett\\Box Sync\\Julia_class\\Problem Sets"
@everywhere cd(dir)

#bootstrap protocol
@everywhere function Simulate(seed::Int64)
    println(seed)
    temp = CSV.read("$dir\\lwage.csv"; datarow=1) #read in data
    data = convert(Array{Float64,2}, temp)

    #random sampling procerdure
    Random.seed!(seed)
    data_boot = Any[] #preallocate random sub-sample
    n = length(data[:,1])
    for i = 1:n
        check = rand()
        if check>0.5 #coin flip
            push!(data_boot, data[i,:]) #if pass coin-flip, add line to subsample
        end
    end
    data_boot = hcat(data_boot...)'  #bootstrapped sample w/o replacement

    #everything else from here is the exact same
    Y = data_boot[:,1]
    n = length(Y)
    X = hcat(ones(n), data_boot[:,2], data_boot[:,3], data_boot[:,3].^2)
    nvar = 4

    function Log_Likelihood(X, Y, β, log_σ)
        σ = exp(log_σ)
        llike = -n/2*log(2π) - n/2* log(σ^2) - (sum((Y - X * β).^2) / (2σ^2))
        llike = -llike
    end

    guess_init = [2.0, 0.5, 0.03, 0.0, -0.693]
    func = TwiceDifferentiable(vars -> Log_Likelihood(X, Y, vars[1:nvar], vars[nvar + 1]), guess_init; autodiff=:forward)
    opt = optimize(func, guess_init; g_tol = 1e-5)
    parameters = Optim.minimizer(opt)
    parameters[nvar+1] = exp(parameters[nvar+1])
    parameters #return
end

#array of parameter guesses that all cores have access to
boot = SharedArray{Float64}(100, 5)

#main distributed loop
@sync @distributed for seed = 1:100
    boot[seed,:] = Simulate(seed)
end

#compute bootstrapped standard errors
se = zeros(5)
for i = 1:5
    se[i] = std(boot[:,i])/10
end


##########
