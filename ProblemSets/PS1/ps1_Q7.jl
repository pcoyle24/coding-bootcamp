# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 05/26/2020
# ps1_Q7.jl
# ------------------------------------------------------------------------------

#clearconsole()
using TickTock
using Plots
## Housekeeping
# Directory
dir = "/Users/philipcoyle/Documents/School/University_of_Wisconsin/SecondYear/Summer_2020/CodingBootcamp/ProblemSets/PS1/";

# Param
β = 0.99;
δ = 0.025;|
θ = 0.36;

tol = 1e-10;
maxit = 10000;

# Shocks
Zg = 1.25;
Zb = 0.2;

# Transition Probabilities
# Pyx = Pr(Z' = Zy | Z = Zx) x∈{g,b} & y∈{g,b}
Pgg = 0.977;
Pbg = 1 - Pgg;
Pbb = 0.926;
Pgb = 1- Pbb;
P = [Pgg Pbg
    Pgb Pbb];

# Grids
k_lb = 1;
k_ub = 75;
n_k = 100;
k_grid = range(k_lb,stop = k_ub,length = n_k);

n_Z = 2;
Z_grid = [Zg, Zb];

## Main Code
converged = 0;
it = 1;
tick()
while converged == 0 && it < maxit
    global converged
    global it

    if it == 1
        # PF Space Allocation
        global pf_k = zeros(n_k,n_Z);
        global pf_c = zeros(n_k,n_Z);
        global pf_v = zeros(n_k,n_Z);
    end

    # To make updating work
    global pf_k_up = zeros(n_k,n_Z);
    global pf_c_up = zeros(n_k,n_Z);
    global pf_v_up = zeros(n_k,n_Z);

    for (i_Z, Z) in enumerate(Z_grid)
        Pr = P[i_Z,:];
        for (i_k, k_today) in enumerate(k_grid)

            # Find optimal investment/consumption given capital level today
            global pf_v_temp = zeros(n_k);
            for (i_k′, k_tomorrow) in enumerate(k_grid)
                y_today = Z*k_today^θ;
                c_today_temp = y_today + (1-δ)*k_today - k_tomorrow;
                v_tomorrow = Pr[1]*pf_v[i_k′,1] + Pr[2]*pf_v[i_k′,2]
                if c_today_temp < 0
                    pf_v_temp[i_k′] = log(0) + β*v_tomorrow;
                else
                    pf_v_temp[i_k′] = log(c_today_temp) + β*v_tomorrow;
                end
            end

            v_today, max_v_inx = findmax(pf_v_temp);
            k′ = k_grid[max_v_inx];
            c_today = Z*k_today^θ + (1-δ)*k_today - k′;

            # Update PFs
            pf_k_up[i_k,i_Z] = k′;
            pf_c_up[i_k,i_Z] = c_today;
            pf_v_up[i_k,i_Z] = v_today;
        end
    end

    diff_v = sum(abs.(pf_v_up - pf_v));
    diff_k = sum(abs.(pf_k_up - pf_k));
    diff_c = sum(abs.(pf_c_up - pf_c));

    max_diff = diff_v + diff_k + diff_c;

    if max_diff < tol
        converged = 1;
        println(" ")
        println("*************************************************")
        println("AT ITERATION = ", it)
        println("MAX DIFFERENCE = ", max_diff)
        println("*************************************************")
    end

    # Update the policy functions
    pf_v = pf_v_up;
    pf_k = pf_k_up;
    pf_c = pf_c_up;

    if mod(it,50) == 0
        println(" ")
        println("*************************************************")
        println("AT ITERATION = ", it)
        println("MAX DIFFERENCE = ", max_diff)
        println("*************************************************")
    end

    it = it + 1;
end
tock()

## Plotting

pf1 = plot(k_grid,pf_v[:,1],legend = false,color=:black, label = "Good State",lw = 2);
pf_1 = plot!(pf1,k_grid,pf_v[:,2],title="Value Function",legend = true,color=:blue, label = "Bad State",lw = 2);

pf2 = plot(k_grid,pf_k[:,1],legend = false,color=:black, lw = 2);
pf_2 = plot!(pf2,k_grid,pf_k[:,2],title="Capital Investment",legend = false,color=:blue, lw = 2);

pf3 = plot(k_grid,pf_c[:,1],legend = false,color=:black, lw = 2);
pf_3 = plot!(pf3,k_grid,pf_c[:,2],title="Consumption",legend = false,color=:blue, lw = 2);

pf = plot(pf_1,pf_2,pf_3,layout=(1,3),size = (600,400)) #Size can be adjusted so don't need to mess around with 'blank space'
xlabel!("Initial Capital Stock")
# savefig(dir*"Q7_PFs.pdf")
