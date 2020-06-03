# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 05/26/2020
# ps1_Q7.jl
# ------------------------------------------------------------------------------

#clearconsole()
using Plots
using TickTock

## Housekeeping
# Param|
β = 0.8;
δ = 1;
θ = 0.5;

tol = 1e-10;
maxit = 10000;

# Shocks
Zg = 1;
Zb = 0.2;

# Transition Probabilities
# Pyx = Pr(Z' = Zy | Z = Zx) x∈{g,b} & y∈{g,b}
Pgg = 1;
Pbg = 1 - Pgg;
Pbb = 0.926;
Pgb = 1- Pbb;
P = [Pgg Pbg
    Pgb Pbb];

# Grids
k_lb = 0.001;
k_ub = 1;
n_k = 100;
k_grid = range(k_lb,stop = k_ub,length = n_k);

n_Z = 1;
Z_grid = Zg;

## Main Code
converged = 0;
it = 1;
tick()
while converged == 0 && it < maxit
    global converged
    global it

    # Housekeeping matters for Dynamic Progamming Algorithm
    if it == 1
        # PF Space Allocation
        global pf_k = zeros(n_k,n_Z);
        global pf_c = zeros(n_k,n_Z);
        global pf_v = zeros(n_k,n_Z);
    end

    global pf_k_up = zeros(n_k,n_Z);
    global pf_c_up = zeros(n_k,n_Z);
    global pf_v_up = zeros(n_k,n_Z);

    for (i_Z, Z) in enumerate(Z_grid)
        for (i_k, k_today) in enumerate(k_grid)

            # Find optimal investment/consumption given capital level today
            global pf_v_temp = zeros(n_k);
            for (i_kpr, k_tomorrow) in enumerate(k_grid)
                y_today = Z*k_today^θ;
                c_today_temp = y_today + (1-δ)*k_today - k_tomorrow;
                v_tomorrow = pf_v[i_kpr,1];
                if c_today_temp < 0
                    pf_v_temp[i_kpr] = log(0) + β*v_tomorrow;
                else
                    pf_v_temp[i_kpr] = log(c_today_temp) + β*v_tomorrow;
                end
            end

            v_today, max_v_inx = findmax(pf_v_temp);
            kpr = k_grid[max_v_inx];
            c_today = Z*k_today^θ + (1-δ)*k_today - kpr;

            # Update PFs
            pf_k_up[i_k,i_Z] = kpr;
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

    if mod(it,20) == 0
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

pf1 = plot(k_grid,pf_v[:,1],legend = false,color=:black, title="Value Function", lw = 2);
pf2 = plot(k_grid,pf_k[:,1],legend = false,color=:black, title="Capital Investment", lw = 2);
pf3 = plot(k_grid,pf_c[:,1],legend = false,color=:black, title="Consumption", lw = 2);

pf = plot(pf1,pf2,pf3,layout=(1,3),size = (600,400), fmt = :pdf) #Size can be adjusted so don't need to mess around with 'blank space'
xlabel!("Initial Capital Stock")
dir = "/Users/philipcoyle/Documents/School/University_of_Wisconsin/SecondYear/Summer_2020/CodingBootcamp/ProblemSets/PS1/"
# savefig(dir*"testfig.pdf")
