#function that solves the model
function Solve_model()
    prim = Primitives()
    val_func = zeros(prim.na, prim.nh, prim.n_age)
    pol_func = zeros(prim.na, prim.nh, prim.n_age)
    res = Results(val_func, pol_func)
    Backward_Induct(prim, res)
    prim, res #return deliverables
end

#backward induction protocol
function Backward_Induct(prim::Primitives, res::Results)
    @unpack age_grid, n_age = prim

    for i = 1:n_age #loop over ages
        age = age_max - i + 1 #now going backward!
        j =
        println("Working on ", age)
        res.val_func[:,:,j] = Bellman(prim, res, j)
    end
end

#bellman equation
function Bellman(prim::Primitives, res::Results, j::Int64)
    @unpack n_age, na, nh, h_grid, ρ = prim
    v_next = zeros(na, nh) #value function to fill

    if j == n_age #check for maximal age
        for i_a = 1:na, i_h = 1:nh #loop over state space
            v_next[i_a, i_h] = utility(h_grid[i_h], ρ) #consumption from I = 0
            res.pol_func[i_a, i_h, j] = 0.0 #update policy fucntion
        end
    elseif j!=n_age #not in terminal age
        for i_a = 1:na, i_h = 1:nh #loop over state space
            h = h_grid[i_h]
            obj(I) = -Obj_func(prim, res, j, h, i_a, I) #objective function
            opt = optimize(obj, 0.0, 1.0) #solve
            v_next[i_a, i_h] = -opt.minimum #update value function next guess
            res.pol_func[i_a, i_h, j] = opt.minimizer #policy function
        end
    end
    v_next #return
end

#objective function: age, hc, ability, and investment given
function Obj_func(prim::Primitives, res::Results, j::Int64, h::Float64, i_a::Int64, I::Float64)
    @unpack a_grid, β, ρ, κ, μ, σ, a_grid, na, h_grid, nh = prim
    ε_grid = discretize_lognormal(μ, σ, 10) #discretized grid of possible shocks
    interp_vp = interpolate(res.val_func[i_a,:,j+1], BSpline(Linear())) #interpolated next-period value function
    a = a_grid[i_a] #ability level

    c = h * (1-I) #consumption
    val = utility(c, ρ) #flow utility

    for i = 1:10 #loop over possible HC shocks
        hprime_val = (h + a*(h*I)^κ) * ε_grid[i] #next period HC value
        hprime = get_index(hprime_val, h_grid)
        val+= (β *  interp_vp(hprime))/10 #add to continuation value, accounting for likelihood of shock
    end
    val #return total value
end


############
