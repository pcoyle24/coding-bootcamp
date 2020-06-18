#function that solves the model
function Solve_model()
    prim = Primitives()
    val_func = zeros(prim.na, prim.nh, prim.n_age)
    pol_func = zeros(prim.na, prim.nh, prim.n_age)
    res = Results(val_func, pol_func)

    # Do the Dynmaic Programming Part
    Backward_Induct(prim, res)


    return prim, res
end

#backward induction protocol
function Backward_Induct(prim::Primitives, res::Results)
    @unpack age_grid, n_age = prim

    for i = 1:n_age
        j = n_age - i + 1 # Flip the order of the loop (Possibly do n_age:-1:1)
        age = age_grid[n_age] - i + 1
        println("Working on age ", age)

        res.val_func[:,:,j] = Bellman(prim, res, j)
    end
end

#bellman equation
function Bellman(prim::Primitives, res::Results, j::Int64)
    @unpack n_age, na, nh, h_grid, ρ = prim

    v_next = zeros(na, nh)

    # Loop over state space
    if j == n_age
        for i_a = 1:na
            for i_h = 1:nh
                v_next[i_a, i_h] = utility(h_grid[i_h], ρ)
                res.pol_func[i_a, i_h, j] = 0.0
            end
        end
    elseif j != n_age
        for i_a = 1:na
            for i_h = 1:nh
                h = h_grid[i_h]

                # Define an anynonomus Objective Function to solve
                obj(I) = -Obj_func(prim, res, j, h, i_a, I)

                #Maximize
                opt = optimize(obj, 0.0, 1.0)
                v_next[i_a, i_h] = -opt.minimum
                res.pol_func[i_a, i_h, j] = opt.minimizer
            end
        end
    end

    return v_next

end

#objective function: age, hc, ability, and investment given
function Obj_func(prim::Primitives, res::Results, j::Int64, h::Float64, i_a::Int64, I::Float64)
    @unpack a_grid, β, ρ, κ, μ, σ, a_grid, na, h_grid, nh = prim

    ε_grid = discretize_lognormal(μ, σ, 10)
    interp_vp = interpolate(res.val_func[i_a, :, j+1], BSpline(Linear()))
    a = a_grid[i_a]

    c = h*(1-I)
    val = utility(c, ρ)

    # Loop over human capital shocks to  get expectation
    for i = 1:10
        hprime_val = (h + a*(h*I)^κ)*ε_grid[i]
        hprime_ind = get_index(hprime_val, h_grid)
        val+= interp_vp(hprime_ind)/10
    end

    return val

end
