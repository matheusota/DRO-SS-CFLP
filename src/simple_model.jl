function solve_simple_model(data::DataSSCFLP, app::Dict{String,Any})
    F = [f for f in 1:data.n_facilities]
    C = [c for c in 1:data.n_customers]
    
    # Sum of demands and capacities.
    total_demand = sum([c.demand for c in data.customers])
    total_capacity = sum([f.capacity for f in data.facilities])
 
    # Formulation.
    model = Model(CPLEX.Optimizer)
    @variable(model, x[f in F, c in C], Bin)
    @variable(model, y[f in F], Bin)
    @objective(model, Min, sum(data.facilities[f].fixed_cost * y[f] for f in F) + sum(data.cost[f, c] * x[f, c] for f in F, c in C))
    @constraint(model, customer_constraint[c in C], sum(x[f, c] for f in F) >= 1)
    @constraint(model, assign_constraint[c in C, f in F], x[f, c] <= y[f])
    @constraint(model, capacity_constraint[f in F], sum(data.customers[c].demand * x[f, c] for c in C) <= data.facilities[f].capacity * y[f])

    # Declare robust variables and constraints.
    delta = app["delta"]
    @variable(model, u[f in F, c in C] >= 0)
    @variable(model, theta[f in F] >= 0)
    @constraint(model, robust1[f in F], sum((data.customers[c].demand * x[f, c]) + u[f, c] for c in C) + delta * theta[f] <= data.facilities[f].capacity * y[f] )
    @constraint(model, robust2[f in F, c in C], theta[f] + u[f, c] >= data.customers[c].deviation * x[f, c])
    
    #println(model)
    set_time_limit_sec(model, 1800.0)
    stats = @timed optimize!(model)
    println(termination_status(model))
    println(objective_value(model))
    println("----------------------------------------------------")
    println("Time to solve: $(stats.time)")
    # x_is_selected = isapprox.(value.(x), 1; atol = 1e-5)
    # y_is_selected = isapprox.(value.(y), 1; atol = 1e-5)
    
    y_values = Dict()
    x_values = Dict()
    for f in F
        y_values[f] = value(y[f])
        for c in C
            x_values[(f, c)] = value(x[f, c])
            # if value(x[f, c]) > 0.0001
            #     println("x_$(f)_$(c) = $(value(x[f, c]))")
            #     println("u_$(f)_$(c) = $(value(u[f, c]))")
            # end
        end
        #println("theta_$(f) = $(value(theta[f]))")
    end
    
    return x_values, y_values
end