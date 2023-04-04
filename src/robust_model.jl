import Random

mutable struct FacilityRO
    capacities::Array{Int}
    thetas::Array{Int}
    demands::Array{Array{Int}}
end

mutable struct DataRO
    facilities::Array{FacilityRO}
end


function set_robust_data(data::DataSSCFLP, app::Dict{String,Any}, dataRO::DataRO)
    # Computes demand deviations.
    preprocess_thetas = app["preprocess_thetas"]
    delta = app["delta"]
    
    # Compute initial theta values.
    thetas = Set([0.0])
    for c in data.customers
        push!(thetas, c.deviation)
    end
    
    println(thetas)
    println(preprocess_thetas)
    # Populate demands and capacities. We filter the thetas according to the procedure of Pessoa (2020).
    for f in data.facilities
        facility_ro = FacilityRO([], [], [])
        
        for theta in thetas
            if preprocess_thetas && theta > 0.0
                allowed_customers = [c for c in data.customers if c.deviation >= theta]
                
                println("Facility $f theta $theta")
                println(allowed_customers)
                
                # Not enough customers.
                if length(allowed_customers) < delta
                    continue
                end
                
                filtered_demands = [c.demand + max(0.0, c.deviation - theta) for c in data.customers]
                sort!(filtered_demands)
                
                println(filtered_demands)
                
                # Set Y_theta is empty.
                println(sum([d for d in filtered_demands[1:Int(ceil(delta))]]))
                println(f.capacity - delta * theta)
                if sum([d for d in filtered_demands[1:Int(ceil(delta))]]) > f.capacity - delta * theta
                    continue
                end
            end
            
            # Set Y_theta is not empty, set capacity and demands.
            new_demands = [c.demand + max(0.0, c.deviation - theta) for c in data.customers]
            push!(facility_ro.thetas, theta)
            push!(facility_ro.capacities, f.capacity - theta * delta)
            push!(facility_ro.demands, new_demands)
        end
        
        push!(dataRO.facilities, facility_ro)
    end
end

function construct_solution(F, C, optimizer, x, y)
    # Construct solution and return it.
    y_values = Dict()
    x_values = Dict()
    for f in F
        y_values[f] = get_value(optimizer, y[f])
        for c in C
            x_values[(f, c)] = get_value(optimizer, x[f, c])
        end
    end
    
    return x_values, y_values
end
    
function solve_robust_model(data::DataSSCFLP, app::Dict{String,Any})
    Random.seed!(0)
    
    F = [f for f in 1:data.n_facilities]
    #F = [1]
    C = [c for c in 1:data.n_customers]
    
    # Sum of demands and capacities.
    total_demand = sum([c.demand for c in data.customers])
    total_capacity = sum([f.capacity for f in data.facilities])
    
    # Set additional robust data.
    dataRO = DataRO([])
    set_robust_data(data, app, dataRO)
 
    # Formulation.
    model = VrpModel()
    @variable(model.formulation, x[f in F, c in C], Int)
    #@variable(model.formulation, w, Int)
    @variable(model.formulation, y[f in F], Int)
    @objective(model.formulation, Min, sum(data.facilities[f].fixed_cost * y[f] for f in F) + sum(data.cost[f, c] * x[f, c] for f in F, c in C))
    @constraint(model.formulation, customer_constraint[c in C], sum(x[f, c] for f in F) >= 1)
    #@constraint(model.formulation, tmp, w >= sum(2 * y[f] for f in F))
    #@constraint(model.formulation, assign_constraint[c in C, f in F], x[f, c] <= y[f])
    #@constraint(model.formulation, w_constraint[f in F], w[f] == 1 - y[f])
    #@constraint(model.formulation, y_bounds[f in F], y[f] >= 0)
    #@constraint(model.formulation, w_constraint[c in C, f in F], w[f, c] == x[f, c])
    #@constraint(model.formulation, capacity_constraint[f in F], sum(data.customers[c].demand * x[f, c] for c in C) <= data.facilities[f].capacity * y[f])
 
    # Build the VRP-Solver RCSP graphs.
    
    # Graph for ensuring facility capacities
    function build_capacity_graph()
        n_nodes = data.n_facilities + 1
        source = 1
        sink = data.n_facilities + 1
        L = 1
        U = 1

        G = VrpGraph(model, [i for i in source:sink], source, sink, (L, U))
        cap_res_id = add_resource!(G, main = true)

        # Add facility arcs.
        for f in F
            arc1 = add_arc!(G, f, f + 1)
            arc2 = add_arc!(G, f, f + 1)
            set_arc_consumption!(G, arc1, cap_res_id, data.facilities[f].capacity)
            set_arc_consumption!(G, arc2, cap_res_id, 0.0)
            set_arc_resource_bounds!(G, arc1, cap_res_id, 0.0, Float64(total_capacity - total_demand))
            set_arc_resource_bounds!(G, arc2, cap_res_id, 0.0, Float64(total_capacity - total_demand))
            add_arc_var_mapping!(G, arc1, w[f])
        end

        return G
    end
    
    # Graph for ensuring the assigned customers respect the capacity of each facility.
    function build_graph(f::Int)
        facility_ro = dataRO.facilities[f]
        n_modes = length(facility_ro.thetas)
        n_nodes = n_modes * (data.n_customers + 1) + 1
        source = 0
        sink = n_nodes
        L = 0
        U = 1

        arc_plus_ids = [[] for c in C]
        G = VrpGraph(model, [i for i in source:sink], source, sink, (L, U))
        cap_res_id = add_resource!(G, main = true)

        # Add arcs between customers.
        for m in 1:n_modes
            for c in 1:data.n_customers
                u = (m - 1) * (data.n_customers + 1) + c
                v = u + 1
                
                arc1 = add_arc!(G, u, v)
                arc2 = add_arc!(G, u, v)
                set_arc_consumption!(G, arc1, cap_res_id, facility_ro.demands[m][c])
                set_arc_consumption!(G, arc2, cap_res_id, 0)
                set_arc_resource_bounds!(G, arc1, cap_res_id, 0.0, Float64(facility_ro.capacities[m]))
                set_arc_resource_bounds!(G, arc2, cap_res_id, 0.0, Float64(facility_ro.capacities[m]))
                add_arc_var_mapping!(G, arc1, x[f, c])
                push!(arc_plus_ids[c], arc1)
            end
        end
        
        # Add arcs incident to source and sink
        for m in 1:n_modes
            first = (m - 1) * (data.n_customers + 1) + 1
            last = first + data.n_customers
            
            first_arc = add_arc!(G, source, first)
            last_arc = add_arc!(G, last, sink)
            set_arc_consumption!(G, first_arc, cap_res_id, 0)
            set_arc_consumption!(G, last_arc, cap_res_id, 0)
            set_arc_resource_bounds!(G, last_arc, cap_res_id, 0.0, Float64(facility_ro.capacities[m]))
            add_arc_var_mapping!(G, first_arc, y[f])
        end

        return G, arc_plus_ids
    end
    
    #facility_graph = build_capacity_graph()
    #add_graph!(model, facility_graph)
    #println(facility_graph)
    
    graphs = []
    packing_sets = [[] for c in C]
    for f in F
        if length(dataRO.facilities[f].thetas) == 0
            println("Error: Facility $f has 0 theta values.")
            return
        end
        
        G, arc_plus_ids = build_graph(f)
        add_graph!(model, G)
        #println(G)
        push!(graphs, G)
        
        for c in C
            for a in arc_plus_ids[c]
                push!(packing_sets[c], (f, a))
            end
        end
    end
    
    #println(packing_sets)
    
    set_arc_packing_sets!(model, [[(graphs[f], a) for (f, a) in packing_sets[c]] for c in C])
    enable_rank1_cuts!(model)
    #@expression(model.formulation, branch_expr, sum(y[f] for f in F))
    #set_branching_priority!(model, branch_expr, "branch_expr", 3)
    #enable_packset_ryanfoster_branching!(model, 4)
    enable_resource_consumption_branching!(model, 1000)
    for f in F
        @expression(model.formulation, branch_expr, sum(y[i] for i in 1:f))
        set_branching_priority!(model, branch_expr, "y_branch_expr_$f", 5 + f)  
    end
    set_branching_priority!(model, "y", 2) 
    set_branching_priority!(model, "x", 1)
    
    # Solve the problem.
    optimizer = VrpOptimizer(model, app["cfg"], split(basename(app["instance"]), ".")[1])
    set_cutoff!(optimizer, app["ub"])
    stats = @timed optimize!(optimizer)
    status, solution_found = stats.value
    println("----------------------------------------------------")
    println("Time to solve: $(stats.time)")
    
    # Construct solution and return it.
    return construct_solution(F, C, optimizer, x, y)
end