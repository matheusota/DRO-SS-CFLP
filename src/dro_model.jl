import Random

function solve_dro_model(data::DataSSCFLP, app::Dict{String,Any})
    Random.seed!(0)
    
    F = [f for f in 1:data.n_facilities]
    #F = [1]
    C = [c for c in 1:data.n_customers]
    S = [s for s in 1:data.n_scenarios]
    
    # Sum of demands and capacities.
    total_demand = sum([c.demand for c in data.customers])
    total_capacity = sum([f.capacity for f in data.facilities])
    
    # Set additional robust data.
    dataRO = DataRO([])
    set_robust_data(data, app, dataRO)
 
    # Formulation.
    model = VrpModel()
    @variable(model.formulation, x[f in F, c in C], Int)
    @variable(model.formulation, y[f in F], Int)
    @variable(model.formulation, w[s in S])
    @variable(model.formulation, alpha[s in S, c in C] >= 0)
    @variable(model.formulation, beta[s in S, c in C] >= 0)
    @variable(model.formulation, eta >= 0)
    @objective(model.formulation, Min, 100.0 * eta 
    + sum((1.0 / data.n_scenarios) * w[s] for s in S)
    + sum(data.facilities[f].fixed_cost * y[f] for f in F) 
    + sum(data.cost[f, c] * x[f, c] for f in F, c in C))
    @constraint(model.formulation, customer_constraint[c in C], sum(x[f, c] for f in F) >= 1)
 
    # Add DRO constraints.
    @constraint(model.formulation, DRO1[s in S], 
    w[s] >= sum((data.cost[f, c] / data.customers[c].demand) * data.customers[c].deviation * data.scenarios[s][c] * x[f, c] for f in F, c in C)
    + sum((1.0 - data.scenarios[s][c]) * alpha[s, c] + (1.0 + data.scenarios[s][c]) * beta[s, c] for c in C))
    
    @constraint(model.formulation, DRO2[s in S, c in C],
    eta >= sum((data.cost[f, c] / data.customers[c].demand) * data.customers[c].deviation * x[f, c] for f in F, c in C) - alpha[s, c] -beta[s, c])
    @constraint(model.formulation, DRO3[s in S, c in C],
    eta >= -sum((data.cost[f, c] / data.customers[c].demand) * data.customers[c].deviation * x[f, c] for f in F, c in C) + alpha[s, c] + beta[s, c])
         
    # Build the VRP-Solver RCSP graphs.
    
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
    
    # enum_paths, complete_form = get_complete_formulation(model, app["cfg"])
    # open("SSCFLP/data/myfile.txt", "w") do io
    #     redirect_stdout(io) do
    #         print_enum_paths(enum_paths)
    #     end
    # end
    # println(complete_form)
    # complete_form.solver = CplexSolver() # set MIP solver
    # solve(complete_form)
    # println("Objective value: $(getobjectivevalue(complete_form))\n")
    
    # Solve the problem.
    optimizer = VrpOptimizer(model, app["cfg"], split(basename(app["instance"]), ".")[1])
    set_cutoff!(optimizer, app["ub"])
    (status, solution_found) = optimize!(optimizer)
    
    # Construct solution and return it.
    return construct_solution(F, C, optimizer, x, y)
end