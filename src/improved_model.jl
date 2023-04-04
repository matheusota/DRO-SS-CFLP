function build_improved_model(data::DataSSCFLP, app::Dict{String,Any})
    F = [f for f in 1:data.n_facilities]
    C = [c for c in 1:data.n_customers]
    
    # Sum of demands and capacities.
    total_demand = sum([c.demand for c in data.customers])
    total_capacity = sum([f.capacity for f in data.facilities])
 
    # Formulation.
    model = VrpModel()
    @variable(model.formulation, x[f in F, c in C], Int)
    @variable(model.formulation, y[f in F], Int)
    @objective(model.formulation, Min, sum(data.facilities[f].fixed_cost * y[f] for f in F) + sum(data.cost[f, c] * x[f, c] for f in F, c in C))
    @constraint(model.formulation, customer_constraint[c in C], sum(x[f, c] for f in F) >= 1)
    @constraint(model.formulation, assign_constraint[c in C, f in F], x[f, c] <= y[f])
    #@constraint(model.formulation, capacity_constraint[f in F], sum(data.customers[c].demand * x[f, c] for c in C) <= data.facilities[f].capacity * y[f])
 
    # Build the VRP-Solver RCSP graph.
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
            if f > 1
                set_resource_bounds!(G, f, cap_res_id, 0.0, Float64(total_capacity - total_demand))
            end
            
            arc1 = add_arc!(G, f, f + 1)
            arc2 = add_arc!(G, f, f + 1)
            set_arc_consumption!(G, arc1, cap_res_id, 0)
            set_arc_consumption!(G, arc2, cap_res_id, data.facilities[f].capacity)
            add_arc_var_mapping!(G, arc1, y[f])
        end
        
        set_resource_bounds!(G, sink, cap_res_id, 0.0, Float64(total_capacity - total_demand))

        return G
    end
    
    function build_graph()
        n_nodes = data.n_facilities * (data.n_customers + 2)
        source = 0
        sink = n_nodes - 1
        L = 1
        U = 1
        println(data.n_facilities)
        println(data.n_customers)
        println(n_nodes)
        println(source)
        println(sink)
        G = VrpGraph(model, [i for i in source:sink], source, sink, (L, U))
        customers_resource = add_resource!(G, main = true, disposable = true)
        facilities_resource = add_resource!(G, main = true, disposable = true)
        packing_sets = [Array{Int64,1}() for c in C]
        acc_capacity = 0

        # Build the graph for each facility.
        for f in F
            first_node = (f - 1) * (data.n_customers + 2)
            last_node = f * (data.n_customers + 2) - 1
            
            if f == length(F)
                set_resource_bounds!(G, last_node, customers_resource, total_capacity, total_capacity)
            else
                set_resource_bounds!(G, last_node, customers_resource, acc_capacity, acc_capacity + data.facilities[f].capacity)
            end
            set_resource_bounds!(G, last_node, facilities_resource, 0, total_capacity - total_demand)
            
            if f > 1
                add_arc!(G, first_node - 1, first_node)
            end
            
            for c1 in 0:data.n_customers
                u = first_node + c1
                set_resource_bounds!(G, u, customers_resource, acc_capacity, acc_capacity + data.facilities[f].capacity)
                set_resource_bounds!(G, u, facilities_resource, 0, total_capacity - total_demand)
                
                for c2 in (c1 + 1):(data.n_customers + 1)
                    v = first_node + c2
                    arc = add_arc!(G, u, v)
                    
                    # x arcs.
                    if c2 < data.n_customers + 1
                        set_arc_consumption!(G, arc, customers_resource, data.customers[c2].demand)
                        add_arc_var_mapping!(G, arc, x[f, c2])
                    end
                    
                    # y arcs.
                    if c1 == 0
                        if c2 < data.n_customers + 1
                            add_arc_var_mapping!(G, arc, y[f])
                        else
                            set_arc_consumption!(G, arc, facilities_resource, data.facilities[f].capacity)
                        end
                    end
                end
                
                if c1 > 0
                    push!(packing_sets[c1], u)
                end
            end
            
            acc_capacity += data.facilities[f].capacity
        end

        return G, packing_sets
    end
    
    #facility_graph = add_graph!(model, build_capacity_graph())
    graph, packing_sets = build_graph()
    add_graph!(model, graph)
    #println(graph)
    
    #println(packing_sets)
    set_vertex_packing_sets!(model, [[(graph, Int64(v)) for v in packing_sets[c]] for c in C])
    enable_rank1_cuts!(model)
    # @expression(model.formulation, branch_expr, sum(y[f] for f in F))
    # set_branching_priority!(model, branch_expr, "branch_expr", 3)
    # enable_packset_ryanfoster_branching!(model, 2)
    enable_resource_consumption_branching!(model, 3) 
    set_branching_priority!(model, "y", 2) 
    set_branching_priority!(model, "x", 1)
 
    return (model, x, y)
end