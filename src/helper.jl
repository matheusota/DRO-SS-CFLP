using Random, Distributions
#import Plots

function get_epsilon_scenario(distribution, customer)
    x = (round(rand(distribution)) - customer.demand) / customer.deviation
    return min(max(-1.0, x), 1.0)
end

function set_random_data(data::DataSSCFLP, app::Dict{String,Any})
    Random.seed!(666)
    
    # # First, set the deviations
    # println("-----------------------------------------")
    # println("deviations")
    # deviation_factor = app["deviation"]
    # for c in data.customers
    #     c.deviation = round((rand(100:500) / 1000.0) * c.demand)
    #     #c.deviation = floor(deviation_factor * c.demand + 0.5)
    #     println(c.deviation)
    # end
    # println("-----------------------------------------")
    
    # Next, generate scenarios following normal distributions.
    data.n_scenarios = app["scenarios"]
    distributions = [Normal(c.demand, 0.5 * c.deviation) for c in data.customers]
    for s in 1:data.n_scenarios
        push!(data.scenarios, [get_epsilon_scenario(distributions[c], data.customers[c]) for c in 1:length(data.customers)])
    end
    #println(data.n_scenarios)
    #println(data.scenarios)
    
    return distributions
end                                

function run_simulation(data::DataSSCFLP, app::Dict{String,Any}, solution::SolutionSSCFLP, distributions)
    Random.seed!(665)
    
    # Compute fixed cost.
    fixed_cost = 0.0
    for f in solution.opened_facilities
        fixed_cost += data.facilities[f].fixed_cost
    end
    
    # Run multiple simulations
    n_failures = 0
    total_cost = 0.0
    n_trials = app["trials"]
    for s in 1:n_trials
        demands = [round(rand(distributions[c])) for c in 1:data.n_customers]
        
        has_failed = false
        curr_cost = 0.0
        acc_demand = [0.0 for f in 1:data.n_facilities]
        for c in 1:data.n_customers
            customer = data.customers[c]
            f = solution.assigned_facility[c]
            facility = data.facilities[f]
            curr_cost += (data.cost[f, c] / customer.demand) * demands[c]
            acc_demand[f] += demands[c]
            
            if acc_demand[f] >= facility.capacity + 0.001
                #println("acc demand $(acc_demand[f]) capacity $(facility.capacity)")
                has_failed = true
            end
        end
        
        if has_failed
            n_failures += 1
        end
        
        total_cost += curr_cost
    end
    
    println("number of failures $(n_failures)")
    return fixed_cost, total_cost / n_trials, n_failures / n_trials
end

function plot_solution(data::DataSSCFLP, solution::SolutionSSCFLP)
    Plots.scatter(
        [c.pos_x for c in data.customers],
        [c.pos_y for c in data.customers];
        label = nothing,
        markershape = :circle,
        markercolor = :blue,
        markersize = 2 .* (2 .+ [c.demand for c in data.customers]),
    )

    Plots.scatter!(
        [f.pos_x for f in data.facilities],
        [f.pos_y for f in data.facilities];
        label = nothing,
        markershape = :rect,
        markercolor = [(f in solution.opened_facilities) ? :red : :white for f in 1:data.n_facilities],
        markersize = [f.capacity for f in data.facilities],
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    
    for c in data.customers
        f = solution.assigned_facility[c]
        Plots.plot!(
            [c.pos_x, f.pos_x],
            [c.pos_y, f.pos_y];
            color = :black,
            label = nothing,
        )
    end
    
    annotate!.([(c.pos_x + 0.01, c.pos_y, text("$(c.demand)\$\\pm\$$(c.deviation)", :red, :left, 11)) for c in data.customers])
end