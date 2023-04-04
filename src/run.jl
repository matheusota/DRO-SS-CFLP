#using VrpSolver
using JuMP 
using ArgParse
using CPLEX

@static if VERSION <= v"1.6"
    using BaPCodVRPSolver
end

include("data.jl")
include("file_handler.jl")
include("robust_model.jl")
include("simple_model.jl")
include("dro_simple_model.jl")
include("dro_model.jl")
include("improved_model.jl")
include("helper.jl")

const STAT_HEADER = "statistics_cols: instance & :Optimal & cutoff & :bcRecRootDb & :bcTimeRootEval & :bcCountNodeProc & :bcRecBestDb & :bcRecBestInc & :bcTimeMain \\"

function parse_commandline(args_array::Array{String,1}, appfolder::String)
    s = ArgParseSettings(usage="##### VRPSolver #####\n\n"*
        "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
    @add_arg_table s begin
        "instance"
        help = "Instance file path"
        "--cfg", "-c"
        help = "Configuration file path"
        default = "$appfolder/../config/SSCFLP.cfg"
        "--ub","-u"
        help = "Upper bound (primal bound)"
        arg_type = Float64
        default = 10000000.0
        "--out","-o"
        help = "Path to write the solution found"
        "--batch","-b" 
        help = "batch file path"
        "--deviation","-q"
        help = "Demand deviation"
        arg_type = Float64
        default = 0.20
        "--delta","-d"
        help = "Delta value for Uncertainty Set"
        arg_type = Float64
        default = 5.0
        "--scenarios","-s"
        help = "Number of scenarios."
        arg_type = Int
        default = 24
        "--trials","-t"
        help = "Number of trials."
        arg_type = Int
        default = 10000
        "--model", "-m"
        help = "Model type"
        default = "robust"
        "--preprocess_thetas","-p" 
        help = "preprocess theta values"
        action = :store_true
    end
   return parse_args(args_array, s)
end

function run(app::Dict{String,Any})
    println("Application parameters:")
    for (arg,val) in app
        println("  $arg  =>  $(repr(val))")
    end
    flush(stdout)

    instance_name = split(basename(app["instance"]), ".")[1]

    data = read_instance_elena(app["instance"])
    #data = read_instance_generated(app["instance"])

    solution_found = false
    distributions = set_random_data(data, app)
    
    model = app["model"]
    if model == "simple"
        x,y = solve_simple_model(data, app)
    elseif model == "dro"
        x,y = solve_dro_model(data, app)
    elseif model == "simple-dro"
        x,y = solve_simple_dro_model(data, app)
    else
        x,y = solve_robust_model(data, app)
    end
    
    solution = construct_solution(data, x, y)
    check_solution(data, x)
    fixed_cost, average_cost, failure_ratio = run_simulation(data, app, solution, distributions)
    println("--------------------------------------")
    println("fixed cost $(fixed_cost) average cost $(average_cost) failure ratio $(failure_ratio)")
    println("--------------------------------------")
    
    if app["out"] != nothing
        write_simulation_results(instance_name, app["delta"], fixed_cost, average_cost, failure_ratio, app["out"])
    end
end

function construct_solution(data::DataSSCFLP, x, y)
    solution = SolutionSSCFLP([], Dict())
    for f in 1:data.n_facilities
        if y[f] >= 0.999
            push!(solution.opened_facilities, f)
        end
        
        for c in 1:data.n_customers
            if x[(f, c)] >= 0.999
                solution.assigned_facility[c] = f
            end
        end
    end
    
    return solution
end

function check_solution(data::DataSSCFLP, x_values)
    for f in 1:data.n_facilities
        acc = 0
        for c in 1:data.n_customers
            if haskey(x_values, (f, c))
                acc += x_values[(f, c)] * data.customers[c].demand
            end
        end
        
        println("Facility $f Capacity: $(data.facilities[f].capacity)")
        println("Accumulated Demand: $acc")
        
        if acc > data.facilities[f].capacity + 0.001
            println("ERROR: Invalid solution, accumulated demand is larger than facility capacity.")
        end
    end
end

function main(args)
    appfolder = dirname(@__FILE__)
    app = parse_commandline(args, appfolder)
    isnothing(app) && return
    if app["batch"] != nothing
        for line in readlines(app["batch"])
            if isempty(strip(line)) || strip(line)[1] == '#'
                continue
            end
            args_array = [String(s) for s in split(line)]
            app_line = parse_commandline(args_array, appfolder)
            run(app_line)
        end
    else
        run(app)
    end
end

if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
