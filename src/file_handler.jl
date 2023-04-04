import Unicode
using Random

function read_instance_generated(filePath::String, fileString::String="")
    if fileString == ""
        str = Unicode.normalize(read(filePath, String))
    else  # fileString contain the actual file data
        str = Unicode.normalize(fileString)
    end

    aux = split(str, '\n')
    
    data = DataSSCFLP(0, 0, 0, 0.0, [], [], zeros((1, 1)), [])

    dim = 0
    for i in 1:length(aux)
        if occursin("#customers", aux[i])
            tmp = [parse(Float64, m.match) for m in eachmatch(r"[+-]?([0-9]*[.])?[0-9]+", aux[i])]
            data.n_customers = Int(tmp[1])
            data.n_facilities = Int(tmp[2])
            data.ratio = tmp[3]
            data.cost = zeros((data.n_facilities, data.n_customers))
        # Facilities
        elseif occursin("capacity fixcost varcost xcoord ycoord name", aux[i])
            j = i + 1
            while occursin(r"\d", aux[j])
                tmp = [parse(Float64, m.match) for m in eachmatch(r"[+-]?([0-9]*[.])?[0-9]+", aux[j])]
                name = last(split(aux[j], " "))
                f = Facility("", 0, 0, 0, 0, 0)
                f.id = name
                f.capacity = tmp[1]
                f.fixed_cost = tmp[2]
                f.variable_cost = tmp[3]
                f.pos_x = tmp[4]
                f.pos_y = tmp[5]
                push!(data.facilities, f)
                j += 1
            end
        # Customers
        elseif occursin("demand xcoord ycoord name", aux[i])
            j = i + 1
            while occursin(r"\d", aux[j]) 
                tmp = [parse(Float64, m.match) for m in eachmatch(r"[+-]?([0-9]*[.])?[0-9]+", aux[j])]
                name = last(split(aux[j], " "))
                c = Customer("", 0, 0, 0, 0)
                c.id = name
                c.demand = tmp[1]
                c.pos_x = tmp[2]
                c.pos_y = tmp[3]
                push!(data.customers, c)
                j += 1
            end
        # Cost
        elseif occursin("[MATRIX]", aux[i])
            j = i + 2
            f = 1
            while occursin(r"\d", aux[j]) 
                tmp = [parse(Float64, m.match) for m in eachmatch(r"[+-]?([0-9]*[.])?[0-9]+", aux[j])]
                for c in 1:length(tmp)
                    data.cost[f, c] = floor(tmp[c] + 0.5)
                end
                j += 1
                f += 1
            end
        else
            i += 1
        end
    end

    return data
end

function read_instance_elena(filePath::String, fileString::String="")
    if fileString == ""
        str = Unicode.normalize(read(filePath, String))
    else  # fileString contain the actual file data
        str = Unicode.normalize(fileString)
    end
    breaks_in = [' '; ':'; '\r'; '\n']
    aux = split(str, breaks_in; limit=0, keepempty=false)
    numbers = [parse(Float64, x) for x in aux]
    
    data = DataSSCFLP(0, 0, 0, 0.0, [], [], zeros((1, 1)), [])
    
    # Number of customers and facilities.
    data.n_customers = Int(numbers[1])
    data.n_facilities = Int(numbers[2])
    data.cost = zeros((data.n_facilities, data.n_customers))
    data.customers = [Customer("", 0, 0, 0, 0) for c in 1:data.n_customers]
    data.facilities = [Facility("", 0, 0, 0, 0, 0) for f in 1:data.n_facilities]
    i = 3
    
    # Cost matrix.
    for c in 1:data.n_customers
        for f in 1:data.n_facilities
            data.cost[f, c] = 10 * Int(numbers[i])
            i += 1
        end
    end
    
    # Customer demands.
    for c in 1:data.n_customers
        data.customers[c].demand = Int(numbers[i])
        i += 1
    end
    
    # Fixed opening costs.
    for f in 1:data.n_facilities
        data.facilities[f].fixed_cost = Int(numbers[i])
        i += 1
    end
    
    # Facilities capacities.
    for f in 1:data.n_facilities
        data.facilities[f].capacity = 2 * Int(numbers[i])
        i += 1
    end
    
    # Customer deviations.
    for c in 1:data.n_customers
        data.customers[c].deviation = numbers[i]
        i += 1
    end

    # Random.seed!(666)
    
    # println(filePath)
    # open(filePath, "a") do file
    #     println("-----------------------------------------")
    #     println("deviations")
    #     write(file, "\n")
    #     for c in data.customers
    #         c.deviation = round((rand(100:500) / 1000.0) * c.demand)
    #         println(c.deviation)
    #         write(file, "$(c.deviation)\n")
    #     end
    #     println("-----------------------------------------")
    # end

    return data
end

function write_simulation_results(instance_name, delta, fixed_cost, average_cost, failure_ratio, output_file)
    open(output_file, "a") do file
        write(file, "$(instance_name) & $(delta) & $(fixed_cost) & $(average_cost) & $(failure_ratio)\n")
    end
end