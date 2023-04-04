import Unicode

mutable struct Facility
    id::String
    pos_x::Int
    pos_y::Int
    capacity::Int
    fixed_cost::Int
    variable_cost::Int
end

mutable struct Customer
    id::String
    pos_x::Int
    pos_y::Int
    demand::Int
    deviation::Float64
end

mutable struct DataSSCFLP
    n_facilities::Int
    n_customers::Int
    n_scenarios::Int
    ratio::Float64
    facilities::Array{Facility}
    customers::Array{Customer}
    cost::Matrix{Float64}
    scenarios::Array{Array{Float64}}
end

mutable struct SolutionSSCFLP
    opened_facilities::Array{Int}
    assigned_facility::Dict{Int,Int}
end

# Euclidian distance
function distance(facility::Facility, customer::Customer)
    x_sq = (facility.pos_x - customer.pos_x)^2
    y_sq = (facility.pos_y - customer.pos_y)^2
    return floor(sqrt(x_sq + y_sq) + 0.5)
end

function distance(c1::Customer, c2::Customer)
    x_sq = (c1.pos_x - c2.pos_x)^2
    y_sq = (c1.pos_y - c2.pos_y)^2
    return floor(sqrt(x_sq + y_sq) + 0.5)
end