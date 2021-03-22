using JuMP, Gurobi

function filtro(input::Vector{<:Number},weights)
    return sum(input.*weights[end:-1:1])
end

function compute_filter_weights(input::Vector{<:Number}, output::Vector{<:Number}, zerovalue_fixed::Boolean, filter_length::Integer)

    model = Model(optimizer_with_attributes(
        Gurobi.Optimizer, "OutputFlag" => 0)
    )

    @variable(model, 0<=w[1:filter_length+1]<=1)
    for i=1:filter_length
        @constraint(model,w[i]>=w[i+1])
    end

    if zerovalue_fixed
        @constraint(model, w[1]==1)
    end

    @constraint(model, w[filter_length+1]==0)

    @objective(model,Min,sum((output[filter_length+1:end] - [filtro(input[i:i+filter_length],w) for i=1:length(input)-filter_length]).^2))

    optimize!(model)
    println(termination_status(model))

    if termination_status(model) == MOI.OPTIMAL
        return value.(w)
    else
        error("The model was not solved correctly.")
    end
end

function compute_active_weights(incidence::Vector{<:Number}, active::Vector{<:Number}; filter_length::Integer=30)
    return compute_filter_weights(incidence, active, true, filter_length)
end

function compute_icu_weights(incidence::Vector{<:Number}, icu_occupancy::Vector{<:Number}; filter_length::Integer=30)
    return compute_filter_weights(incidence, icu_occupancy, false, filter_length)
end
