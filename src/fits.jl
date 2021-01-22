using JuMP, Ipopt

function filtro(input::Vector{<:Number},weights)
    return sum(input.*weights[end:-1:1])
end

function compute_filter_weights(input::Vector{<:Number}, output::Vector{<:Number}, filter_length::Integer)

    model = Model(optimizer_with_attributes(
        Ipopt.Optimizer, "print_level" => 0)
    )

    @variable(model, 0<=w[1:filter_length]<=1)
    for i=1:filter_length-1
        @constraint(model,w[i]>=w[i+1])
    end

    @objective(model,Min,sum((output[filter_length+1:end] - [filtro(input[i+1:i+filter_length],w) for i=1:length(input)-filter_length]).^2))

    optimize!(model)

    return value.(w)
end

function compute_active_weights(incidence::Vector{<:Number}, active::Vector{<:Number}; filter_length::Integer=30)
    return compute_filter_weights(incidence, active, filter_length)
end

function compute_icu_weights(incidence::Vector{<:Number}, icu_occupancy::Vector{<:Number}; filter_length::Integer=30)
    return compute_filter_weights(incidence, icu_occupancy, filter_length)
end
