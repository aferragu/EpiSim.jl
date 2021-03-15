#utilities

function compute_harvard_levels(population::Number)

    pop = population/100000
    levels = (1,10,25,45)

    return levels.*pop

end

function compute_quantiles(result::Array{<:Number,2})

    n=size(result,2)

    med=zeros(n)
    low=zeros(n)
    up=zeros(n)

    for j=1:n
        res = quantile(result[:,j],[0.025,0.5,0.975])
        low[j] =res[1]
        med[j] =res[2]
        up[j] =res[3]
    end

    return med,low,up
end
