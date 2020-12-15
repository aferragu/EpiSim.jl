#simulation code
using Distributions

function simula_Epiestim(R::Vector{<:Real},incidence_0::Vector{<:Number},dias::Integer, si_distr::Function; T=30)

    periodo=length(R)

    pesos = si_distr.((1:T))
    h = length(pesos)

    m=length(incidence_0)

    incidence=zeros(dias+m)
    Lambda = zeros(dias+m)

    incidence[1:m] = incidence_0;

    for i=m:dias+m-1
        idx = i%periodo + 1
        Lambda[i+1] = sum(incidence[i:-1:max(i-h+1,1)].*pesos[1:min(i,h)])
        incidence[i+1] = rand(Poisson(R[idx]*Lambda[i+1]))
    end

    return incidence
end

function calculo_prediccion(R::Vector{<:Real},incidence_0::Vector{<:Real},dias::Integer,reps::Integer, si_distr::Function; T=30)

    incidencias = zeros(reps, dias)

    for j=1:reps
        incidencia = simula_Epiestim(R,incidence_0,dias, si_distr; T=T)
        incidencias[j,:] = incidencia[end-dias+1:end]
    end

    median = zeros(dias)
    lower = zeros(dias)
    upper = zeros(dias)

    for i=1:dias
        median[i] = quantile(incidencias[:,i],0.5)
        lower[i] = quantile(incidencias[:,i],0.025)
        upper[i] = quantile(incidencias[:,i],0.975)
    end

    return median,lower,upper

end
