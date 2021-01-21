#simulation code
using Distributions

function simulate(R::Vector{<:Real},incidence_0::Vector{<:Number},dias::Integer, si_distr::Function; T=30)

    periodo=length(R)

    pesos = si_distr.(1:T)
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

function compute_active(incidence::Array{<:Real}, activos_iniciales::Array{<:Real}, periodo_actividad::Integer)

    activos = sum(incidence[i:-1:max(i-periodo_actividad+1,1)])
    return activos

end

function compute_prediction(R::Vector{<:Real},incidence_0::Vector{<:Real},dias::Integer,reps::Integer, si_distr::Function; T=30)

    incidencias = zeros(reps, dias)
    activos = zeros(reps, dias)

    for j=1:reps
        incidencia = simulate(R,incidence_0,dias, si_distr; T=T)
        incidencias[j,:] = incidencia[end-dias+1:end]

        activo = compute_active(incidencia, [0.0], 12)
        activos[j,:] = activo[end-dias+1:end]

    end

    medianI = zeros(dias)
    lowerI = zeros(dias)
    upperI = zeros(dias)

    medianA = zeros(dias)
    lowerA = zeros(dias)
    upperA = zeros(dias)

    for i=1:dias
        medianI[i] = quantile(incidencias[:,i],0.5)
        lowerI[i] = quantile(incidencias[:,i],0.025)
        upperI[i] = quantile(incidencias[:,i],0.975)

        medianA[i] = quantile(activos[:,i],0.5)
        lowerA[i] = quantile(activos[:,i],0.025)
        upperA[i] = quantile(activos[:,i],0.975)

    end

    return medianI,lowerI,upperI, medianA, lowerA, upperA

end
