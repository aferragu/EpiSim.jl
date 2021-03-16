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

function compute_active_fixed_length(incidence::Array{<:Real}, periodo_actividad::Integer)

    activos = [sum(incidence[i:-1:max(i-periodo_actividad+1,1)]) for i=1:length(incidence)]
    return activos

end

function simulate_active(incidence::Array{<:Real}, activity_distribution::Vector{<:Number})

    dias = length(incidence)

    activity = DiscreteNonParametric(collect(1:length(activity_distribution)),activity_distribution/sum(activity_distribution))

    activos = zeros(dias)

    for i=1:dias

        m = rand(activity,incidence[i])

        for j=1:length(m)
            for k=i:min(i+m[j],dias)
                activos[k] = activos[k]+1
            end
        end
    end

    return activos

end

function simulate_mortality(incidence::Vector{<:Number}, death_prob::Real, mortality_weights::Vector{<:Number})

    dias = length(incidence)

    mortality_dist = DiscreteNonParametric(collect(1:length(mortality_weights)),mortality_weights/sum(mortality_weights))

    mortality = zeros(dias)

    for i=1:dias

        num_deaths = rand(Binomial(incidence[i],death_prob))

        m = rand(mortality_dist,num_deaths)

        for j=1:length(m)
            if i+m[j]<=dias
                mortality[i+m[j]] = mortality[i+m[j]]+1
            end
        end
    end

    return mortality

end

function simulate_icu(incidence::Vector{<:Number}, icu_prob::Real, icu_weights::Vector{<:Number})

    dias = length(incidence)

    icu_distr = DiscreteNonParametric(collect(1:length(icu_weights)),icu_weights/sum(icu_weights))

    icu = zeros(dias)

    for i=1:dias

        new_icu = rand(Binomial(incidence[i],icu_prob))


        m = rand(icu_distr,new_icu)

        for j=1:length(m)
            for k=i:min(i+m[j],dias)
                icu[k] = icu[k]+1
            end
        end
    end

    return icu

end

function compute_prediction(R::Vector{<:Real},incidence_0::Vector{<:Real},dias::Integer,reps::Integer, si_distr::Function; T=30)

    incidencias = zeros(reps, dias)
    activos = zeros(reps, dias)

    for j=1:reps
        incidencia = simulate(R,incidence_0,dias, si_distr; T=T)
        incidencias[j,:] = incidencia[end-dias+1:end]

        activo = compute_active_fixed_length(incidencia, 12)
        activos[j,:] = activo[end-dias+1:end]

    end

    medianI,lowerI,upperI =  compute_quantiles(incidencias)
    medianA, lowerA, upperA = compute_quantiles(activos)

    return medianI,lowerI,upperI, medianA, lowerA, upperA

end
