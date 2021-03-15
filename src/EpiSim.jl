module EpiSim

export  si_cori, si_geom, si_bn, si_covid, epi_estim_R,
        simulate, compute_prediction,
        compute_harvard_levels,
        compute_active_weights,
        compute_icu_weights,
        compute_active_fixed_length,
        simulate_active,
        simulate_mortality,
        simulate_icu,
        compute_quantiles


include("si_distributions.jl")
include("epi_estim_cori.jl")
include("simulation.jl")
include("utilities.jl")
include("fits.jl")

end
