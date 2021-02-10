module EpiSim

export  si_cori, si_geom, si_bn, si_covid, epi_estim_R,
        simulate, compute_active, compute_prediction,
        compute_harvard_levels,
        compute_active_weights,
        compute_icu_weights,
        simulate_mortality

include("si_distributions.jl")
include("epi_estim_cori.jl")
include("simulacion.jl")
include("utilities.jl")
include("fits.jl")

end
