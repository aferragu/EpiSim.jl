module EpiSim

export si_cori,si_geom,si_bn,si_covid, epi_estim_R, simula_Epiestim, calculo_prediccion, compute_harvard_levels

include("si_distributions.jl")
include("epi_estim_cori.jl")
include("simulacion.jl")
include("utilities.jl")

end
