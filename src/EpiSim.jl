module EpiSim

export si_cori,si_geom,si_bn,si_covid, epi_estim_R

include("si_distributions.jl")
include("epi_estim_cori.jl")

end
