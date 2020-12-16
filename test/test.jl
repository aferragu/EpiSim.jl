using Plots, CSV, HTTP, JSON, Dates, DataFrames, LaTeXStrings, EpiSim
default(size=(1000,500))
#GUIAD
req = HTTP.get("https://raw.githubusercontent.com/GUIAD-COVID/datos-y-visualizaciones-GUIAD/master/datos/estadisticasUY.csv")
dataguiad = CSV.read(IOBuffer(req.body),missingstring="N/A", DataFrame)
dataguiad[!,:fecha] = Dates.Date.(dataguiad[!,:fecha],"dd/mm/yyyy")

fechauy = dataguiad[!,:fecha]
Iuy = dataguiad[!,:cantCasosNuevos]; Iuy[1] = 4; Iuy = max.(Iuy,0);
activos_uy = dataguiad[!,:cantPersonasConInfeccionEnCurso];

cti_uy = dataguiad[!,:cantCTI];
cti_uy[ismissing.(cti_uy)].=0
cti_uy=collect(skipmissing(cti_uy));


ventana=7
start = Date(2020,11,15)
startdate = Date(2020,12,1)
incidence_0 = Iuy[fechauy.<=startdate]
#filtro una pasada mas del R
R0,Rl,Ru,a0,b0, Lambda = epi_estim_R(incidence_0)
Rfuturos = R0[end-ventana+1:end]

dias = 30
reps = 10000

median,lower,upper = calculo_prediccion(Rfuturos,incidence_0,dias,reps, si_covid);

verde, amarillo, naranja, rojo = compute_harvard_levels(3505984)

p=plot(legend=:topleft)
bar!(p,fechauy[fechauy.>=start],Iuy[fechauy.>=start], alpha=0.6, label="Incidencia observada UY", xticks = start:Dates.Day(2):Date(2020,12,31))

fechas = start:Dates.Day(1):startdate+Dates.Day(dias)
plot!(fechas, ones(size(fechas))*verde, fillrange=0, color=:green, fillalpha=0.1, label=:none)
plot!(fechas, ones(size(fechas))*amarillo, fillrange=verde, color=:yellow, fillalpha=0.1, label=:none)
plot!(fechas, ones(size(fechas))*naranja, fillrange=amarillo, color=:orange, fillalpha=0.1, label=:none)
plot!(fechas, ones(size(fechas))*rojo, fillrange=naranja, color=:red, fillalpha=0.1, label=:none)

plot!(p,(startdate+Dates.Day(1):Dates.Day(1):startdate+Dates.Day(dias)),median, label="Proyeccion al $(startdate)", ribbon=(median-lower,upper-median),lw=2, color=:green, fillalpha=0.2)
plot!(p,title="Incidencia observada y proyección - Todo el país")
plot!(xformatter = x -> Dates.format(Date(Dates.UTD(x)), "dd/mm"))
