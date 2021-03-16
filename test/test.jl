using Plots, CSV, HTTP, Dates, DataFrames, StatsBase, LaTeXStrings, EpiSim
default(size=(1000,500))
const UYPOP = 3505984

#GUIAD
req = HTTP.get("https://raw.githubusercontent.com/GUIAD-COVID/datos-y-visualizaciones-GUIAD/master/datos/estadisticasUY.csv")
dataguiad = CSV.read(IOBuffer(req.body),missingstring="N/A", DataFrame)
dataguiad[!,:fecha] = Dates.Date.(dataguiad[!,:fecha],"dd/mm/yyyy")

fecha = dataguiad[!,:fecha]
incidencia = dataguiad[!,:cantCasosNuevos]; incidencia[1] = 4; incidencia = max.(incidencia,0);
activos = dataguiad[!,:cantPersonasConInfeccionEnCurso];

cti = dataguiad[!,:cantCTI];
cti[ismissing.(cti)].=0
cti=collect(skipmissing(cti));


ventana=7
start = Date(2021,2,21)
start_date = Date(2021,3,8)
incidence_0 = incidencia[fecha.<=start_date]
#filtro una pasada mas del R con ventana 14 dias
R0,Rl,Ru,a0,b0, Lambda = epi_estim_R(incidence_0, window=14)
Rfuturos = R0[end-ventana+1:end]

println("R efectivo= $(geomean(Rfuturos))")

end_date = Date(2021,3,31)
dias = Dates.value(end_date-start_date)

reps = 10000

incidencias = zeros(Int64,reps, dias)

for j=1:reps
    incidencia_sim = simulate(Rfuturos,incidence_0,dias, EpiSim.si_covid)
    incidencias[j,:] = incidencia_sim[end-dias+1:end]
end

Imedian,Ilower,Iupper = compute_quantiles(incidencias);

verde, amarillo, naranja, rojo = compute_harvard_levels(UYPOP)

p=plot(legend=:topleft)
bar!(p,fecha[fecha.>=start],incidencia[fecha.>=start], alpha=0.6, label="Incidencia observada UY", xticks = start:Dates.Day(2):Date(2021,4,1))

fechas = start:Dates.Day(1):start_date+Dates.Day(dias)
plot!(fechas, ones(size(fechas))*verde, fillrange=0, color=:green, fillalpha=0.1, label=:none)
plot!(fechas, ones(size(fechas))*amarillo, fillrange=verde, color=:yellow, fillalpha=0.1, label=:none)
plot!(fechas, ones(size(fechas))*naranja, fillrange=amarillo, color=:orange, fillalpha=0.1, label=:none)
plot!(fechas, ones(size(fechas))*rojo, fillrange=naranja, color=:red, fillalpha=0.1, label=:none)

plot!(p,(start_date+Dates.Day(1):Dates.Day(1):start_date+Dates.Day(dias)),Imedian, label="Proyeccion al $(start_date)", ribbon=(Imedian-Ilower,Iupper-Imedian),lw=2, color=:green, fillalpha=0.2, xticks=start:Dates.Day(4):start_date+Dates.Day(dias))
plot!(p,title="Incidencia observada y proyección - Todo el país")
plot!(xformatter = x -> Dates.format(Date(Dates.UTD(x)), "dd/mm"))


activos_simulados = zeros(reps, dias)

pact = zeros(20)
pact[10]=1.0;

for j=1:reps
    activos_sim = simulate_active(incidencias[j,:], pact)
    activos_simulados[j,:] = activos_sim[end-dias+1:end]
end

Amedian,Alower,Aupper = compute_quantiles(activos_simulados)

plot(Amedian)
plot!(Alower)
plot!(Aupper)

cti_simulados = zeros(reps, dias)

pact = zeros(20)
pact[10]=1.0;

for j=1:reps
    cti_sim = simulate_icu(incidencias[j,:],0.01, pact)
    activos_simulados[j,:] = cti_sim[end-dias+1:end]
end

CTImedian,CTIlower,CTIupper = compute_quantiles(activos_simulados)

plot(CTImedian)
plot!(CTIlower)
plot!(CTIupper)
