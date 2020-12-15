using DSP

function epi_estim_R(incidence; window::Integer=7, si_distr::Function=si_covid, T::Integer=30, a0::Real=1, b0::Real=5)

    @assert a0>=0
    @assert b0>0

    w=si_distr.((0:T));
    Lambda=filt(w,[1.0],incidence)

    n=length(incidence)

    a = a0 .+ filt(ones(window), incidence)
    b = 1.0./(1/b0 .+ filt(ones(window),[1.0], Lambda))
    R = a.*b

    Rl = [quantile(Gamma(a[i],b[i]),0.025) for i in 1:n]
    Ru = [quantile(Gamma(a[i],b[i]),0.975) for i in 1:n]


    return R, Rl,Ru, a, b, Lambda
end;
