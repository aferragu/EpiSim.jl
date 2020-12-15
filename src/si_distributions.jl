## Defines serial interval distributions

using Distributions

#shifted gamma from Cori et al.
function shifted_gamma(k::Integer,a::Real,b::Real)

    @assert a>=0
    @assert b>=0

    res = k * cdf(Gamma(a,b),k) +  (k - 2) * cdf(Gamma(a,b),k-2) - 2 * (k - 1) * cdf(Gamma(a,b),k - 1)

    res = res + a * b * (2 * cdf(Gamma(a+1,b),k - 1) -  cdf(Gamma(a + 1, b),k - 2) - cdf(Gamma(a+1,b),k))

    return res

end

#Serial Interval from Cori et al.
function si_cori(k::Integer, mu::Real,sigma::Real)

    @assert mu>=0
    @assert sigma>0
    @assert k>=0

    a = ((mu - 1) / sigma)^2
    b = sigma^2 / (mu - 1)

    return shifted_gamma(k,a,b)

end

#Geometric
function si_geom(k::Integer, mean::Real)

    @assert mean>0
    @assert k>=0

    p=1.0/mean
    dist=Geometric(p)
    return pdf(dist,k-1)

end

#Negative Binomial
function si_bn(k::Integer,mean::Real,sigma::Real)

    @assert mean>=0
    @assert sigma>0

    p = mean/sigma
    r = mean*p/(1-p)
    dist=NegativeBinomial(r,p)
    return pdf(dist,k-1)
end

#Covid distribution
function si_covid(k::Integer, mean::Real=3.95, sigma::Real=4.75)

    return si_cori(k,mean,sigma)

end
