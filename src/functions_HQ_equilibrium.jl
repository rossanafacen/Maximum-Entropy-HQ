abstract type Lagr_Multiplier end

struct Lagr_Multiplier_1D{A<:Real} <: Lagr_Multiplier
    mult_nu::A
end


function Base.zero(::Type{MaximumEntropy.Lagr_Multiplier_1D{Float64}})
    MaximumEntropy.Lagr_Multiplier_1D(0.0)
end

function Base.:+(a::Lagr_Multiplier_1D{Float64}, b::Real) 
    MaximumEntropy.Lagr_Multiplier_1D(a.mult_nu + float(b))
end 

"""charm distribution function that depends on one lagrange multiplier 
"""
function f_ME(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi) #negative defined
    arg = udotp/T+lm.mult_nu  

    return exp(arg)
end


"""charm quark density, defined with negative sign due to the (-,+,+,+) metric"""
function charm_density_integrand(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D; m)
    
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(T,ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt

end



"""integrate the distribution function over the momentum space, to get the charm quark fields. 
Factor 2 comes from the rapidity integration, made only for positive rapidity.
In the integration of the distribution function, it is possible to impose eta_p = and phi_p = 0, since a shift for these coordinates does not change the integral result
"""
function charm_density(T,ur,lm::Lagr_Multiplier_1D;m,ccbar,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0,pt_max=8.0,rtol=10E-4) 
    eta = 0 
    phi = 0
    fact = Fluidum.besseli(1, ccbar/2)./Fluidum.besseli(0, ccbar/2)
    hcubature( b->fact*2*fmGeV^3*charm_density_integrand(T,ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

function lagrangian_multipliers_system(T,ur,density,unknown;m,ccbar)
    lm = Lagr_Multiplier_1D(unknown)
    eq_charm = charm_density(T,ur,lm;m=m,ccbar = ccbar)[1]-density 
    return (eq_charm)
end



