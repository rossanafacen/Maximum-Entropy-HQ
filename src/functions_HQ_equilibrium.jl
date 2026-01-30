struct Lagr_Multiplier_1D_eq{A<:Real} <: Lagr_Multiplier
    mult_n::A
end

struct Lagr_Multiplier_2D_eq{A<:Real, B<:Real} <: Lagr_Multiplier
    mult_e::A
    mult_n::B
end


"""charm distribution function that depends on one lagrange multiplier, the fugacity 
"""
function f_ME(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D_eq; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi) #negative defined
    arg = udotp/T+lm.mult_n  

    return exp(arg)
end


"""charm quark density, defined with negative sign due to the (-,+,+,+) metric"""
function charm_density_integrand(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D_eq; m)
    
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(T,ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt

end



"""integrate the distribution function over the momentum space, to get the charm quark fields. 
Factor 2 comes from the rapidity integration, made only for positive rapidity.
In the integration of the distribution function, it is possible to impose eta_p = and phi_p = 0, since a shift for these coordinates does not change the integral result
"""
function charm_density(T,ur,lm::Lagr_Multiplier_1D_eq;m,ccbar,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0,pt_max=8.0,rtol=10E-4) 
    eta = 0 
    phi = 0
    fact = Fluidum.besseli(1, ccbar/2)./Fluidum.besseli(0, ccbar/2)
    hcubature( b->fact*2*fmGeV^3*charm_density_integrand(T,ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

function lagrangian_multipliers_system_eq(T,ur,n,unknown;m,ccbar)
    lm = Lagr_Multiplier_1D_eq(unknown)
    eq_charm = charm_density(T,ur,lm;m=m,ccbar = ccbar)-n 
    return (eq_charm)
end

#####consider two lagrange multipliers now

"""charm distribution function that depends on two lagrange multiplier: the fugacity and the temperature
"""
function f_ME(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D_eq; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi) #negative defined
    arg = lm.mult_e*udotp+lm.mult_n  
    return exp(arg)
end

"""
Integrand of the energy density in p_T*dp_T, as defined in kinetic theory
"""
function energy_integrand(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D_eq; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return udotp^2*f_ME(ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt
end

"""charm quark density, defined with negative sign due to the (-,+,+,+) metric"""
function charm_density_integrand(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D_eq; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt

end

"""
Energy density integrate p_T*dp_T, as defined in kinetic theory
"""
function charm_energy_density(ur,lm::Lagr_Multiplier_2D_eq;m,ccbar,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0.,pt_max=8.0,rtol=1e-8)
    eta = 0
    phi = 0 
    fact = Fluidum.besseli(1, ccbar/2)./Fluidum.besseli(0, ccbar/2)
    hcubature(b->fact*2*fmGeV^3*energy_integrand(ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

"""integrate the distribution function over the momentum space, to get the charm quark fields. 
Factor 2 comes from the rapidity integration, made only for positive rapidity.
In the integration of the distribution function, it is possible to impose eta_p = and phi_p = 0, since a shift for these coordinates does not change the integral result
"""
function charm_density(ur,lm::Lagr_Multiplier_2D_eq;m,ccbar,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0,pt_max=8.0,rtol=1e-8) 
    eta = 0 
    phi = 0
    fact = Fluidum.besseli(1, ccbar/2)./Fluidum.besseli(0, ccbar/2)
    hcubature( b->fact*2*fmGeV^3*charm_density_integrand(ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

"""lagrange multipliers system in 2D
"""
function lagrangian_multipliers_system_eq_2d(ur,epsilon,n,unknown;m,ccbar)
    lm = Lagr_Multiplier_2D_eq(unknown[1],unknown[2])
    eq_charm1 = charm_energy_density(ur,lm;m=m, ccbar=ccbar)-epsilon 
    eq_charm2 = charm_density(ur,lm;m=m,ccbar = ccbar)-n 
    return SVector{2}(eq_charm1, eq_charm2)
end