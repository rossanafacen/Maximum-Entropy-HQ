"""
Maximum Entropy distribution function expressed considering a symmetric 1+1 system, with energy and bulk constraints 
"""
function f_ME_LF(ur, eta, phi, etap, phip, pt, lm::Lagr_Multiplier_2D_LF; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta) + ur*pt*cos(phip-phi)
    pp = pt^2+mt^2/2*(cosh(2*etap-2*eta)-1)+ur^2*(pt^2*cos(phip-phi)^2+mt^2*cosh(etap-eta)^2)-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta) #p^⟨μ⟩ p_⟨μ⟩ 
    
    arg = lm.mult_e*udotp-lm.mult_P/udotp*pp
    return exp(arg)
end


function LF_moments_integrand_one_part(ur, etap, phip, pt, part, lm::Lagr_Multiplier_2D)
    m = part.Mass
    deg = part.Degeneracy
    eta = 0.0
    phi = 0.0
        
    f = f_ME_LF(ur, eta, phi, etap, phip, pt, lm; m = m, q = q)
    
    ut = sqrt(1 + ur^2)
    mt = sqrt(pt^2 + m^2)

    udotp = -ut*mt*cosh(etap-eta) + ur*pt*cos(phip - phi)
    pp = pt^2+mt^2/2*(cosh(2*etap-2*eta)-1)+ur^2*(pt^2*cos(phip-phi)^2+mt^2*cosh(etap-eta)^2)-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)  
    

    energy_integrand = udotp^2*f/(2*pi)^3*pt
    pressure_integrand = 1/3*pp*f/(2*pi)^3*pt #+1/3 because of our mostly-plus metric 


    return SVector(energy_integrand, pressure_integrand)
end


function LF_moments_one_part(ur, part, lm::Lagr_Multiplier_2D_LF; etap_min = 0.0, etap_max = 10.0, phip_min = 0.0, phip_max = 2π, pt_min = 0.0, pt_max = 10.0, rtol = 1e-4)
        val = hcubature(b -> 2*fmGeV^3*LF_moments_integrand_one_part(ur, b[1], b[2], b[3], part, lm; ccbar = ccbar),
        (etap_min, phip_min, pt_min),
        (etap_max, phip_max, pt_max);
        rtol = rtol)[1]

    epsilon = val[1]
    P_tot = val[2]

    return SVector{2}(epsilon, P_tot)
end


