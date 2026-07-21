function charm_moments_integrand_one_part(T, ur, etap, phip, pt, part, lm::Lagr_Multiplier_2D; ccbar = 1)
    m = part.Mass
    q = 1
    deg = part.Degeneracy
    eta = 0.0
    phi = 0.0
    
    fact = 1
    
    f = f_ME(T, ur, eta, phi, etap, phip, pt, lm; m = m, q = q)
    
    ut = sqrt(1 + ur^2)
    mt = sqrt(pt^2 + m^2)

    udotp = -ut*mt*cosh(etap-eta) + ur*pt*cos(phip - phi)
    p_transv = -ut*ur*mt*cosh(etap-eta) + ut^2*pt*cos(phip - phi)
    

    density_integrand = q*fact*deg*(-udotp)*f/(2π)^3 * pt
    diff_integrand = q*fact*deg*p_transv*f/(2π)^3 * pt


    return SVector(density_integrand, diff_integrand)
end

function charm_moments_one_part(T, ur, part, lm::Lagr_Multiplier_2D; ccbar = 1, etap_min = 0.0, etap_max = 10.0, phip_min = 0.0, phip_max = 2π, pt_min = 0.0, pt_max = 10.0, rtol = 1e-4)
    

    val = hcubature(b -> 2*fmGeV^3*charm_moments_integrand_one_part(T, ur, b[1], b[2], b[3], part, lm; ccbar = ccbar),
        (etap_min, phip_min, pt_min),
        (etap_max, phip_max, pt_max);
        rtol = rtol)[1]

    n0 = val[1]
    nur0 = val[2]

    return SVector{2}(n0, nur0)
end

