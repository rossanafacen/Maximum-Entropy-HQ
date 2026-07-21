"""define a charm distribution function that depends on two lagrange multipliers  
"""
function f_ME(T, ur, eta, phi, etap, phip, pt, lm::Lagr_Multiplier_2D; m, q)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp = -ut*mt*cosh(etap-eta) + ur*pt*cos(phip-phi)
    p_transv = -ut*ur*mt*cosh(etap-eta) + (ut)^2*pt*cos(phip-phi)
    arg = udotp/T + q*lm.mult_n - q*lm.mult_nu*p_transv/udotp
    return exp(arg)
end



function charm_moments_integrand(T, ur, etap, phip, pt, lambda_n, lambda_nu, degeneracies, masses, charges; ccbar)
    
    eta = 0.0
    phi = 0.0
    ut = sqrt(1 + ur^2)
    total_density_integrand = 0.0
    total_diff_integrand = 0.0

    for (deg, m, q) in zip(degeneracies, masses, charges)
        #fact = Fluidum.canonical_factor(ccbar, q)
        fact = 1
        mt = sqrt(pt^2 + m^2)
        udotp = -ut*mt*cosh(etap-eta) + ur*pt*cos(phip - phi)
        p_transv = -ut*ur*mt*cosh(etap-eta) + ut^2*pt*cos(phip - phi)

        arg = udotp/T + q*lambda_n - q*lambda_nu*p_transv/udotp
        f = exp(arg)

        total_density_integrand += q*fact*deg*(-udotp)*f/(2π)^3 * pt
        total_diff_integrand += q*fact*deg*p_transv*f/(2π)^3 * pt
    end

    return SVector(total_density_integrand, total_diff_integrand)
end

function charm_moments(T, ur, lambda_n, lambda_nu, degeneracies, masses, charges; ccbar, etap_min = 0.0, etap_max = 10.0, phip_min = 0.0, phip_max = 2π, pt_min = 0.0, pt_max = 10.0, rtol = 1e-4)
    val = hcubature(b -> 2*fmGeV^3*charm_moments_integrand(T, ur, b[1], b[2], b[3], lambda_n, lambda_nu, degeneracies, masses, charges; ccbar = ccbar),
        (etap_min, phip_min, pt_min),
        (etap_max, phip_max, pt_max);
        rtol = rtol)[1]
    n0 = val[1]
    nur0 = val[2]
    return SVector{2}(n0, nur0)
end

function lagrangian_multipliers_system_2D(T, ur, n_hydro, nur_hydro, unknown, degeneracies, masses, charges; ccbar, rtol = 1e-4)
    lambda_n = unknown[1]
    lambda_nu = unknown[2]
    moments = charm_moments(T, ur, lambda_n, lambda_nu, degeneracies, masses, charges; ccbar, rtol = rtol)
    
    n_ME = moments[1]
    nur_ME = moments[2]
    return SVector(n_ME - n_hydro, nur_ME - nur_hydro)
end

function reconstruct_moments_2D(T, ur, lambda_n, lambda_nu, degeneracies, masses, charges; ccbar, rtol = 1e-4)
    moments = charm_moments(T, ur, lambda_n, lambda_nu, degeneracies, masses, charges; ccbar, rtol = rtol)
    density_ME = moments[1]
    diffusion_ME = moments[2]
    return density_ME, diffusion_ME
end


function lagrange_multiplier_2D(fo, params, degeneracies, masses, charges; ccbar, alpha_step)
    x = fo.x
    lb = Fluidum.leftbounds(x) 
    rb = Fluidum.rightbounds(x)
    alpha_array = range(lb[1], rb[1]; length = alpha_step)
    
    density_ME = zeros(length(alpha_array))
    diffusion_ME = copy(density_ME);

    lm = [Lagr_Multiplier_2D(0.0, 0.0) for _ in alpha_array]

    for (i, alpha) in enumerate(alpha_array)

        V = fo.fields(alpha)

        T = V[1]
        ur = V[2]
        fug = V[6]
        nur = V[7]
        n = densities(T, fug, params.eos.hq_eos).density

        if i == 1
            guess = SVector(fug, nur/T^3)
        else
            old_lm = lm[i-1]
            guess = SVector(old_lm.mult_n, old_lm.mult_nu)
        end
     
        _f(unknown, p) = MaximumEntropy.lagrangian_multipliers_system_2D(T, ur, n, nur, unknown, degeneracies, masses, charges; ccbar = ccbar, rtol = 1e-4)

        problem = NonlinearProblem{false}(_f, guess)

        lambda_sol = solve(problem, TrustRegion(); reltol = 1e-7)

        lm[i] = Lagr_Multiplier_2D(lambda_sol.u[1], lambda_sol.u[2])
        
        density_ME[i], diffusion_ME[i] = MaximumEntropy.reconstruct_moments_2D(T, ur, lambda_sol.u[1], lambda_sol.u[2], degeneracies, masses, charges; ccbar = ccbar, rtol = 1e-4)

        
        if abs(diffusion_ME[i] - nur) > 1e-5
            @warn "Starting diverging at alpha = $alpha"
        end
    end
    return lm, density_ME, diffusion_ME
end
