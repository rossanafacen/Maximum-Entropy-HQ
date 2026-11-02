"""
Maximum Entropy distribution function expressed considering a symmetric 1+1 system, with energy and bulk constraints 
"""
function f_ME_LF(eta, etap, pt, lm::Lagr_Multiplier_2D_LF; m)
    mt = sqrt(pt^2+m^2)
    p0 = mt*cosh(etap-eta) 
    p2 = p0^2-m^2 #p = sqrt(E^2 - m^2)
    arg = -lambda_e*p0 + lambda_P*p2/p0
    return exp(arg)
end


function LF_moments_integrand(etap, pt, lambda_e, lambda_P, degeneracies, masses)
    eta = 0.0
    total_energy_integrand = 0.0
    total_pressure_integrand = 0.0

    for (deg, m) in zip(degeneracies, masses)
        mt = sqrt(pt^2 + m^2)
        p0 = mt*cosh(etap-eta) 
        p2 = p0^2-m^2 #p = sqrt(E^2 - m^2)
    
        arg = -lambda_e*p0 + lambda_P/p0*p2
        f = exp(arg)

        total_energy_integrand += deg*p0^2*f/(2π)^3 * pt #GeV3
        total_pressure_integrand += deg*1/3*p2*f/(2π)^3 * pt #GeV3
    end

    return SVector(total_energy_integrand, total_pressure_integrand)
end

function LF_moments(lambda_e, lambda_P, degeneracies, masses; etap_min = 0.0, etap_max = 10.0, pt_min = 0.0, pt_max = 10.0, rtol = 1e-4)
    val = hcubature(b -> 2*2π*LF_moments_integrand(b[1], b[2], lambda_e, lambda_P, degeneracies, masses),
        (etap_min, pt_min),
        (etap_max, pt_max);
        rtol = rtol)[1]
    epsilon0 = val[1]
    P0 = val[2]
    return SVector{2}(epsilon0, P0)
end

function lagrangian_multipliers_system_2D(epsilon_hydro, P_hydro, unknown, degeneracies, masses; rtol = 1e-4)
    lambda_e = unknown[1]
    lambda_P = unknown[2]
    moments = LF_moments(lambda_e, lambda_P, degeneracies, masses; etap_min = 0.0, etap_max = 10.0, pt_min = 0.0, pt_max = 10.0, rtol = rtol)
    
    epsilon_ME = moments[1]
    P_ME = moments[2]
    return SVector(epsilon_ME - epsilon_hydro, P_ME - P_hydro)
end

function reconstruct_moments_2D(lambda_e, lambda_P, degeneracies, masses; rtol = 1e-4)
    moments = LF_moments(lambda_e, lambda_P, degeneracies, masses; etap_min = 0.0, etap_max = 10.0, pt_min = 0.0, pt_max = 10.0, rtol = rtol)
    epsilon_ME = moments[1]
    P_ME = moments[2]
    return epsilon_ME, P_ME
end


function lagrange_multiplier_2D_LF_LRF(fo, hrg::Fluidum.HadronResonanceGas_LF; alpha_step)
    
    degeneracies = hrg.particle_list.Degeneracy
    masses = hrg.particle_list.Mass

    x = fo.x
    lb = Fluidum.leftbounds(x) 
    rb = Fluidum.rightbounds(x)
    alpha_array = range(lb[1], rb[1]; length = alpha_step)
    
    epsilon_ME = zeros(length(alpha_array))
    P_ME = copy(epsilon_ME);

    lm = [Lagr_Multiplier_2D_LF(0.0, 0.0) for _ in alpha_array]

    for (i, alpha) in enumerate(alpha_array)

        V = fo.fields(alpha)

        T = V[1]
        piB = V[5]/fmGeV^3
        
        thermo = Fluidum.thermodynamic(T,hrg)
        
        epsilon = T*thermo.pressure_derivative[1] - thermo.pressure
        P = thermo.pressure + piB

        if i == 1
            guess = SVector(1/T, 0.0)
        else
            old_lm = lm[i-1]
            guess = SVector(old_lm.mult_e, old_lm.mult_P)
        end
     
        _f(unknown, p) = MaximumEntropy.lagrangian_multipliers_system_2D(epsilon, P, unknown, degeneracies, masses; rtol = 1e-4)

        problem = NonlinearProblem{false}(_f, guess)

        lambda_sol = solve(problem, TrustRegion(); reltol = 1e-7)

        lm[i] = Lagr_Multiplier_2D_LF(lambda_sol.u[1], lambda_sol.u[2])
        
        epsilon_ME[i], P_ME[i] = MaximumEntropy.reconstruct_moments_2D(lambda_sol.u[1], lambda_sol.u[2], degeneracies, masses; rtol = 1e-4)

        
        if abs(P_ME[i] - P) > 1e-5
            @warn "Starting diverging at alpha = $alpha"
        end
    end
    return lm, epsilon_ME, P_ME
end
