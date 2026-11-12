"""
Maximum Entropy distribution function expressed considering a symmetric 1+1 system, with energy and bulk constraints
In this case, since we are using this for freeze-out, we cannot consider the lrf     
"""
function f_ME(ur, eta, phi, etap, phip, pt, lm::Lagr_Multiplier_1D_LF; m, theta)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta) + ur*pt*cos(phip-phi)
    
    arg = lm.mult_e*udotp
    z = exp(arg)
    return z/(1 + theta*z)
end


function LF_moments_integrand(etap, pt, lambda_e, degeneracies, masses)
    eta = 0.0
    total_energy_integrand = 0.0

    for (deg, m) in zip(degeneracies, masses)
        theta = (-1)^deg #bosons have odd degeneracy, fermions have even degeneracy
        mt = sqrt(pt^2 + m^2)
        p0 = mt*cosh(etap-eta) 
        
        arg = -lambda_e*p0     
        z = exp(arg)
        f = z/(1 + theta*z)
        total_energy_integrand += deg*p0^2*f/(2π)^3 * pt #GeV3
    end

    return SVector(total_energy_integrand)
end

function LF_moments(lambda_e, degeneracies, masses; etap_min = 0.0, etap_max = 10.0, pt_min = 0.0, pt_max = 10.0, rtol = 1e-4)
    val = hcubature(b -> 2*2π*LF_moments_integrand(b[1], b[2], lambda_e, degeneracies, masses),
        (etap_min, pt_min),
        (etap_max, pt_max);
        rtol = rtol)[1]
    epsilon0 = val[1]
    return epsilon0
end

function lagrangian_multipliers_system_1D(epsilon_hydro, unknown, degeneracies, masses; rtol = 1e-4)
    lambda_e = unknown[1]
    moments = LF_moments(lambda_e, degeneracies, masses; etap_min = 0.0, etap_max = 10.0, pt_min = 0.0, pt_max = 10.0, rtol = rtol)
    
    epsilon_ME = moments[1]
    #return SVector(epsilon_ME - epsilon_hydro)
    return (epsilon_ME/epsilon_hydro - 1)
end

function reconstruct_moments_1D(lambda_e, degeneracies, masses; rtol = 1e-7)
    moments = LF_moments(lambda_e, degeneracies, masses; etap_min = 0.0, etap_max = 10.0, pt_min = 0.0, pt_max = 10.0, rtol = rtol)
    epsilon_ME = moments[1]
    return epsilon_ME
end

function lagrange_multiplier_1D_LF_LRF(fo, hrg::Fluidum.HadronResonanceGas_LF, params; alpha_step)
    
    degeneracies = hrg.particle_list.Degeneracy
    masses = hrg.particle_list.Mass

    x = fo.x
    lb = Fluidum.leftbounds(x) 
    rb = Fluidum.rightbounds(x)
    alpha_array = range(lb[1], rb[1]; length = alpha_step)
    
    epsilon_ME = zeros(length(alpha_array))

    lm = [Lagr_Multiplier_1D_LF(0.0) for _ in alpha_array]

    for (i, alpha) in enumerate(alpha_array)

        V = fo.fields(alpha)
        T = V[1]
        
        thermo = Fluidum.thermodynamic(T, params.eos)
        epsilon = (T*thermo.pressure_derivative[1])/fmGeV^3 - thermo.pressure/fmGeV^3
        
        if i == 1
            guess = 1/T
        else
            old_lm = lm[i-1]
            guess = old_lm.mult_e
        end
     
        _f(unknown, p) = MaximumEntropy.lagrangian_multipliers_system_1D(epsilon, unknown, degeneracies, masses; rtol = 1e-4)

        problem = NonlinearProblem{false}(_f, guess)

        lambda_sol = solve(problem, TrustRegion(); abstol = 1e-6, reltol = 1e-7)
        lm[i] = Lagr_Multiplier_1D_LF(lambda_sol.u[1])
            
        epsilon_ME[i] = MaximumEntropy.reconstruct_moments_1D(lambda_sol.u[1], degeneracies, masses; rtol = 1e-4)
        
        if abs(epsilon_ME[i] - epsilon) > 1e-5
            @warn "Starting diverging at alpha = $alpha"
        end
    end
    return lm, epsilon_ME
end


