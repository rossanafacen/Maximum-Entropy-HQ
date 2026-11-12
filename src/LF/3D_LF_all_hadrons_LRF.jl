####fix also the shear viscosity
"""
Maximum Entropy distribution function expressed considering a symmetric 1+1 system, with energy, bulk and shear constraints 
In this case, since we are using this for freeze-out, we cannot consider the lrf    
"""
function f_ME(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_3D_LF; m, theta)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp = -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    
    pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2*cos(phip-phi)^2+mt^2*cosh(etap-eta)^2)-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)  
    
    g1 = lm.mult_piϕϕ
    g2 = lm.mult_piηη

    gpp = -(g1+g2)*(ur*mt*cosh(etap-eta)-ut*pt*cos(phip-phi))^2 + g1*pt^2*sin(phip - phi)^2 + g2*mt^2*sinh(etap - eta)^2 #p^⟨μ⟩ p^⟨ν⟩ γ_⟨μν⟩
    arg = lm.mult_e*udotp - lm.mult_P/udotp*pp + gpp/udotp
    z = exp(arg)
    return z/(1 + theta*z)
end




function LF_moments_integrand(etap, phip, pt, lambda_e, lambda_P, lambda_piϕϕ, lambda_piηη, degeneracies, masses)
    
    eta = 0.0
    phi = 0.0
    total_energy_integrand = 0.0
    total_pressure_integrand = 0.0
    total_shearϕϕ_integrand = 0.0
    total_shearηη_integrand = 0.0
    
    g1 = lambda_piϕϕ
    g2 = lambda_piηη
    
    for (deg, m) in zip(degeneracies, masses)
        theta = (-1)^deg #bosons have odd degeneracy, fermions have even degeneracy
        
        mt = sqrt(pt^2 + m^2)
        p0 = mt*cosh(etap-eta) 
        p2 = pt^2 + mt^2*sinh(etap-eta)^2 #
        gpp = -(g1+g2)*(pt*cos(phip-phi))^2 + g1*pt^2*sin(phip - phi)^2 + g2*mt^2*sinh(etap - eta)^2 #p^⟨μ⟩ p^⟨ν⟩ γ_⟨μν⟩
    
        pp_ϕϕ = 1/6*(mt^2 + pt^2 - 3*pt^2 *cos(2phi - 2phip) - mt^2 *cosh(2eta - 2etap))
        
        
        pp_ηη = -1/3*(mt^2 + pt^2 - mt^2 *cosh(2eta - 2etap))      

        arg = -lambda_e*p0 + lambda_P*p2/p0 - gpp/p0
        z = exp(arg)
        f = z/(1 + theta*z)

        total_energy_integrand += deg*p0^2*f/(2π)^3 * pt #GeV3
        total_pressure_integrand += deg*1/3*p2*f/(2π)^3 * pt #GeV3
        total_shearϕϕ_integrand += deg*pp_ϕϕ*f/(2π)^3 * pt
        total_shearηη_integrand += deg*pp_ηη*f/(2π)^3 * pt
    end

    return SVector(total_energy_integrand, total_pressure_integrand, total_shearϕϕ_integrand, total_shearηη_integrand)
end

function LF_moments(lambda_e, lambda_P, lambda_piϕϕ, lambda_piηη, degeneracies, masses; etap_min = 0.0, etap_max = 10.0, phip_min = 0.0, phip_max = 2π, pt_min = 0.0, pt_max = 10.0, rtol = 1e-4)
    val = hcubature(b -> 2*LF_moments_integrand(b[1], b[2], b[3], lambda_e, lambda_P, lambda_piϕϕ, lambda_piηη, degeneracies, masses),
        (etap_min, phip_min, pt_min),
        (etap_max, phip_max, pt_max);
        rtol = rtol)[1]
    epsilon0 = val[1]
    P0 = val[2]
    piϕϕ0 = val[3]
    piηη0 = val[4]
    return SVector{4}(epsilon0, P0, piϕϕ0, piηη0)
end

function lagrangian_multipliers_system_3D(epsilon_hydro, P_hydro, piϕϕ_hydro, piηη_hydro, unknown, degeneracies, masses; rtol = 1e-4)
    lambda_e = unknown[1]
    lambda_P = unknown[2]
    lambda_piϕϕ = unknown[3]
    lambda_piηη = unknown[4]
    moments = LF_moments(lambda_e, lambda_P, lambda_piϕϕ, lambda_piηη, degeneracies, masses; etap_min = 0.0, etap_max = 10.0, phip_min = 0.0, phip_max = 2π, pt_min = 0.0, pt_max = 10.0, rtol = rtol)
    
    epsilon_ME = moments[1]
    P_ME = moments[2]
    piϕϕ_ME = moments[3]
    piηη_ME = moments[4]
    return SVector(epsilon_ME - epsilon_hydro, P_ME - P_hydro, piϕϕ_ME - piϕϕ_hydro, piηη_ME - piηη_hydro)
end

function reconstruct_moments_3D_LF(lambda_e, lambda_P, lambda_piϕϕ, lambda_piηη, degeneracies, masses; rtol = 1e-4)
    moments = LF_moments(lambda_e, lambda_P, lambda_piϕϕ, lambda_piηη, degeneracies, masses; etap_min = 0.0, etap_max = 10.0, phip_min = 0.0, phip_max = 2π, pt_min = 0.0, pt_max = 10.0, rtol = rtol)
    epsilon_ME = moments[1]
    P_ME = moments[2]
    piϕϕ_ME = moments[3]
    piηη_ME = moments[4]
    return epsilon_ME, P_ME, piϕϕ_ME, piηη_ME
end


function lagrange_multiplier_3D_LF_LRF(fo, hrg::Fluidum.HadronResonanceGas_LF, params; alpha_step)
    
    degeneracies = hrg.particle_list.Degeneracy
    masses = hrg.particle_list.Mass

    x = fo.x
    lb = Fluidum.leftbounds(x) 
    rb = Fluidum.rightbounds(x)
    alpha_array = range(lb[1], rb[1]; length = alpha_step)
    #alpha_array = range(0.32, rb[1]; length = alpha_step)
    
    epsilon_ME = zeros(length(alpha_array))
    P_ME = copy(epsilon_ME);
    piϕϕ_ME = copy(epsilon_ME);
    piηη_ME = copy(epsilon_ME);

    lm = [Lagr_Multiplier_3D_LF(0.0, 0.0, 0.0, 0.0) for _ in alpha_array]

    for (i, alpha) in enumerate(alpha_array)

        V = fo.fields(alpha)

        T = V[1]
        piϕϕ = V[3]/fmGeV^3
        piηη = V[4]/fmGeV^3
        piB = V[5]/fmGeV^3
        
        thermo = Fluidum.thermodynamic(T, params.eos)
        
        epsilon = (T*thermo.pressure_derivative[1] - thermo.pressure)/fmGeV^3
        P = thermo.pressure/fmGeV^3 + piB

        if i == 1
            guess = SVector(1/T, 0.0, 0.0, 0.0)
        else
            old_lm = lm[i-1]
            guess = SVector(old_lm.mult_e, old_lm.mult_P, old_lm.mult_piϕϕ, old_lm.mult_piηη)
        end
     
        _f(unknown, p) = MaximumEntropy.lagrangian_multipliers_system_3D(epsilon, P, piϕϕ, piηη, unknown, degeneracies, masses; rtol = 1e-4)

        problem = NonlinearProblem{false}(_f, guess)

        lambda_sol = solve(problem, TrustRegion(); reltol = 1e-7)

        lm[i] = Lagr_Multiplier_3D_LF(lambda_sol.u[1], lambda_sol.u[2], lambda_sol.u[3], lambda_sol.u[4])
        
        epsilon_ME[i], P_ME[i], piϕϕ_ME[i], piηη_ME[i] = MaximumEntropy.reconstruct_moments_3D_LF(lambda_sol.u[1], lambda_sol.u[2], lambda_sol.u[3], lambda_sol.u[4], degeneracies, masses; rtol = 1e-4)

        
        if abs(P_ME[i] - P) > 1e-5
            @warn "Starting diverging at alpha = $alpha"
        end
    end
    return lm, epsilon_ME, P_ME, piϕϕ_ME, piηη_ME
end
