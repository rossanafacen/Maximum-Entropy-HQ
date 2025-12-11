struct Lagr_Multiplier_LF{A<:Real} <: Lagr_Multiplier
    mult_e::A
end

"""
Maximum Entropy distribution function expressed considering a symmetric 1+1 system, with energy and bulk constraints. The udopt/T is positive because we are employing the mostly plus metric 
"""
function f_ME(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_LF;m = 1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    
    arg = lm.mult_e*udotp
    return exp(arg)
end

"""
Integrand of the energy density in p_T*dp_T, as defined in kinetic theory
"""
function energy_integrand(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_LF; m = 1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return udotp^2*f_ME(ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt
end


"""
Energy density integrate p_T*dp_T, as defined in kinetic theory
"""
function energy_density(ur,lm::Lagr_Multiplier_LF;m=1.5,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0.,pt_max=8.0,rtol=10E-6)
    eta = 0
    phi = 0 
    hcubature(b->2*fmGeV^3*energy_integrand(ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end


function lagrangian_multipliers_system_LF_1D(ur,epsilon,unknown; m = 1.5)
    lm = Lagr_Multiplier_LF(unknown)
    eq_energy = energy_density(ur,lm;m=m)-epsilon
    return (eq_energy)
end


"""
Consider the case of energy and bulk viscosity only
"""
function lambda_res_1d_LF(result, discretization::CartesianDiscretization,fluidproperty,t; guess_n = 0., guess_nu = 0., diff_n = 0.0,diff_nu = 0.0, m = 1.5)  
    lambda_array = []
    
    for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]

        epsilon = T*pressure_derivative(T,Val(1),fluidproperty.eos)-pressure(T,fluidproperty.eos)
        
        _f(unknown,p) = lagrangian_multipliers_system_LF_1D(ur,epsilon,unknown)
        problem = NonlinearProblem{false}(_f,-1/T)
        
        lambda_sol = solve(problem,NewtonRaphson(); reltol = 1e-6) 
        lm = Lagr_Multiplier_LF(lambda_sol.u)

        push!(lambda_array,lm)
           
    end
    
    return (lambda_array,t)
end


