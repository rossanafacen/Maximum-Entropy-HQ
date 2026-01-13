struct Lagr_Multiplier_3D{A<:Real, B<:Real, C<:Real} <: Lagr_Multiplier
    mult_e::A
    mult_n::B
    mult_nu::C
end

"""define a charm distribution function that depends on two lagrange multipliers  
"""
function f_ME(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_3D; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)

    arg = lm.mult_e*udotp+lm.mult_n-lm.mult_nu*pt*cos(phip-phi)/udotp
    
    return exp(arg)
end

function charm_energy_density_integrand(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_3D; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return udotp^2*f_ME(ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt
end

"""charm quark density, defined with negative sign due to the (-,+,+,+) metric"""
function charm_density_integrand(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_3D; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt
end

function charm_current_integrand(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_3D; m)
    return pt*cos(phip-phi)*f_ME(ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt
end


function charm_energy_density(ur,lm::Lagr_Multiplier_3D;m,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0.,pt_max=8.0,rtol=10E-6)
    eta = 0
    phi = 0 
    hcubature(b->2*fmGeV^3*charm_energy_density_integrand(ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

function charm_density(ur,lm::Lagr_Multiplier_3D;m,ccbar,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0.,pt_max=10.0,rtol=10E-6)
    eta = 0
    phi = 0
    fact = Fluidum.besseli(1, ccbar/2)./Fluidum.besseli(0, ccbar/2)
    hcubature( b->fact*2*fmGeV^3*charm_density_integrand(ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

function charm_diff_current(ur,lm::Lagr_Multiplier_3D;m,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0.,pt_max=10.0,rtol=10E-6)
    eta = 0
    phi = 0
    hcubature(b->2*fmGeV^3*charm_current_integrand(ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end



function lagrangian_multipliers_system_3(ur,epsilon,n,nur,unknown;m,ccbar)
    lm = Lagr_Multiplier_3D(unknown[1],unknown[2],unknown[3])
    eq_charm_1 = charm_energy_density(ur,lm; m = m)-epsilon
    eq_charm_2 = charm_density(ur,lm; m = m,ccbar = ccbar)-n
    eq_charm_3 = charm_diff_current(ur,lm; m = m)-nur
    return SVector{3}(eq_charm_1,eq_charm_2,eq_charm_3)
end

function lambdae_lambdan_lambdar(result, discretization::CartesianDiscretization, fluidproperty, t; guess_n = 0., guess_nu = 0., diff_n = 0.0,diff_nu = 0.0, m)  
    
    lambda_array = []
    #for i in (2)
    for i in div(lastindex(discretization.grid)-1,6):div(lastindex(discretization.grid)-1,3)
        T = result(t)[1,i]
        ur = result(t)[2,i]
        fug = result(t)[6,i]
        nur = result(t)[7,i]

        n = thermodynamic(T,fug,fluidproperty.eos.hadron_list).pressure
        epsilon = T*pressure_derivative(T,Val(1),fluidproperty.eos)-pressure(T,fluidproperty.eos)
        
        _f(unknown,p) = lagrangian_multipliers_system_3(ur,epsilon,n,nur,unknown;m,ccbar)
        
        problem = NonlinearProblem{false}(_f,SVector{3}(guess_e,guess_n,guess_nu))
        lambda_sol = solve(problem, LevenbergMarquardt();reltol = 1e-8, abstol = 1e-8)      
        
        lm = Lagr_Multiplier_3D(lambda_sol.u[1], lambda_sol.u[2], lambda_sol.u[3])

        push!(lambda_array,lm)
           
    end
    return (lambda_array,t) 
end