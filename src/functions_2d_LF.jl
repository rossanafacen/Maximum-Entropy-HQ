struct Lagr_Multiplier_2D_LF{A<:Real, B<:Real} <: Lagr_Multiplier
    mult_e::A
    mult_P::B
end

"""
Maximum Entropy distribution function expressed considering a symmetric 1+1 system, with energy and bulk constraints 
"""
function f_ME(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D_LF;m = 1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pp = pt^2+mt^2/2*(cosh(2*etap-2*eta)-1)+ur^2*(pt^2*cos(phip-phi)^2+mt^2*cosh(etap-eta)^2)-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta) #p^⟨μ⟩ p_⟨μ⟩ 
    
    arg = lm.mult_e*udotp-lm.mult_P/udotp*pp
    return exp(arg)
end

"""
Integrand of the energy density in p_T*dp_T, as defined in kinetic theory
"""
function energy_integrand(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D_LF; m = 1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return udotp^2*f_ME(ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt
end

function pressure_integrand(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D_LF; m = 1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    pp = pt^2+mt^2/2*(cosh(2*etap-2*eta)-1)+ur^2*(pt^2*cos(phip-phi)^2+mt^2*cosh(etap-eta)^2)-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)  
    return 1/3*pp*f_ME(ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt #+1/3 because of our mostly-plus metric 
end

"""
Energy density integrate p_T*dp_T, as defined in kinetic theory
"""
function energy_density(ur,lm::Lagr_Multiplier_2D_LF;m=1.5,etap_min=0,etap_max=1,phip_min=0,phip_max=2pi,pt_min=0.2,pt_max=8.0,rtol=10E-5)
    eta = 0
    phi = 0 
    hcubature(b->2*fmGeV^3*energy_integrand(ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

function pressure_total(ur,lm::Lagr_Multiplier_2D_LF;m=1.5,etap_min=0,etap_max=1,phip_min=0,phip_max=2pi,pt_min=0.2,pt_max=8.0,rtol=10E-5)
    eta = 0
    phi = 0
    hcubature( b->2*fmGeV^3*pressure_integrand(ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

function lagrangian_multipliers_system_LF(ur,epsilon,P_tot,unknown; m = 1.5)
    lm = Lagr_Multiplier_2D_LF(unknown[1],unknown[2])
    eq_energy = energy_density(ur,lm;m=m)-epsilon
    eq_pressure = pressure_total(ur,lm;m=m)-P_tot
    return SVector{2}(eq_energy,eq_pressure)
end


"""
Consider the case of energy and bulk viscosity only
"""
function lambda_res_LF(result, discretization::CartesianDiscretization,fluidproperty,t; guess_n = 0., guess_nu = 0., diff_n = 0.0,diff_nu = 0.0, m = 1.5)  
    lambda_array = []
    
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    
    for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]
        piB = result(t)[5,i]

        epsilon = T*pressure_derivative(T,Val(1),fluidproperty.eos)-pressure(T,fluidproperty.eos)
        P_tot = pressure(T,fluidproperty.eos)+piB
        
        _f(unknown,p) = lagrangian_multipliers_system_LF(ur,epsilon,P_tot,unknown)
        problem = NonlinearProblem{false}(_f,SVector{2}(-1/T,0.))
        
        
        lambda_sol = solve(problem,NewtonRaphson(); reltol = 1e-6) 
        lm = Lagr_Multiplier_2D_LF(lambda_sol.u[1], lambda_sol.u[2])

        push!(lambda_array,lm)
           
    end
    
    return (lambda_array,t)
end


####fix also the shear viscosity
"""
Maximum Entropy distribution function expressed considering a symmetric 1+1 system, with energy, bulk and shear constraints 
"""
function f_ME(ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D_LF,g1,g2; m = 1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp = -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    
    pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2*cos(phip-phi)^2+mt^2*cosh(etap-eta)^2)-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)  
    
    gpp = -ur^2*(g1+g2)*(mt*cosh(etap-eta))^2-2*ut*ur*(g1+g2)*mt*pt*cosh(etap-eta)*cos(phip-phi)-ut^2*(g1+g2)*(pt*cos(phip-phi))^2 #p^⟨μ⟩ p^⟨ν⟩ γ_⟨μν⟩
    arg = -lm.mult_e*udotp-lm.mult_P/udotp*pp-1.0/udotp*gpp
    return exp(arg)
end



function energy_integrand(ur::Real,pt,m,etap,phip,mult_E,mult_P,g1,g2)
    eta = 0
    phi = 0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return udotp^2*f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)/(2*pi)^2*pt
end



"""
Integrand of ε*u^⟨μ⟩, as defined in kinetic theory
"""
function velocity_integrand(ur::Real,pt,m,etap,phip,mult_E,mult_P)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pr=pt*cos(phip-phi)
    return udotp^2*f_ME(ur,pt,m,etap,phip,mult_E,mult_P)/(2*pi)^2*pt*pr
end

function velocity_integrand(ur::Real,pt,m,etap,phip,mult_E,mult_P,g1,g2)
    eta = 0
    phi = 0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pr=pt*cos(phip-phi)
    return udotp^2*f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)/(2*pi)^2*pt*pr
end




function pressure_integrand(ur::Real,pt,m,etap,phip,mult_E,mult_P,g1,g2)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2*cos(phip-phi)^2+mt^2*cosh(2etap-2eta))-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)
    return -1/3*pp*f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)/(2*pi)^2*pt
end

function shear_integrand(ur::Real,pt,m,etap,phip,mult_E,mult_P,g1,g2)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    pp_phiphi = 1/6*(mt^2 + pt^2 - 3*pt^2 *cos(2phi - 2phip) - mt^2 *cosh(2eta - 2etap) -     
    2*(pt^2 *cos(phi - phip)^2 + mt^2*cosh(eta - etap)^2)*ur^2 +     
    4*mt* pt *cos(phi - phip) *cosh(eta - etap)* ur*ut)
    pp_etaeta = -1/3*(mt^2 + pt^2 - mt^2 *cosh(2eta - 2etap) +
     (pt^2 *cos(phi - phip)^2 + mt^2 *cosh(eta - etap)^2)*ur^2 - 
    2*mt* pt *cos(phi - phip) *cosh(eta - etap)* ur*ut)
    return (piphiphi=f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)*pp_phiphi/(2*pi)^2*pt,pietaeta=f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)*pp_etaeta/(2*pi)^2*pt)
end


function energy_density(ur::Real,m,mult_E,mult_P,g1,g2;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->energy_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end



function pressure_density(ur::Real,m,mult_E,mult_P,g1,g2;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->pressure_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end

function velocity(ur::Real,m,mult_E,mult_P;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    e::Float64 = energy_density(ur,m,mult_E,mult_P)[1]
    hcubature( b->velocity_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)./e
end

function velocity(ur::Real,m,mult_E,mult_P,g1,g2;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    e::Float64=energy_density(ur,m,mult_E,mult_P,g1,g2)[1]
    hcubature( b->velocity_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)./e
end


function shear(ur::Real,m,mult_E,mult_P,g1,g2;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    (hcubature( b->shear_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2)[:piphiphi],(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol),
    hcubature( b->shear_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2)[:pietaeta],(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol))
end




function lagrangian_multipliers_system(ur::Real,m,mult_E,mult_P,g1,g2,E,P_tot,piϕϕ,piηη)
    eq_energy = energy_density(ur,m,mult_E,mult_P,g1,g2)[1]-E
    eq_pressure = pressure_density(ur,m,mult_E,mult_P,g1,g2)[1]-P_tot
    eq_shear_ϕϕ = shear(ur,m,mult_E,mult_P,g1,g2)[1][1]-piϕϕ
    eq_shear_ηη = shear(ur,m,mult_E,mult_P,g1,g2)[2][1]-piηη
    return SVector{4}(eq_energy,eq_pressure,eq_shear_ϕϕ,eq_shear_ηη)
end





