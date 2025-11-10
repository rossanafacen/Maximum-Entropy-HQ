"""
transverse projector 
projects transverse to fluid velocity u^μ, indices up up 
"""
@inline function projector(u,g)
    @tullio Δ[i,j] := inv(g)[i,j]+u[i]*u[j]
end


@inline function deldel(u,g)
    @tullio dd[i,j,k,l] := projector(u,g)[i,j]*projector(u,g)[k,l]
end

"""
transverse projector rank 4
projects transverse to fluid velocity u^μ, indices up up up up 
"""
@inline function del4(u,g)
    @tullio d4[i,j,k,l] := 1/2*projector(u,g)[i,k]*projector(u,g)[j,l]+1/2*projector(u,g)[i,l]*projector(u,g)[j,k]-
    1/3*projector(u,g)[i,j]*projector(u,g)[k,l] 
end


"""
Lagrange Multipliers for energy, bulk and shear viscosity 
"""
struct multipliers{T,V}
    λE::T
    mult_P::T
    λshear::V
end

#Most general case in which no symmetries are considered

"""
Maximum Entropy distribution function with energy and bulk constraints 
"""
function f_ME(u::AbstractVector,p,g,mult_E,mult_P)
    #@spack a, b, c = multipliers(...)
    udotp = dot(u,g,p)
    @tullio transv_p_up[i] := projector(u,g)[i,j]*g[j,k]*p[k] #p^⟨μ⟩
    arg = -mult_E*udotp-mult_P/udotp*dot(transv_p_up,g,transv_p_up)
    return exp(arg)
end

"""
Maximum Entropy distribution function with energy, bulk and shear constraints
"""
function f_ME(u::AbstractVector,p,g,mult_E,mult_P,mult_shear)
    #@spack a, b, c = multipliers(...)
    udotp = dot(u,g,p)
    @tullio transv_p_up[i] := projector(u,g)[i,j]*g[j,k]*p[k] 
    @tullio transv_pp_upup[i,j] := del4(u,g)[i,j,a,b]*g[a,m]*g[b,n]*p[m]*p[n] 
    @tullio mult_transv := mult_shear[i,j]*transv_pp_upup[i,j]
    arg = -mult_E*udotp-mult_P/udotp*dot(transv_p_up,g,transv_p_up)-mult_transv/udotp
    return exp(arg)
end

"""
Integrand of the energy density in dP, as defined in kinetic theory
"""
function energy_integrand(u::AbstractVector,p,g,mult_E,mult_P)
    udotp = dot(u,g,p)
    return udotp^2*f_ME(u,p,g,mult_E,mult_P)/(2*pi)^3/abs(p[1])
end

function energy_integrand(u::AbstractVector,p,g,mult_E,mult_P,mult_shear)
    udotp = dot(u,g,p)
    return udotp^2*f_ME(u,p,g,mult_E,mult_P,mult_shear)/(2*pi)^3/abs(p[1])
end

function pressure_integrand(u::AbstractVector,p,g,mult_E,mult_P)
    @tullio transv_p_up[i] := projector(u,g)[i,j]*g[j,k]*p[k]
    return -1/3*dot(transv_p_up,transv_p_up,g)*f_ME(u,p,g,mult_E,mult_P)/(2*pi)^3/abs(p[1])
end

function pressure_integrand(u::AbstractVector,p,g,mult_E,mult_P,mult_shear)
    @tullio transv_p_up[i] := projector(u,g)[i,j]*g[j,k]*p[k]
    return -1/3*dot(transv_p_up,transv_p_up,g)*f_ME(u,p,g,mult_E,mult_P,mult_shear)/(2*pi)^3/abs(p[1])
end

function shear_integrand(u::AbstractVector,p,g,mult_E,mult_P)
    @tullio transv_pp_upup[i,j] := del4(u,g)[i,j,a,b]*g[a,m]*g[b,n]*p[m]*p[n]
    return @tullio pi_shear[i,j] := transv_pp_upup[i,j]*f_ME(u,p,g,mult_E,mult_P)/(2*pi)^3*abs(p[1])
end

function shear_integrand(u::AbstractVector,p,g,mult_E,mult_P,mult_shear)
    @tullio transv_pp_upup[i,j] := del4(u,g)[i,j,a,b]*g[a,m]*g[b,n]*p[m]*p[n]
    return @tullio pi_shear[i,j] := transv_pp_upup[i,j]*f_ME(u,p,g,mult_E,mult_P,mult_shear)/(2*pi)^3*abs(p[1])
end


function energy_density(u::AbstractVector,g,mult_E,mult_P;pt_min=0.2,pt_max=10.,rtol =1000*sqrt(eps()))
    hcubature(b->4*π*b^2*energy_integrand(u,b,g,mult_E,mult_P),pt_min,pt_max;rtol=rtol)
end

function energy_density(u::AbstractVector,g,mult_E,mult_P,mult_shear;pt_min=0.2,pt_max=10.,rtol =1000*sqrt(eps()))
    hcubature(b->4*π*b^2*energy_integrand(u,b,g,mult_E,mult_P,mult_shear),pt_min,pt_max;rtol=rtol)
end

function pressure_density(u::AbstractVector,g,mult_E,mult_P;pt_min=0.2,pt_max=10.,rtol =1000*sqrt(eps()))
    hcubature(b->4*π*b^2*pressure_integrand(u,b,g,mult_E,mult_P),pt_min,pt_max;rtol=rtol)
end

function pressure_density(u::AbstractVector,g,mult_E,mult_P,mult_shear;pt_min=0.2,pt_max=10.,rtol =1000*sqrt(eps()))
    hcubature(b->4*π*b^2*pressure_integrand(u,b,g,mult_E,mult_P,mult_shear),pt_min,pt_max;rtol=rtol)
end


function lagrangian_multipliers_system(u::AbstractVector,g,mult_E,mult_P,E,P_tot)
    eq_energy = energy_density(u,g,mult_E,mult_P)-E
    eq_pressure = pressure_density(u,g,mult_E,mult_P)-P_tot
    return SVector{2}(eq_energy,eq_pressure)
end

function lagrangian_multipliers_system(u::AbstractVector,g,mult_E,mult_P,mult_shear,E,P_tot,piϕϕ,piηη)
    eq_energy = energy_density(u,g,mult_E,mult_P,mult_shear)-E
    eq_pressure = pressure_density(u,g,mult_E,mult_P,mult_shear)-P_tot
    eq_shear_ϕϕ = shear_phiphi(u,g,mult_E,mult_P,mult_shear)-piϕϕ
    eq_shear_ηη = shear_etaeta(u,g,mult_E,mult_P,mult_shear)-piηη
    return SVector{4}(eq_energy,eq_pressure,eq_shear_ϕϕ,eq_shear_ηη)
end

function solve_lambda(result, discretization::CartesianDiscretization, t,fluidproperty)  
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]

    mult_E = []
    mult_P = []
    
    for i in x
        T = result(t)[1,i]
        ur = result(t)[2,i]
        ut = sqrt(ur^2+1)
        u = SVector{4}(ut,ur,0,0)
        piB = result(t)[5,i]

        g = SMatrix{4,4}([-1 0 0 0; 0 1 0 0; 0 0 1/i^2 0; 0 0 0 1/t^2]) #indices up

        P_tot = pressure(T,fluidproperty.eos)+piB
        E = T*pressure_derivative(T,Val(1),fluidproperty.eos)-pressure(T,fluidproperty.eos)

        _f = lagrangian_multipliers_system(u,g,mult_E,mult_P,E,P_tot)
        problem = NonlinearProblem{false}(_f, SVector{2,N}(mult_E,mult_P))
        mult_E_temp,mult_P_temp = solve(problem,NewtonRaphson(), tol = 1e-6) 
        push!(mult_E,mult_E_temp)
        push!(mult_P,mult_P_temp)
    end
    return (mult_E,mult_P)
end


#From now on, we will consider our case of 1+1 system

"""
Maximum Entropy distribution function expressed considering a symmetric 1+1 system, with energy and bulk constraints 
"""
function f_ME(ur::Real,pt,m,etap,phip,mult_E,mult_P)
    #@spack a, b, c = multipliers(...)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    eta = 0.0
    phi = 0.0
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2*cos(phip-phi)^2+mt^2*cosh(etap-eta)^2)-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta) #p^⟨μ⟩ p_⟨μ⟩ 
    arg = -mult_E*udotp-mult_P/udotp*pp
    return exp(arg)
end

"""
Maximum Entropy distribution function expressed considering a symmetric 1+1 system, with energy, bulk and shear constraints 
"""
function f_ME(ur::Real,pt,m,etap,phip,mult_E,mult_P,g1,g2)
    #@spack a, b, c = multipliers(...)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    eta = 0.0
    phi = 0.0
    udotp = -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2*cos(phip-phi)^2+mt^2*cosh(etap-eta)^2)-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)  
    
    gpp = -ur^2*(g1+g2)*(mt*cosh(etap-eta))^2-2*ut*ur*(g1+g2)*mt*pt*cosh(etap-eta)*cos(phip-phi)-ut^2*(g1+g2)*(pt*cos(phip-phi))^2 #p^⟨μ⟩ p^⟨ν⟩ γ_⟨μν⟩
    arg = -mult_E*udotp-mult_P/udotp*pp-1.0/udotp*gpp
    return exp(arg)
end

"""
Integrand of the energy density in p_T*dp_T, as defined in kinetic theory
"""
function energy_integrand(ur::Real,pt,m,etap,phip,mult_E,mult_P)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return udotp^2*f_ME(ur,pt,m,etap,phip,mult_E,mult_P)/(2*pi)^3*pt
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


function pressure_integrand(ur::Real,pt,m,etap,phip,mult_E,mult_P)
    eta = 0
    phi = 0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2*cos(phip-phi)^2+mt^2*cosh(etap-eta)^2)-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)  
    return -1/3*pp*f_ME(ur,pt,m,etap,phip,mult_E,mult_P)/(2*pi)^2*pt
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

function energy_density(ur::Real,m,mult_E,mult_P;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->energy_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end

function energy_density(ur::Real,m,mult_E,mult_P,g1,g2;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->energy_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end

function pressure_density(ur::Real,m,mult_E,mult_P;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->pressure_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
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



function lagrangian_multipliers_system(unknown,ur::Real,m,E,P_tot)
    mult_E=unknown[1]
    mult_P=unknown[2]
    eq_energy = energy_density(ur,m,mult_E,mult_P)[1]-E
    eq_pressure = pressure_density(ur,m,mult_E,mult_P)[1]-P_tot
    return SVector{2}(eq_energy,eq_pressure)
end

function lagrangian_multipliers_system(ur::Real,m,mult_E,mult_P,g1,g2,E,P_tot,piϕϕ,piηη)
    eq_energy = energy_density(ur,m,mult_E,mult_P,g1,g2)[1]-E
    eq_pressure = pressure_density(ur,m,mult_E,mult_P,g1,g2)[1]-P_tot
    eq_shear_ϕϕ = shear(ur,m,mult_E,mult_P,g1,g2)[1][1]-piϕϕ
    eq_shear_ηη = shear(ur,m,mult_E,mult_P,g1,g2)[2][1]-piηη
    return SVector{4}(eq_energy,eq_pressure,eq_shear_ϕϕ,eq_shear_ηη)
end


"""
Consider the case of energy and bulk viscosity only
"""
function distr_function_symmetric(result, discretization::CartesianDiscretization, t,fluidproperty)  
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    x = discretization.grid[10][1]
    sol1 = []
    sol2 = []
    for i in (1)
    #for i in eachindex(x)
        T = result(t)[1,i]
        ur = result(t)[2,i]
        ut = sqrt(ur^2+1)
        piB = result(t)[5,i]

        P_tot = pressure(T,fluidproperty.eos)+piB
        E = T*pressure_derivative(T,Val(1),fluidproperty.eos)-pressure(T,fluidproperty.eos)
        m = 0.14
        _f(unknown,p) = lagrangian_multipliers_system(unknown,ur,m,E,P_tot)
        problem = NonlinearProblem{false}(_f,SVector{2}(-1/T,0.))
        #@show problem
        mult_E_temp,mult_P_temp = solve(problem,NewtonRaphson(); reltol = 1e-6,maxiters = 100) 
        push!(sol1,mult_E_temp) 
        push!(sol2,mult_P_temp)    
    end
    
    return (sol1,sol2)
end




