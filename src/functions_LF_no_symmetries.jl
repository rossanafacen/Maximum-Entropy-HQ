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


