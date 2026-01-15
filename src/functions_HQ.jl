abstract type Lagr_Multiplier end

struct Lagr_Multiplier_1D{A<:Real} <: Lagr_Multiplier
    mult_nu::A
end


function Base.zero(::Type{MaximumEntropy.Lagr_Multiplier_1D{Float64}})
    MaximumEntropy.Lagr_Multiplier_1D(0.0)
end

function Base.:+(a::Lagr_Multiplier_1D{Float64}, b::Real) 
    MaximumEntropy.Lagr_Multiplier_1D(a.mult_nu + float(b))
end 

"""charm distribution function that depends on one lagrange multiplier 
"""
function f_ME(T,ur,fug,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D;m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi) #negative defined
    p_transv = -ut*ur*mt*cosh(etap-eta)+(ut)^2*pt*cos(phip-phi)
    
    arg = udotp/T+fug+lm.mult_nu*p_transv/udotp  

    return exp(arg)
end


"""charm quark density, defined with negative sign due to the (-,+,+,+) metric"""
function charm_density_integrand(T,ur,fug,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(T,ur,fug,eta,phi,etap,phip,pt,lm;m)/(2*pi)^3*pt

end

"""radial diffusion current"""
function charm_diff_current_integrand(T,ur,fug,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    p_transv = -ut*ur*mt*cosh(etap-eta)+(ut)^2*pt*cos(phip-phi)
    return p_transv*f_ME(T,ur,fug,eta,phi,etap,phip,pt,lm;m)/(2*pi)^3*pt
end

"""total radial charm current"""
function charm_total_current_integrand(T,ur,fug,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return (pt*cos(phip-phi)-udotp*ur)*f_ME(T,ur,fug,eta,phi,etap,phip,pt,lm;m)/(2*pi)^3*pt
end

"""integrate the distribution function over the momentum space, to get the charm quark fields. 
Factor 2 comes from the rapidity integration, made only for positive rapidity.
In the integration of the distribution function, it is possible to impose eta_p = and phi_p = 0, since a shift for these coordinates does not change the integral result
"""
function charm_density(T,ur,fug,lm::Lagr_Multiplier_1D;m,ccbar,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0,pt_max=8.0,rtol=1e-5) 
    eta = 0 
    phi = 0
    fact = Fluidum.besseli(1, ccbar/2)./Fluidum.besseli(0, ccbar/2)
    hcubature( b->fact*2*fmGeV^3*charm_density_integrand(T,ur,fug,eta,phi,b[1],b[2],b[3],lm;m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

function charm_diff_current(T,ur,fug,lm::Lagr_Multiplier_1D;m,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0,pt_max=8.0,rtol=1e-5)
    eta = 0 
    phi = 0
    hcubature( b->2*fmGeV^3*charm_diff_current_integrand(T,ur,fug,eta,phi,b[1],b[2],b[3],lm;m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

function charm_total_current(T,ur,fug,lm::Lagr_Multiplier_1D;m,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0,pt_max=8.0,rtol=1e-5)
    eta = 0 
    phi = 0
    hcubature( b->2*fmGeV^3*charm_total_current_integrand(T,ur,fug,eta,phi,b[1],b[2],b[3],lm;m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end


function lagrangian_multipliers_system(T,ur,fug,nur,unknown;m)
    lm = Lagr_Multiplier_1D(unknown)
    eq_charm = charm_diff_current(T,ur,fug,lm;m)[1]-nur 
    return (eq_charm)
end


"""compute the lagrange multiplier over the discretized grid, equating the ME ν with the hydro ν. Store the time information too
"""
function lambdar(result, discretization::CartesianDiscretization, t; guess = 0.0, diff = 0.0, m = 1.5)  
    lambda_array = []    

    for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]
        fug = result(t)[6,i]
        nur = result(t)[7,i]

        _f(unknown,p) = lagrangian_multipliers_system(T,ur,fug,nur,unknown;m)
        problem = NonlinearProblem{false}(_f,guess)
        punctual_lambda_r = solve(problem,NewtonRaphson()) 
        lm = Lagr_Multiplier_1D(only(punctual_lambda_r.u))
        
        push!(lambda_array,lm)   
        
        guess = lambda_array[i-1].mult_nu + diff

        if (i>2)
            diff = abs(lambda_array[i-1].mult_nu-lambda_array[i-2].mult_nu)
        end
       
    end
    return (lambda_array,t)
end


"""compute the ME distribution function, using the Lagrange multiplier found in lambdar function 
"""
function distr_function_vs_r(result, sol, discretization::CartesianDiscretization; etap = 0, phip = 0, pt = 0, m = 1.5)  
    eta = 0
    phi = 0
    t = sol[2]
    distr_function = []
    for i in 2:lastindex(discretization.grid)-1
        
        T = result(t)[1,i]
        ur = result(t)[2,i]
        fug = result(t)[6,i]
        lm = sol[1][i-1]
        push!(distr_function,f_ME(T,ur,fug,eta,phi,etap,phip,pt,lm;m))   
    end
    return  distr_function
end

function regularize!(sol1,discretization)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]

    sol = similar(sol1)

    j = 0
    for i in eachindex(x)
        sol[i]=sol1[i].u
        if isnan(sol1[i].u)
            j = i #sol[i]=sol1[73].u
            @show j
            break
        end
    end

    for i in j:lastindex(discretization.grid)-2
        sol[i]=sol1[j-1].u
    end

    sol1.=sol
end

function not_regularize!(sol1,discretization)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    
    sol = similar(sol1)
    
    for i in eachindex(x)
        sol[i]=sol1[i].u
    end
    
    sol1.=sol
end



