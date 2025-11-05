"""charm distribution function that depends only on one lagrange multiplier 
"""
function f_ME(ur,T,fug,pt,m,eta,phi,etap,phip,mult_diff)
    #@spack a, b, c = multipliers(...)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    p0 = -mt*cosh(etap-eta) 

    arg = udotp/T+fug+mult_diff*pt*cos(phip-phi)/udotp
    return exp(arg)
end


"""charm quark density, defined with negative sign due to the (-,+,+,+) metric"""
function charm_density_integrand(ur,T,fug,pt,m,eta,phi,etap,phip,mult_diff)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(ur,T,fug,pt,m,eta,phi,etap,phip,mult_diff)/(2*pi)^3*pt 
    
end

"""radial diffusion current"""
function charm_diff_current_integrand(ur,T,fug,pt,m,eta,phi,etap,phip,mult_diff)
    return pt*cos(phip-phi)*f_ME(ur,T,fug,pt,m,eta,phi,etap,phip,mult_diff)/(2*pi)^3*pt   
end

"""total radial charge current"""
function charm_total_current_integrand(ur,T,fug,pt,m,eta,phi,etap,phip,mult_diff)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return (pt*cos(phip-phi)-udotp*ur)*f_ME(ur,T,fug,pt,m,eta,phi,etap,phip,mult_diff)/(2*pi)^3*pt 
end

"""integrate the distribution function over the momentum space, to get the charm quark fields. 
Factor 2 comes from the rapidity integration, made only for positive rapidity.
In the integration of the distribution funciton, it is possible to impose eta_p = and phi_p = 0, since a sihft for these coordinates does not change the integral result
"""
function charm_density(ur,T,m,fug,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-4) 
    eta = 0 
    phi = 0
   hcubature( b->2*fmGeV^3*charm_density_integrand(ur,T,fug,b[1],m,eta,phi,b[3],b[2],mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))
end

function charm_diff_current(ur,T,m,fug,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-4)
    eta = 0 
    phi = 0
    hcubature( b->2*fmGeV^3*charm_diff_current_integrand(ur,T,fug,b[1],m,eta,phi,b[3],b[2],mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))

end

function charm_total_current(ur,T,m,fug,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-4)
    eta = 0 
    phi = 0
    hcubature( b->2*fmGeV^3*charm_total_current_integrand(ur,T,fug,b[1],m,eta,phi,b[3],b[2],mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))
end


function lagrangian_multipliers_system(unknown,ur,T,fug,m,nur)
    mult_diff=unknown
    eq_charm = charm_diff_current(ur,T,m,fug,mult_diff)[1]-nur 
    return (eq_charm)
end


function distr_function(result, discretization::CartesianDiscretization, t; guess = 0.0, diff = 0.0, m = 1.5)  
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    sol1 = []    

    for i in 2:50
    #for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]
        fug = result(t)[6,i]
        nur = result(t)[7,i]

           
        _f(unknown,p) = lagrangian_multipliers_system(unknown,ur,T,fug,m,nur)
        problem = NonlinearProblem{false}(_f,guess)
        mult_diff_temp= solve(problem,NewtonRaphson()) 
            
        push!(sol1,mult_diff_temp)   
        guess = sol1[i-1].u + diff
        if (i>2)
            diff = abs(sol1[i-1].u-sol1[i-2].u)
        end
         
    end
    return (sol1)
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


function distr_function_wrap(result, discretization::CartesianDiscretization, t,sol)  
    eta = 0
    phi = 0
    #x = discretization.grid[10][1]
    sol1 = []
    for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]
        fug = result(t)[6,i]
        mult_diff = sol[i]
        push!(sol1,f_ME(ur,T,fug,pt,m,eta,phi,etap,phip,mult_diff))   
    end
    return  sol1
end
