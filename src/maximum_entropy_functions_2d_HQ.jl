"""define a charm distribution function that depends on two lagrange multipliers  
"""
function f_ME(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    eta = 0.0
    phi = 0.0
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    
    arg = udotp/T+mult_n-mult_diff*pt*cos(phip-phi)/udotp 
    
    return exp(arg)
end


function charm_density_integrand(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff)
    eta = 0.0 #fix
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff)/(2*pi)^3*pt
end

function charm_current_integrand(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return pt*cos(phip-phi)*f_ME(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff)/(2*pi)^3*pt
end


function charm_density(ur,T,m,fug,mult_n,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-2)
    hcubature( b->2*fmGeV^3*charm_density_integrand(ur,T,fug,b[1],m,b[3],b[2],mult_n,mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end

function charm_current(ur,T,fug,m,mult_n,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-2)
    hcubature( b->2*fmGeV^3*charm_current_integrand(ur,T,fug,b[1],m,b[3],b[2],mult_n,mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end





function lagrangian_multipliers_system_2(unknown,ur,T,fug,m,n,nur)
    mult_n=unknown[1]
    mult_diff=unknown[2]
    eq_charm_1 = charm_density(ur,T,fug,m,mult_n,mult_diff)[1]-n
    eq_charm_2 = charm_current(ur,T,fug,m,mult_n,mult_diff)[1]-nur
    
    return SVector{2}(eq_charm_1,eq_charm_2)
end

function distr_function_2(result, discretization::CartesianDiscretization, t; guess_n = 5.0, guess_nu = -3, diff_n = 0.0,diff_nu = 0.0, m = 1.5)  
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    sol1 = []
    sol2 = []
    #for i in (2)
    for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]
        fug = result(t)[6,i]
        nur = result(t)[7,i]

       n = thermodynamic(T,fug,eos.hadron_list).pressure
     
        _f(unknown,p) = lagrangian_multipliers_system_2(unknown,ur,T,fug,m,n,nur)
        problem = NonlinearProblem{false}(_f,SVector{2}(guess_n,guess_nu))
       
        mult_n_temp,mult_diff_temp= solve(problem,NewtonRaphson())      
        
        push!(sol1,mult_n_temp)
        push!(sol2,mult_diff_temp)  
           
        guess_n = sol1[i-1].u[1] + diff_n
        guess_nu = sol1[i-1].u[2] - diff_nu
        
        if (i>2)
            diff_n = abs(sol1[i-1].u[1]-sol1[i-2].u[1])
            diff_nu = abs(sol1[i-1].u[2]-sol1[i-2].u[2])  
        end        
         
    end
    return (sol1,sol2) #(sol1)
end

#mult_diff_temp = solve(problem,NewtonRaphson();maxiters=2) 
#problem = NonlinearProblem{false}(_f,0)
#mult_diff_temp= solve(problem,NewtonRaphson()) 
#guess = mult_diff_temp[1]