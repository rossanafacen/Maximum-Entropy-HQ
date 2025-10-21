function charm_density_integrand(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff)
    eta = 0.0
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


function f_ME(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff)
    #@spack a, b, c = multipliers(...)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    eta = 0.0
    phi = 0.0
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    #pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2cos(phip-phi)^2+mt^2cosh(2etap-2eta))-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)
    #gpp::Float64 = -ur^2*(g1+g2)*(mt*cosh(etap-eta))^2-2*ut*ur*(g1+g2)*mt*pt*cosh(etap-eta)*cos(phip-phi)-ut^2*(g1+g2)*(pt*cos(phip-phi))^2
    arg = udotp/T+mult_n*mt*cosh(etap-eta)/udotp-mult_diff*pt*cos(phip-phi)/udotp
    #@show mult_diff*pt*cos(phip-phi)/udotp
    return exp(arg)
end

#ur = 0.1

function lagrangian_multipliers_system_2(unknown,ur,T,fug,m,n,nur)
    mult_n=unknown[1]
    mult_diff=unknown[2]
    eq_charm_1 = charm_density(ur,T,fug,m,mult_n,mult_diff)[1]-n
    eq_charm_2 = charm_current(ur,T,fug,m,mult_n,mult_diff)[1]-nur-n*ur
    #@show nur,n
    #@show charm_current(ur,T,fug,m,mult_n,mult_diff)[1], charm_density(ur,T,fug,m,mult_n,mult_diff)[1]
    return SVector{2}(eq_charm_1,eq_charm_2)
end

function distr_function_2(result, discretization::CartesianDiscretization, t,fluidpropery; guess_n = 5.0, guess_nu = -3, diff_n = 0.0,diff_nu = 0.0)  
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    #x = discretization.grid[10][1]
    sol1 = []
    sol2 = []
    #for i in (2)
    for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]
        ut = sqrt(ur^2+1)
        piϕϕ = result(t)[3,i]
        piηη = result(t)[4,i]
        piB = result(t)[5,i]
        fug = result(t)[6,i]
        nur = result(t)[7,i]

        n = federica(T,fug,fluidpropery.eos)[1]
        m = 1.5
        #guess_n = fug*5
        _f(unknown,p) = lagrangian_multipliers_system_2(unknown,ur,T,fug,m,n,nur)
        problem = NonlinearProblem{false}(_f,SVector{2}(guess_n,guess_nu))
        #mult_diff_temp = solve(problem,NewtonRaphson();maxiters=2) 
            #problem = NonlinearProblem{false}(_f,0)
            #mult_diff_temp= solve(problem,NewtonRaphson()) 
        #@show i, discretization.grid[i]
        
        #guess = mult_diff_temp[1]
        
        push!(sol1,solve(problem,NewtonRaphson()))   
        guess_n = sol1[i-1].u[1] + diff_n
        guess_nu = sol1[i-1].u[2] - diff_nu
        
        if (i>2)
        diff_n = abs(sol1[i-1].u[1]-sol1[i-2].u[1])
        diff_nu = abs(sol1[i-1].u[2]-sol1[i-2].u[2])
        
        end
        #mult_n_temp,mult_diff_temp= solve(problem,NewtonRaphson()) 
        #@show solve(problem,NewtonRaphson()).u
        
        #push!(sol2,mult_diff_temp)  
         
    end
    return (sol1)#(sol1,sol2)
end

