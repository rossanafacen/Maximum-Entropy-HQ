struct Lagr_Multiplier_2D{A<:Real, B<:Real} <: Lagr_Multiplier
    mult_n::A
    mult_nu::B
end


"""define a charm distribution function that depends on two lagrange multipliers  
"""
function f_ME(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)

    p_transv = -ut*ur*mt*cosh(etap-eta)+(ut)^2*pt*cos(phip-phi)
    
    #arg = udotp+lm.mult_n/exp(T)-lm.mult_nu*pt*cos(phip-phi)/udotp/exp(T)
    arg = udotp/T+lm.mult_n-lm.mult_nu*p_transv/udotp
    #arg = lm.mult_n-lm.mult_nu*p_transv/udotp

    return exp(arg)
end


"""charm quark density, defined with negative sign due to the (-,+,+,+) metric"""
function charm_density_integrand(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(T,ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt
end

function charm_diff_current_integrand(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D; m)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    p_transv = -ut*ur*mt*cosh(etap-eta)+(ut)^2*pt*cos(phip-phi)
    return p_transv*f_ME(T,ur,eta,phi,etap,phip,pt,lm;m=m)/(2*pi)^3*pt
end


function charm_density(T,ur,lm::Lagr_Multiplier_2D;m,ccbar,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0.,pt_max=10.0,rtol=1e-6)
    eta = 0
    phi = 0
    fact = Fluidum.besseli(1, ccbar/2)./Fluidum.besseli(0, ccbar/2)
    hcubature( b->fact*2*fmGeV^3*charm_density_integrand(T,ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end

function charm_diff_current(T,ur,lm::Lagr_Multiplier_2D;m,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0.,pt_max=10.0,rtol=1e-6)
    eta = 0
    phi = 0
    hcubature(b->2*fmGeV^3*charm_diff_current_integrand(T,ur,eta,phi,b[1],b[2],b[3],lm;m=m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)[1]
end



function lagrangian_multipliers_system_2(T,ur,nur,n,unknown;m,ccbar)
    lm = Lagr_Multiplier_2D(unknown[1],unknown[2])
    eq_charm_1 = (charm_density(T,ur,lm; m = m,ccbar = ccbar)-n) / max(abs(n), 1e-30)
    eq_charm_2 = (charm_diff_current(T,ur,lm; m = m)-nur) / max(abs(nur), 1e-30)
    return SVector{2}(eq_charm_1,eq_charm_2)
end

function lambdan_lambdar(result, discretization::CartesianDiscretization, fluidproperty, t; guess_n = 0., guess_nu = 0., diff_n = 0.0,diff_nu = 0.0, m)  
    
    lambda_array = []
    #for i in (2)
    for i in div(lastindex(discretization.grid)-1,6):div(lastindex(discretization.grid)-1,3)
        T = result(t)[1,i]
        ur = result(t)[2,i]
        fug = result(t)[6,i]
        nur = result(t)[7,i]

        n = thermodynamic(T,fug,fluidproperty.eos.hadron_list).pressure
     
        _f(unknown,p) = lagrangian_multipliers_system_2(T,ur,nur,n,unknown;m)
        problem = NonlinearProblem{false}(_f,SVector{2}(guess_n,guess_nu))
        lambda_sol = solve(problem, LevenbergMarquardt();reltol = 1e-8, abstol = 1e-8)      
        
        lm = Lagr_Multiplier_2D(lambda_sol.u[1], lambda_sol.u[2])

        push!(lambda_array,lm)
           
        # guess_n = lambda_array[i-1].mult_n - diff_n
        # guess_nu = lambda_array[i-1].mult_nu - diff_nu

        # if (i>2)
        #     diff_n = abs(lambda_array[i-1].mult_n-lambda_array[i-2].mult_n)
        #     diff_nu = abs(lambda_array[i-1].mult_nu-lambda_array[i-2].mult_nu)
        # end

    end
    print("ciao Edu")
    return (lambda_array,t) 
end


function lambdan_lambdar_foo(fo, fluidproperty; alpha_step = 2, guess_nu = 0., guess_n = 0., diff_n = 0.0,diff_nu = 0.0, m)  
    print("starting function")
    @show fo.x    
    # lb=leftbounds(fo.x)
    # print("lb evaluated")
    
    # rb=rightbounds(fo.x)
    # @show rb, lb
    for alpha in range(lb...,rb...,alpha_step)
        @show alpha
        V = fo.fields(alpha)
        T = V[1] 
        ur = V[2]
        fug = V[6]
        nur = V[7]
        
        n = thermodynamic(T,fug,fluidproperty.eos.hadron_list).pressure
     
        _f(unknown,p) = lagrangian_multipliers_system_2(T,ur,nur,n,unknown;m)
        problem = NonlinearProblem{false}(_f,SVector{2}(guess_n,guess_nu))
        lambda_sol = solve(problem, NewtonRaphson();reltol = 1e-8, abstol = 1e-8)      
        
        
        lm = Lagr_Multiplier_2D(lambda_sol.u[1], lambda_sol.u[2])

        push!(lambda_array,lm)
           
        guess_n = lambda_array[i-1].mult_n - diff_n
        guess_nu = lambda_array[i-1].mult_nu - diff_nu

        if (i>2)
            diff_n = abs(lambda_array[i-1].mult_n-lambda_array[i-2].mult_n)
            diff_nu = abs(lambda_array[i-1].mult_nu-lambda_array[i-2].mult_nu)
        end

    end
    print("ciao Edu")
    return (lambda_array) 
end


function regularize_2d!(sol1,discretization)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    sol = [zeros(2) for _ in 1:length(discretization.grid)-2]
    j = 0
    for i in eachindex(x)
        sol[i][1]=sol1[i].u[1]
        if isnan(sol1[i].u[1])
            j = i 
            @show j
            break
        end
    end
    k = 0
    for i in eachindex(x)
        sol[i][2]=sol1[i].u[2]
        if isnan(sol1[i].u[2])
            k = i 
            @show k
            break
        end
    end

    for i in j:lastindex(discretization.grid)-2
        sol[i][1]=sol1[j-1].u[1]
    end
    for i in k:lastindex(discretization.grid)-2
        sol[i][2]=sol1[k-1].u[2]
    end
        sol1.=sol
end


function distr_function_2d(result, discretization::CartesianDiscretization, t,sol;m)  
    eta = 0
    phi = 0

    sol1 = []
    sol2 = []
    for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]
        
        lm.mult_n = sol[1][i]
        lm.mult_nu = sol[2][i]

        push!(sol1,f_ME(T,ur,eta,phi,etap,phip,pt,lm;m))   
        push!(sol2,f_ME(ur,T,eta,phi,etap,phip,pt,lm;m))   
    end
    return  (sol1,sol2)
end




