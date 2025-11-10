struct Lagr_Multiplier_2D{A<:Float64, B<:Float64} <: Lagr_Multiplier
    mult_n::A
    mult_diff::B
end

"""define a charm distribution function that depends on two lagrange multipliers  
"""
function f_ME(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D; m = 1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)

    arg = udotp/T+lm.mult_n-lm.mult_diff*pt*cos(phip-phi)/udotp
    return exp(arg)
end

"""charm quark density, defined with negative sign due to the (-,+,+,+) metric"""
function charm_density_integrand(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D; m = 1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(T,ur,eta,phi,etap,phip,pt,lm;m)/(2*pi)^3*pt
end

function charm_current_integrand(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D; m = 1.5)
    return pt*cos(phip-phi)*f_ME(T,ur,eta,phi,etap,phip,pt,lm;m)/(2*pi)^3*pt
end


function charm_density(T,ur,lm::Lagr_Multiplier_2D;m=1.5,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0.,pt_max=8.0,rtol=10E-2)
    eta = 0
    phi = 0
    hcubature( b->2*fmGeV^3*charm_density_integrand(T,ur,eta,phi,b[1],b[2],b[3],lm;m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)
end

function charm_current(T,ur,lm::Lagr_Multiplier_2D;m=1.5,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0.,pt_max=8.0,rtol=10E-2)
    eta = 0
    phi = 0
    hcubature( b->2*fmGeV^3*charm_current_integrand(T,ur,eta,phi,b[1],b[2],b[3],lm;m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max);rtol=rtol)
end




function lagrangian_multipliers_system_2(T,ur,nur,n,unknown;m=1.5)
    lm.mult_n=unknown[1]
    lm.mult_diff=unknown[2]
    eq_charm_1 = charm_density(T,ur,lm; m)[1]-n
    eq_charm_2 = charm_current(T,ur,lm; m)[1]-nur

    return SVector{2}(eq_charm_1,eq_charm_2)
end

function distr_function_2d(result, discretization::CartesianDiscretization, t; guess_n = 5.0, guess_nu = -3, diff_n = 0.0,diff_nu = 0.0, m = 1.5)  
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
     
        _f(unknown,p) = lagrangian_multipliers_system_2(T,ur,nur,n,unknown;m)
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
    return (sol1,sol2) 
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


function distr_function_wrap_2d(result, discretization::CartesianDiscretization, t,sol;m = 1.5)  
    eta = 0
    phi = 0

    sol1 = []
    sol2 = []
    for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]
        
        lm.mult_n = sol[1][i]
        lm.mult_diff = sol[2][i]

        push!(sol1,f_ME(T,ur,eta,phi,etap,phip,pt,lm;m))   
        push!(sol2,f_ME(ur,T,eta,phi,etap,phip,pt,lm;m))   
    end
    return  (sol1,sol2)
end




