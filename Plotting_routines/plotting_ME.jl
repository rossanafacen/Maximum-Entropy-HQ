"""Plot the r component of the charm diffusion current at different times."""
function plot_charm_diff_current_ME(field_result,lambda_sol,discretization::CartesianDiscretization,tspan;t2 = 2., t3 = 6.)
    t1 = tspan[1]
    t4 = tspan[2]

    #times = [t1, t2, t3, t4]
    times = [t2]
    #x = [discretization.grid[i][1] for i in 2:50]
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    
    fig, ax = subplots()
    set_fields_axis(ax, 7)
    for time in times
         
        T = field_result(time)[1,:]
        ur = field_result(time)[2,:]
        mu = field_result(time)[6,:]

        nur = field_result(time)[7,:]
        
        ax.plot(x,[nur[i] for i in 2:lastindex(discretization.grid)-1], label = string(L"$\mathrm{\tau =", time, L"\ \mathrm{fm}}/c$"))
        ax.plot(x,[charm_diff_current(T[i],ur[i],mu[i],lambda_sol[i-1]) for i in 2:lastindex(discretization.grid)-1], ls ="--",label = string(L"$\mathrm{\tau=", time, L"\ \mathrm{fm}}/c$"))
         
    end
    ax.legend()
    display(fig)
end

function plot_lambda_ME(field_result,discretization::CartesianDiscretization,fluidproperty,tspan;t2 = 4, t3 = 8)
    t1 = tspan[1] 
    t4 = tspan[2]
    #times = [t1, t2, t3, t4]
    times = [t2]
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    fig, ax = subplots()
    set_fields_axis(ax, 7)
    for time in times
        lambda=lambdar(field_result,discretization,time)
        ax.plot(x,[lambda[i].mult_diff for i in 1:lastindex(discretization.grid)-2],label = string(L"$\mathrm{\tau=", time, L"\ \mathrm{fm}}/c$"))
        
        #ax.plot(x,lambda[eachindex(x)],label = string(L"$\mathrm{\tau=", time, L"\ \mathrm{fm}}/c$"))
    end
    ax.legend()
    ax.set_ylabel(L"\lambda_r")
    display(fig)
end

"""Plot the r and tau components of the charm total conserved current at different times."""
function plot_charm_total_current_ME(res,discretization::CartesianDiscretization,fluidproperty,tspan;t2 = 4., t3 = 6., m = 1.5)
    t1 = tspan[1]
    t4 = tspan[2]
    sol = distr_function(res,discretization,tstep,fluidproperty) 
    x = [discretization.grid[i][1] for i in 2:50]
    
    times = [t1, t2, t3, t4]
    for time in times
        T = res(time)[1,:]
        ur = res(time)[2,:]
        ut = sqrt.(1 .+ur.^2)
        mu = res(time)[6,:]
        
        nur = res(time)[7,:]
        nut = nur.*ur./ut
        
        #r component
        fig1, ax1 = subplots()
        ax1.plot(x,[thermodynamic(T,fug,eos.hadron_list).pressure*ur[i]+nur[i] for i in 2:50],label=string(time))
        ax1.plot(x,[charm_total_current(T[i],ur[i],mu[i],sol[i-1].u)[1] for i in 2:50], ls ="--",label=string(time)*"_mult")
        ax1.set_xlabel("r [fm]")
        ax1.set_ylabel(L"N^r [fm^{-3}]")
        ax1.legend()
        

        #tau component
        fig2, ax2 = subplots()
        ax2.plot(x,[thermodynamic(T,fug,eos.hadron_list).pressure*ut[i]+nut[i] for i in 2:50],label=string(time))
        ax2.plot(x,[charm_density(T[i],ur[i],mu[i],sol[i-1].u)[1] for i in 2:50], ls ="--",label=string(time)*"_mult")
        ax1.set_xlabel("r [fm]")
        ax1.set_ylabel(L"N^\tau [fm^{-3}]")
        ax2.legend()
    end
    display(fig1)
    display(fig2)
end





