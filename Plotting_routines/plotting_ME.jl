"""Plot the r component of the charm diffusion current at different times."""
function plot_charm_diff_current_ME(field_result,lambda_sol,discretization::CartesianDiscretization)
    
    #times = [t1, t2, t3, t4]
    time = lambda_sol[2]
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    
    fig, ax = subplots()
    set_fields_axis(ax, 7)
         
    T = field_result(time)[1,:]
    ur = field_result(time)[2,:]
    fug = field_result(time)[6,:]

    nur = field_result(time)[7,:]
    
    ax.plot(x,[nur[i] for i in 2:lastindex(discretization.grid)-1], label = string(L"$\mathrm{hydro}, \mathrm{\tau =", time, L"\ \mathrm{fm}}/c$"))
    ax.plot(x,[charm_diff_current(T[i],ur[i],fug[i],lambda_sol[1][i-1]) for i in 2:lastindex(discretization.grid)-1], ls ="--",label = string(L"$\mathrm{ME}, \mathrm{\tau=", time, L"\ \mathrm{fm}}/c$"))

    
    ax.legend()
    display(fig)
end

"""Plot the density of charm at different times."""
function plot_charm_density_ME(field_result,lambda_sol,discretization::CartesianDiscretization,fluidproperty)
    time = lambda_sol[2]
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    
    fig, ax = subplots()
    set_fields_axis(ax, 8)
         
    T = field_result(time)[1,:]
    ur = field_result(time)[2,:]
    fug = field_result(time)[6,:]
    
    ax.plot(x,[[thermodynamic(T[i],fug[i],fluidproperty.eos.hadron_list).pressure] for i in 2:lastindex(discretization.grid)-1], label = string(L"$\mathrm{hydro}, \mathrm{\tau =", time, L"\ \mathrm{fm}}/c$"))
    ax.plot(x,[charm_density(T[i],ur[i],fug[i],lambda_sol[1][i-1]) for i in 2:lastindex(discretization.grid)-1], ls ="--",label = string(L"$\mathrm{ME}, \mathrm{\tau=", time, L"\ \mathrm{fm}}/c$"))
    
    ax.legend()
    display(fig)
end



function plot_lambda_ME(lambda_sol,discretization::CartesianDiscretization)
    time = lambda_sol[2]
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    
    fig, ax = subplots()
    set_fields_axis(ax, 7)
    ax.plot(x,[lambda_sol[1][i].mult_nu for i in 1:lastindex(discretization.grid)-2],label = string(L"$\mathrm{\tau=", time, L"\ \mathrm{fm}}/c$"))

    ax.legend()
    ax.set_ylabel(L"\lambda_r")
    display(fig)
end

function plot_f_ME(field_result,lambda_sol,discretization::CartesianDiscretization;etap = 0, phip = 0, pt = 0)
    time = lambda_sol[2]
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    fME=MaximumEntropy.distr_function_vs_r(field_result,lambda_sol,discretization; etap=etap, phip=phip, pt=pt)
    
    fig, ax = subplots()
    plot_params()

    ax.plot(x,[fME[i] for i in 1:lastindex(discretization.grid)-2],label = string(L"$\mathrm{\tau=", time, L"\ \mathrm{fm}}/c$"))       
    
    ax.text(0.6, 0.6, string(L"\eta_p =",etap, L", \phi_p =", phip, L", p_T =", pt, L"\ \mathrm{GeV}/c"))
    ax.legend()
    ax.set_xlabel("r [fm]")
    ax.set_ylabel(L"f_\mathrm{ME}")
    display(fig)
end
