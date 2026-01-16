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

"""plot the ME spectra and compare it with the multiple species"""
function spectra_plot_ME(fo::FreezeOutResult{A,B}, lm_funct, particle, detector, fluidproperty; data = false, ccbar, pt_min = 0.2, pt_max = 6.,step = 100) where {A<:SplineInterp,B<:SplineInterp, C<:SplineInterp, D<:SplineInterp}        
    fig1, ax1 = subplots()
    plot_params()
    set_spectra_axis(ax1, detector; pt_max = pt_max)
    
    pt = range(pt_min,pt_max,step) 
    det_params = detector_params(pt)    
    
    f_ME = MaximumEntropy.spectra_ME_2d(fo, lm_funct, particle;ccbar = ccbar, pt_min = pt_min, pt_max = pt_max, step = step)
    f_diff_ooe = Fluidum.spectra_analytic(fo,particle,fluidproperty.eos,pt_max=pt_max,pt_min=pt_min;delta_f=true)
    
    
    cons = det_params[detector.name].cons
    
    detector_name = string(detector.name)
    
    ax1.plot(pt,cons.*getindex.(f_ME[:],1),"#0C7BDC",label = string("ME"))
    ax1.plot(pt,cons.*getindex.(f_diff_ooe[:],1),"#FFC20A",label = string("Multiple Species"))
        
        
    if data == true       
        HC_data, HC_pT, HC_xerr, HC_yerrLow = read_hf(particle.name, string(exp_dir,detector_name),detector)
        ax1.errorbar(HC_pT,HC_data,xerr=0,yerr=HC_yerrLow,fmt="o",markersize=3,c="black",label= "$(detector_name) data")
    end

    ax1.legend(title = particle_legends[particle.name], loc="upper right") 
    return fig1 
end
