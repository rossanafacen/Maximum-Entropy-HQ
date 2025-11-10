

using LaTeXStrings

const particle_legends = Dict(
    "Dc1865zer" => L"$\mathrm{D^0}$",
    "Dc1869plu" => L"$\mathrm{D^+}$",
    "Dc2007zer" => L"$\mathrm{D^{*0}}$",
    "Dc2010plu" => L"$\mathrm{D^{*+}}$",
    "Ds1968plu" => L"$\mathrm{D_s^+}$",
    "Lc2286plu" => L"$\mathrm{\Lambda_c^+}$",
    "Om2695zer" => L"$\mathrm{\Omega_c^0}$",
    "Xic2470zer" => L"$\mathrm{\Xi_c^+}$",
    "jp3096zer" => L"$\mathrm{J/\Psi}$",
    "pi0139plu" => L"$\mathrm{\pi^+-\pi^-}$",
    "Ka0492plu" => L"$\mathrm{K^+-K^-}$",
    "pr0938plu" => L"$\mathrm{p-\bar{p}}$"
)



const fields_plotting = Dict(
    1 => L"$\mathrm{T \ [GeV]}$",
    2 => L"$u^r [GeV]$",
    3 => L"$\pi^\phi_\phi/s$",
    4 => L"$\pi^\eta_\eta/s$",
    5 => L"$\Pi/s$",
    6 => L"$\alpha$",
    7 => L"$\nu_q^r\,\mathrm{[fm^{-3}]}$",
    8 => L"$n_q\,\mathrm{[fm^{-3}]}$",
    9 => L"$n_q \times \tau \mathrm{[fm^{-2}]}$",
    10 => L"$\nu_q^r/(T \cdot norm)exp(\mu)$"
)

const initial_fields_plotting = Dict(
    1 => L"$\mathrm{T \ [GeV]}$",
    2 => L"$\alpha_q$",
    3 => L"$n_q\,\mathrm{[fm^{-3}]}$",
    4 => L"$\nu_q^r\,\mathrm{[fm^{-3}]}$"
)


function detector_params(pt)
return Dict(
    :ALICE => (
        ylabel = L"$dN/dp_\mathrm{T} dy \ [\mathrm{GeV}/c]^{-1}$",
        cons = 2*2π *pt,
    ),
    :RHIC => (
        ylabel = L"$d^2 N/2 \pi p_\mathrm{T} dp_\mathrm{T} dy \ [\mathrm{GeV}/c]^{-2}$",
        cons = 1.0,
    )
)
end

function plot_params(;gui = false)
    close("all")
    pygui(gui) 
    rc("text", usetex ="true") 
    rc("font", family = "serif") 
    rc("font", serif = "centaur") 
    rc("font", size = 17)
    rc("xtick", labelsize=20)
    rc("ytick", labelsize=20)
    rc("lines", markersize = 2)
    rc("legend", fontsize = 16)
    rc("axes", titlesize = 16)
    rc("axes", titlepad = 10)
    return nothing
end 

function set_spectra_axis(ax, detector; pt_min = 0.0, pt_max = 6.0, npt = 100, log = true)
    plot_params()
    pt = range(pt_min, pt_max, npt)
    params = detector_params(pt)    
    ax.set_xlim([pt_min,pt_max])

    if log == true
        ax.set_yscale("log")
    end
    ax.set_xlabel(L"$p_\mathrm{T} \ \mathrm{[GeV}/c]$")
    ax.set_ylabel(params[detector.name].ylabel)
    return nothing
end


function set_fields_axis(ax, n_fields)
    plot_params()
    ax.set_xlabel(L"\mathrm{r [fm]}")
    ax.set_ylabel(fields_plotting[n_fields])
    ax.tick_params(axis="both", which="major", labelsize=15)
    return nothing
end

function set_initial_fields_axis(ax, fig, n_fields; grid_params, label = false)
    plot_params()
    ax.set_xlabel(L"\mathrm{r [fm]}")
    ax.set_ylabel(initial_fields_plotting[n_fields])
    ax.tick_params(axis="both", which="major", labelsize=15)
    if label == true
        ax.text(grid_params.rmax/2.5, ax.get_ylim()[2]/2, L"Pb-Pb, $\sqrt{s_\mathrm{NN}} = 5.02 \ \mathrm{TeV}$", fontsize=18)
        ax.text(grid_params.rmax/2.5, ax.get_ylim()[2]/2.5, L"cc: 0-10\%", fontsize=18)
    end  
    
    if n_fields == 3
        ax.set_xlim(0, grid_params.rmax)
    end

    if n_fields == 4
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    end
    fig.subplots_adjust(left=0.18, bottom = 0.18) 
    display(fig)
    return nothing
end

function set_fo_fields_axis(ax, n_fields)
    plot_params()
    ax.set_xlabel(L"\alpha")
    ax.set_ylabel(fields_plotting[n_fields])
    ax.tick_params(axis="both", which="major", labelsize=15)
    return nothing
end


function negative_spectrum_check(f,pt)
    
for i in 1:size(f,1) 
        if f[i][1] < 0
            @warn("Spectrum becomes negative")
            print("pT after which the spectrum becomes negative:")
            @show pt[i]
            break
        end
    end 
end 
