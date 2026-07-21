abstract type Lagr_Multiplier end

struct Lagr_Multiplier_1D{A<:Real} <: Lagr_Multiplier
    mult_nu::A
end

struct Lagr_Multiplier_1D_eq{A<:Real} <: Lagr_Multiplier
    mult_n::A
end

struct Lagr_Multiplier_2D_eq{A<:Real, B<:Real} <: Lagr_Multiplier
    mult_e::A
    mult_n::B
end


struct Lagr_Multiplier_2D{A<:Real, B<:Real} <: Lagr_Multiplier
    mult_n::A
    mult_nu::B
end

struct Lagr_Multiplier_3D{A<:Real, B<:Real, C<:Real} <: Lagr_Multiplier
    mult_e::A
    mult_n::B
    mult_nu::C
end



function Base.zero(::Type{MaximumEntropy.Lagr_Multiplier_1D{Float64}})
    MaximumEntropy.Lagr_Multiplier_1D(0.0)
end

function Base.:+(a::Lagr_Multiplier_1D{Float64}, b::Real) 
    MaximumEntropy.Lagr_Multiplier_1D(a.mult_nu + float(b))
end 

function ratio_bounds_sampled(T, ur; m, nη=120, nφ=120, npt=120,
                              etap_min=0.0, etap_max=10.0,
                              phip_min=0.0, phip_max=2π,
                              pt_min=0.0, pt_max=10.0)

    ut = sqrt(1 + ur^2)

    vals = Float64[]

    for etap in range(etap_min, etap_max, length=nη)
        for phip in range(phip_min, phip_max, length=nφ)
            for pt in range(pt_min, pt_max, length=npt)

                mt = sqrt(pt^2 + m^2)

                udotp = -ut * mt * cosh(etap) + ur * pt * cos(phip)

                ptransv = -ut * ur * mt * cosh(etap) +
                           ut^2 * pt * cos(phip)

                push!(vals, ptransv / (-udotp))
            end
        end
    end

    return extrema(vals)
end

function get_list(cfg; ccbar)
    hq = Fluidum.HadronResonanceGas_HQ(; name_file = cfg.particle_list, ccbar = ccbar)
    degs = hq.particle_list.Degeneracy
    masses = hq.particle_list.Mass
    qs = (hq.particle_list.Nc + hq.particle_list.Nac); 
    return degs, masses, qs
end



function fluidum_diffusion_fraction(T, fug, eos, part)
    q = 1
    n_i = Fluidum.free_hadron(T, fug, part.Degeneracy, q; m = part.Mass)
    norm = Fluidum.normalization(T, fug, eos)
    return n_i/norm
end


function total_ME_moments(T, ur, lm, eos; ccbar, pt_max = 20.0, rtol = 1e-5)
    n_tot = zeros(length(T))
    nu_tot = copy(n_tot)
    for a in eachindex(T)
        for part in eos.hq_eos.particle_list
            ME = charm_moments_one_part(T[a], ur[a], part, lm[a]; ccbar = ccbar, pt_max = pt_max, rtol = rtol)
            n_tot[a] += ME[1]
            nu_tot[a] += ME[2]
        end
    end
    return n_tot, nu_tot 
end 