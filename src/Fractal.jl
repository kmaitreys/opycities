function mean_scat_t(
    lmd::Float64,
    R0::Float64,
    PN::Float64,
    df::Float64,
    k0::Float64,
    refrel::Complex{Float64},
    iqsca::Int64,
    iqcor::Int64,
    iqgeo::Int64,
    nang::Int64,
    iquiet::Bool,
)::Tuple{
    Float64, # c_ext: Extinction cross section
    Float64, # c_sca: Scattering cross section
    Float64, # c_abs: Absorption cross section
    Float64, # g_asym: Asymmetry parameter
    Float64, # dphi: Phase shift induced by aggregate
    Array{Float64,1}, # angs: Angle grid from 0 to 2π
    Array{Float64,2}, # smat: Scattering matrix elements
    Array{Float64,1}, # phase_function
}
    k = 2π / lmd
    Rg = R0 * (PN / k0)^(1.0 / df)
    Rc = sqrt(5.0 / 3.0) * Rg
    xg = k * Rg
    x0 = k * R0

    xstop = x0 + 4.0 * x0^(1.0 / 3.0) + 2.0
    nstop = round(Int, xstop)
    numax = nstop
    nmax = nstop

    if !iquiet
        println("Aggr. size param. = ", xg)
        println("Monomer size param. = ", x0)
        println("Characteristic radius = ", Rc)
        println("x_stop = ", xstop)
        println("nstop = ", nstop)
        println("numax = ", numax)
        println("nmax = ", nmax)
    end

    if iqsca != 1 && iqsca != 2 && iqsca != 3
        error("Method must be 1, 2 or 3")
    end

    if iqcor != 1 && iqcor != 2 && iqcor != 3
        error("Correlation function must be 1, 2 or 3")
    end

    if iqgeo != 1 && iqgeo != 2 && iqgeo != 3
        error("Geometric cross section must be 1, 2 or 3")
    end

    if nang <= 1
        error("Number of angles must be greater than 1")
    end

    if PN < 1.0
        error("Number of monomers must be greater than 1")
    end

    if df > 3.0
        error("Fractal dimension must be greater than 1")
    end

    if (numax + nmax) >= 500 && !iquiet
        println("Warning: numax + nmax >= 500")
    end

    ff = PN * (R0 / Rc)^3.0
    mgmref = mg_mixing(refrel, ff)
    dphic = 2.0 * k * Rc * abs(mgmref - 1.0)
    dphi0 = 2.0 * x0 * abs(refrel - 1.0)
    dphi = max(dphic, dphi0)







    c_ext = 0
    c_sca = 0
    c_abs = 0
    g_asym = 0
    dphi = 0
    angs = range(0, stop=2π, length=nang)
    smat = zeros(Float64, nang, 4)
    phase_function = zeros(Float64, nang)




    return c_ext, c_sca, c_abs, g_asym, dphi, angs, smat, phase_function
end


function mg_mixing(refrel::Complex{Float64}, f1::Float64)
    eps_1 = refrel * refrel
    eps_2 = complex(1.0)
    mg = eps_2 * (
        2.0 * f1 * (eps_1 - eps_2) + eps_1 + 2.0 * eps_2
    ) / (eps_1 + 2.0 * eps_2 - f1 * (eps_1 - eps_2))
    mgav = sqrt(mg)
    return mgav
end
