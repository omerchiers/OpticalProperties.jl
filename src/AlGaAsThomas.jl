module AlGaAsThomas

using OpticalProperties
using MyPhysicalConstants
const OP = OpticalProperties
const hbar = ħ


function epsAlGaAs(x, w, T, N, Nsign, mx, mu_max)
    # EPSGAAS Computes the dielectric function of GaAs
    #mx = m* masse effective du porteur
    #ones(1,length(N)) Multicouches et dopage
    #weq=w/(2*pi*100*c) Transition spectroscopique en rad/struct
    #mu : mobilité - environ arctan avec valeur saturation max et min
    #x fraction molaire en alu
    N = 1 #Une seule couche
    Convspectrotorad = (2.0 * pi * 100 * c0)
    weq = w / Convspectrotorad * ones(1, length(N))
    #weq = weq[1]
    # Lorentz
    epsInf = 10.9 - 2.41 * x
    wLOG = 292.2 - 36.4 * x + 2.5 * x^2
    wTOG = 269.2 - 9.2 * x - 11.6 * x^2
    wLOA = 358.2 - 79.7 * x - 40.3 * x^2
    wTOA = 358.2 + 8.9 * x - 4.2 * x^2
    gG = 3.3 + 11.5 * x
    gA = 9.9 - 2.35 * x
    # Drude
    particle_mass = mx[1] * (1 - Nsign) / 2 + mx[2] * (Nsign + 1) / 2
    wp = ones(length(weq), 1) * sqrt(N * electron^2 / (particle_mass * epsInf * epsilon0)) / Convspectrotorad
    #wp = wp[1]
    gammap = ones(length(weq), 1) * (electron / (particle_mass .* mu_max)) / Convspectrotorad
    eps1 = epsInf * ((wLOG^2 - wTOG^2) ./ (wTOG^2 .- (weq .^ 2) - 1im * weq * gG) .+ (wLOA^2 - wTOA^2) ./ (wTOA^2 .- (weq .^ 2) - 1im * weq * gA) .- (wp .^ 2) ./ (weq .* (weq .+ 1im * gammap)))
    # Interband
    (eps2) = epsAlGaAsSb_GonzalezCuevasrevised(x, 1, w, T)
    #(eps2,Eg)= (0.0,0.0)
    # Total
    eps = eps1[1] #+eps2*ones(1,length(N));
    return (eps, eps2, eps + eps2)
end

function mobilityCaugheyThomas(N, T, mu_min, mu_max, Nref, alpha, theta1, theta2)
    # MOBILITYCAUGHEYTHOMAS Computes the particle mobility, using Caughey-Thomas model
    return mu_min + (mu_max .* (300.0 / T) .^ theta1 - mu_min) ./ (1.0 + (N ./ (Nref .* (T / 300.0) .^ theta2)) .^ alpha)
end

struct AlGaAsStruc         #AlGaAs(x)
    depsr::Float64
    epsinf::Float64
    mx::Vector{Float64}
    mu_min::Vector{Float64}
    mu_max::Vector{Float64}
    Nref::Vector{Float64}
    alpha::Vector{Float64}
    theta1::Vector{Float64}
    theta2::Vector{Float64}
end

function initializestructAlGaAs(x)
    depsr = 12.9 - 2.84 * x
    epsinf = 10.9 - 2.41 * x
    mx = [0.063 + 0.083 * x; 0.51 + 0.25 * x] * m0
    mu_min = [(500 * ((0.063 * m0 / mx[1]))^1.5) * (1.0 / 10.89 - 1.0 / 12.9) / (1.0 / epsinf - 1.0 / depsr); 57.1 * x^2 - 67.1 * x + 20]
    mu_max = [(9400 * ((0.063 * m0 / mx[1]))^1.5) * (1.0 / 10.89 - 1.0 / 12.9) / (1.0 / epsinf - 1.0 / depsr); 781 * x^2 - 1073 * x + 491.5]
    Nref = [(6.0e16)^(1 - x) * (5.46e17)^x, 10^(1.40 * x^2 - 0.988 * x + 17.17)]
    alpha = [0.394 + (1 - 0.394) * x, 0.421 * x .^ 2 - 0.313 * x + 0.38]
    theta1 = ((1 - x) * [2.1; 2.2] + x * [2.1; 2.24]) ./ (1 + x .* (1 - x))
    theta2 = [3.0; 3.0]
    #Structures ont besoin d'être initialisées
    #mu=mobilityCaugheyThomas(d.Ndop'/1e6,d.T,mu_min,mu_max,Nref,alpha,theta1,theta2)/1e4
    AlGaAs = AlGaAsStruc(depsr,
        epsinf,
        mx,
        mu_min,
        mu_max,
        Nref,
        alpha,
        theta1,
        theta2
    )
    #testeyy.epsr = 12.9-2.84*x
    return AlGaAs
end






function epsAlGaAsSb_GonzalezCuevasrevised(x, y, w, T)
    #EPSALGAASSB_GONZALEZCUEVAS Computes the permittivity of Al(x)Ga(1-x)As(y)Sb(1-y), as a function of x and y (Gonzalez-Cuevas, 2017)

    ## Section 1: Loading data
    #GaAs
    dGaAs0 = [1.52 1.85 3.04 3.27 5.13 1.82 5.638 5.52 2.89 21.32 1e-3 1e-6 0.6631 0.037
        5.5e-4 3.5e-4 7.2e-4 7.2e-4 6.6e-4 6.1e-4 3.88e-5 0 0 0 0 6.6e-5 4.1e-4 0
        225 225 205 205 43 204 0 0 0 0 0 0 0 0]
    dGaAsT = @view(dGaAs0[1, :])
    dGaAsT[1:6] = @view(dGaAs0[1, 1:6]) - @view(dGaAs0[2, 1:6]) * T^2 ./ (T .+ @view(dGaAs0[3, 1:6]))
    dGaAsT[7] = dGaAs0[1, 7] + T * dGaAs0[2, 7]
    dGaAsT[12:13] = @view(dGaAs0[1, 12:13]) + T * @view(dGaAs0[2, 12:13])



    #AlAs
    dAlAs0 = [3.09 3.37 3.98 4.18 4.86 2.46 5.6513 23.74 2.24 0.705 1e-3 6e-6 0.7421 0.49
        8.8e-4 8.8e-4 6.7e-4 6.7e-4 3.2e-4 6.1e-4 2.9e-5 0 0 0 0 2.9e-5 1.64e-4 0
        530 530 0 0 0 204 0 0 0 0 0 0 0 0]

    dAlAsT = @view(dAlAs0[1, :])
    dAlAsT[1:6] = @view(dAlAs0[1, 1:6]) - @view(dAlAs0[2, 1:6]) * T^2 ./ (T .+ @view(dAlAs0[3, 1:6]))
    dAlAsT[7] = dAlAs0[1, 7] + T * dAlAs0[2, 7]
    dAlAsT[12:13] = @view(dAlAs0[1, 12:13]) + T * @view(dAlAs0[2, 12:13])

    #=
    #GaSb
    dGaSb0=[0.81 1.57 2.19 2.62 4.32 0.88 6.086 1.102 3.34 4.93 1e-3 1e-6 0.6804 0.02;
            4.2e-4 4.2e-4 6.8e-4 6.7e-4 9e-4 6e-4 4.72e-5 0 0 0 0 1.24e-4 4.8e-4 0;
            140 140 147 176 376 140 0 0 0 0 0 0 0 0]
    dGaSbT = dGaSb0[1,:]
    dGaSbT[1:6] = dGaSb0[1,1:6] - dAlAs0[2,1:6] * T^2 ./(T .+ dGaSb0[3,1:6])
    dGaSbT[7] = dGaSb0[1,7] + T * dGaSb0[2,7]
    dGaSbT[12:13] = dGaSb0[1,12:13] + T * dGaSb0[2,12:13]

    #AlSb
    dAlSb0=[2.39 3.06 2.94 3.43 4.18 1.7 6.122 36.58 1.6 1.19 1e-3 1e-6 0.693 1.94;
            4.2e-4 4.2e-4 4.7e-4 4.3e-4 4.7e-4 3.9e-4 2.6e-5 0 0 0 0 5.7e-5 3.06e-4 0;
            140 140 0 0 0 140 0 0 0 0 0 0 0 0]
    dAlSbT = dAlSb0[1,:]
    dAlSbT[1:6] = dAlSb0[1,1:6] - dAlSb0[2,1:6] * T^2 ./(T .+ dAlSb0[3,1:6])
    dAlSbT[7] = dAlSb0[1,7] + T*dAlSb0[2,7]
    dAlSbT[12:13] = dAlSb0[1,12:13] + T*dAlSb0[2,12:13]

    =#
    ## Section 2: Interpolation
    # GaSb
    if x == 0 && y == 0
        dAlGaAsSbT = dGaSbT

        # GaAs
    elseif x == 0 && y == 1
        dAlGaAsSbT = dGaAsT

        # AlSb
    elseif x == 1 && y == 0
        dAlGaAsSbT = dAlSbT

        # AlAs
    elseif x == 1 && y == 1
        dAlGaAsSbT = dAlAsT

    else

        # Interpolation : energy levels
        C = [0.37; 0.07; 0.45; 0.00; 0.02; 0.06]
        dAlGaAsT = x * dAlAsT[1:6] + (1 - x) * dGaAsT[1:6] - x * (1 - x) * C
        #=
        C = [0.72; 0.15; 0.00; 0.00; 0.00; 0.28]
        dAlAsSbT = y * dAlAsT[1:6] + (1-y) * dAlSbT[1:6] - y * (1-y) * C
        C=[1.20; 0.61; 0.00; 0.00; 0.00; 1.09]
        dGaAsSbT = y * dGaAsT[1:6] + (1-y) * dGaSbT[1:6] - y * (1-y) * C
        C = [0.69; 0.30; 0.28; 0.32; 0.00; 0.55]
        dAlGaSbT = x * dAlSbT[1:6] + (1-x) * dGaSbT[1:6] - x * (1-x) *C
        =#
        dAlGaAsSbT = 1.0 / (x * (1 - x) + y * (1 - y)) * (x * (1 - x) * (y * dAlGaAsT))

        # Interpolation : other parameters
        C = [0.002 0.012 0.019 0.010 -0.004 0.033 -0.031 0.008
            0.010 0.010 0.010 0.010 0.010 0.009 0.009 0.010]
        dAlGaAsSbT = [
            dAlGaAsSbT
            x * (y * @view(dAlAsT[7:14])) +
            (1 - x) * (y * @view(dGaAsT[7:14])) +
            @view(C[1, :]) * x * (1 - x)
        ]

    end

    ## Section 3: Computation of the dielectric function
    E = ħ * w / electron
    Eg = dAlGaAsSbT[1] * electron

    X0 = (E + im * dAlGaAsSbT[11]) / dAlGaAsSbT[1]
    Xs0 = (E + im * dAlGaAsSbT[11]) / dAlGaAsSbT[2]
    X1 = (E + im * dAlGaAsSbT[12]) / dAlGaAsSbT[3]
    Xs1 = (E + im * dAlGaAsSbT[12]) / dAlGaAsSbT[4]
    Xg = (E + im * dAlGaAsSbT[14]) / dAlGaAsSbT[6]
    fX0 = X0 .^ (-2) .* (2 - sqrt(1 + X0) - sqrt(1 - X0))
    fXs0 = Xs0 .^ (-2) .* (2 - sqrt(1 + Xs0) - sqrt(1 - Xs0))

    B1 = 44 * (dAlGaAsSbT[3] + 1 / 3 * (dAlGaAsSbT[4] - dAlGaAsSbT[3])) / (dAlGaAsSbT[7] * dAlGaAsSbT[3]^2)
    B2 = 44 * (dAlGaAsSbT[3] + 2 / 3 * (dAlGaAsSbT[4] - dAlGaAsSbT[3])) / (dAlGaAsSbT[7] * dAlGaAsSbT[4]^2)

    eps1 = dAlGaAsSbT[8] * dAlGaAsSbT[1]^(-1.5) * (fX0 + 0.5 * (dAlGaAsSbT[1] / dAlGaAsSbT[2])^(1.5) * fXs0)
    eps2 = -B1 .* X1 .^ (-2) .* log(1 - X1 .^ 2) - B2 .* Xs1 .^ (-2) .* log(1 - Xs1 .^ 2)
    eps3 = dAlGaAsSbT[9] * dAlGaAsSbT[5]^2 ./ (dAlGaAsSbT[5]^2 - E .^ 2 - im * E * dAlGaAsSbT[13])
    eps4 = 2 * dAlGaAsSbT[10] / pi * (-Xg .^ (-2) * log(dAlGaAsSbT[3] / dAlGaAsSbT[6]) +
                                      1 / 2 * (1 + Xg .^ (-1)) .^ 2 .* log((Xg + dAlGaAsSbT[3] / dAlGaAsSbT[6]) ./ (Xg + 1)) +
                                      1 / 2 * (1 - Xg .^ (-1)) .^ 2 .* log((Xg - dAlGaAsSbT[3] / dAlGaAsSbT[6]) ./ (Xg - 1)))
    eps = eps1 + eps2 + eps3 + eps4

    return eps
end



function PermittivityAlGaAs(x, w)
    AlgaaaaaAs = initializestructAlGaAs(x)
    (tabintra, tabinter, toconv) = epsAlGaAs(x, w, 300, 1, 1, AlgaaaaaAs.mx, AlgaaaaaAs.mu_max)
end

struct AlGaAs <: OptProp
    x :: Float64
    T :: Float64
end


function initializemolTAlGaAs(x, T)
    return (AlGaAs(x, T))
end

#Layer(initializemolTAlGaAs(0.0,300.0), 1e-4)

function permittivity(material::AlGaAs, w)
    AlgaaaaaAs = initializestructAlGaAs(material.x)
    (intra, inter, tot) = epsAlGaAs(material.x, w, material.T, 1, 1, AlgaaaaaAs.mx, AlgaaaaaAs.mu_max)
    return  tot
end

end #module