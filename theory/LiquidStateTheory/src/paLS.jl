struct paLSParam <: EoSParam
    N::SingleParam
    N2::SingleParam
    N3::SingleParam
    Z::SingleParam
end

abstract type paLSModel <: EoSModel end

struct paLS <: paLSModel
    components::Vector{String}
    params::paLSParam
end

paLS
export paLS

function paLS(components;
    userlocations=String[])
    params = getparams(components, [""]; userlocations=userlocations)
    N = params["N"]
    N2 = params["N2"]
    N3 = params["N3"]
    Z = params["Z"]
    packagedparams = paLSParam(N,N2,N3,Z)
    model = paLS(components, packagedparams)
    return model
end

@registermodel paLS

function a(model::paLSModel,lB,ρ)
    N           = model.params.N.values
    N2          = model.params.N2.values
    N3          = model.params.N3.values
    α           = @. N2/(N-1)
    β           = @. N3/(N-2)

    f0          = sum(@. ρ./N*(log(ρ./N)-1))
    η           = (π/6)*sum(ρ)
    fhs         = 6η^2*(4-3η)/(π*(1-η)^2)
    
    κ           = sqrt(4π*lB*sum(ρ))
    Γ           = (-1+sqrt(1+2κ))/(2)
    fel         = -Γ^3*(2/3+Γ)/π
    
    ypp         = (2+η)/(2*(1-η)^2)*exp.(-lB/(1+Γ)^2+lB)
    ypm         = (2+η)/(2*(1-η)^2)*exp.(lB/(1+Γ)^2-lB)
    fch2        = -sum(ρ.*(N .-1)/N.*((1 .-α).*log.(ypp) .+ α.*log.(ypm)))

    ωpp         = bond_angle(Γ,lB)
    ξpp         = sin(ωpp/2)
    ωpm         = bond_angle(Γ,-lB)
    ξpm         = sin(ωpm/2)
    p1          = (-23.97695745797514,  195.4172655243969, -671.548442953145, 1210.8414381853227, -1263.8677912611697,  732.0871407433209, -181.9169471540193)
    p2          = (780.8573801610571, -6500.320617503649, 22545.43504610945, -41636.48654322127, 43139.52673720643, -23729.840872111734, 5403.381526901659)
    app, bpp    = 0.0, 0.0
    apm, bpm    = 0.0, 0.0
    for i in 0:6
        app += p1[i+1] * ξpp^i
        bpp += p2[i+1] * ξpp^i
        apm += p1[i+1] * ξpm^i
        bpm += p2[i+1] * ξpm^i
    end
    y3pp        = (1+η*app+η^2*bpp)/(1-η)^3#*exp(lB/(2*ξpp)*exp(-Γ*2*ξpp))
    y3pm        = (1+η*apm+η^2*bpm)/(1-η)^3#*exp(-lB/(2*ξpm)*exp(-Γ*2*ξpm))
    fch3        = -sum(ρ.*(N .-2)/N.*((1 .-β).*log.(y3pp) .+ β.*log.(y3pm)))
    return (f0+fhs+fel+fch2+fch3)
end

function bond_angle(Γ,lB)
    _I1(x) = exp(-lB/(2*sin(x/2))*exp(-Γ*2*sin(x/2)))
    _I2(x) = x*exp(-lB/(2*sin(x/2))*exp(-Γ*2*sin(x/2)))

    # Integrate both equations using Gauss-Lagrange quadrature from Pi/3 to Pi

    w = (0.0666713443086882, 0.1494513491505809, 0.2190863625159823, 0.26926671930999574, 0.29552422471475276, 0.2955242247147535, 0.26926671930999657, 0.21908636251598212, 0.1494513491505805, 0.06667134430868825)
    x = (1.0745225706356332, 1.1885028631666061, 1.3829190662109194, 1.6405445069611637, 1.9384942591756193, 2.250295945610772, 2.548245697825227, 2.805871138575471, 3.0002873416197846, 3.114267634150757)
    I1 = sum(_I1(x[i]) * w[i] for i in 1:10)
    I2 = sum(_I2(x[i]) * w[i] for i in 1:10)
    return I2/I1
end