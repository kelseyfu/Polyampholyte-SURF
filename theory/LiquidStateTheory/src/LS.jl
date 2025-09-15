struct LSParam <: EoSParam
    N::SingleParam
    Z::SingleParam
end

abstract type LSModel <: EoSModel end

struct LS <: LSModel
    components::Vector{String}
    params::LSParam
end

LS
export LS

function LS(components;
    userlocations=nothing)
    params = getparams(components, [""]; userlocations=userlocations)
    N = params["N"]
    Z = params["Z"]

    packagedparams = LSParam(N,Z)
    model = LS(components, packagedparams)
    return model
end

@registermodel LS

function a(model::LSModel,T,ρ)
    # println(T)
    lB          = 1/T
    N           = model.params.N.values
    Z           = model.params.Z.values
        
    f0          = sum(@. ρ./N*(log(ρ./N)-1))
    η           = (π/6)*sum(ρ)
    fhs         = 6η^2*(4-3η)/(π*(1-η)^2)
    
    κ           = sqrt(4π*lB*sum(ρ.*Z.^2))
    Γ           = (-1+sqrt(1+2κ))/(2)
    fel         = -Γ^3*(2/3+Γ)/π
    
    ypp         = (2+η)/(2*(1-η)^2)*exp.(-lB.*Z/(1+Γ)^2+lB.*Z)
    fch         = sum(ρ./N.*((1 .- N).*log.(ypp)))
    return (f0+fhs+fel+fch)
end