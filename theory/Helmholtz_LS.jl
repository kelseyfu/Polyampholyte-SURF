using Clapeyron, Measurements, ForwardDiffOverMeasurements
import Clapeyron: @registermodel, AssocOptions, has_sites, getparams
import Clapeyron: N_A, k_B
import Clapeyron: sigma_LorentzBerthelot, epsilon_LorentzBerthelot, assoc_mix

struct LSParam <: EoSParam
    N::SingleParam
    Nu::SingleParam
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
    userlocations=String[])
    params = getparams(components, [""]; userlocations=userlocations)
    N = params["N"]
    Nu = params["Nu"]
    Z = params["Z"]
    packagedparams = LSParam(N,Nu,Z)
    model = LS(components, packagedparams)
    return model
end

@registermodel LS

function a(model::LSModel,lB,ρ)
    N           = model.params.N.values
    Nu          = model.params.Nu.values
        
    f0          = sum(@. ρ./N*(log(ρ./N)-1))
    η           = (π/6)*sum(ρ)
    fhs         = 6η^2*(4-3η)/(π*(1-η)^2)
    
    κ           = sqrt(4π*lB*sum(ρ))
    Γ           = (-1+sqrt(1+2κ))/(2)
    fel         = -Γ^3*(2/3+Γ)/π
    
    ypp         = (2+η)/(2*(1-η)^2)*exp.(-lB/(1+Γ)^2+lB)
    ypm         = (2+η)/(2*(1-η)^2)*exp.(lB/(1+Γ)^2-lB)
    fch         = sum(ρ./N.*((1 .+ Nu .- N).*log.(ypp) .- Nu.*log.(ypm)))
    return (f0+fhs+fel+fch)
end