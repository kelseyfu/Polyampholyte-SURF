struct paLSParam <: EoSParam
    N::SingleParam
    Z::SingleParam
end

abstract type paLSModel <: EoSModel end

struct paLS <: paLSModel
    components::Vector{String}
    sequence::Vector{Float64}
    params::paLSParam
end

paLS
export paLS

function paLS(components, sequence;
    userlocations=String[])
    params = getparams(components, [""]; userlocations=userlocations)
    N = params["N"]
    Z = params["Z"]
    packagedparams = paLSParam(N,Z)
    model = paLS(components, sequence, packagedparams)
    return model
end

@registermodel paLS

function a(model::paLSModel,lB,ρ; verbose=false)
    N           = model.params.N.values
    seq         = model.sequence
    α           = sum(seq[1:end-1].*seq[2:end].<0)/(N[1]-1)

    f0          = sum(@. ρ./N*(log(ρ./N)-1))
    η           = (π/6)*sum(ρ)
    fhs         = 6η^2*(4-3η)/(π*(1-η)^2)
    
    κ           = sqrt(4π*lB*sum(ρ))
    Γ           = (-1+sqrt(1+2κ))/(2)
    fel         = -Γ^3*(2/3+Γ)/π
    
    ypp         = (2+η)/(2*(1-η)^2)*exp.(-lB/(1+Γ)^2+lB)
    ypm         = (2+η)/(2*(1-η)^2)*exp.(lB/(1+Γ)^2-lB)
    fch        = -sum(ρ.*(N .-1)/N.*((1 .-α).*log.(ypp) .+ α.*log.(ypm)))

    ωGauss = [1.0597717569215825,1.112884701069448,1.2061552453369633,1.3357866280694308,1.4964744414361324,1.6816407769327313,1.883705139843864,2.0943951023931957,2.305085064942527,2.50714942785366,2.6923157633502592,2.8530035767169606,2.9826349594494275,3.0759055037169425,3.1290184478648086]
    wGauss = [0.01537662099805875,0.03518302374405389,0.05357961023358613,0.06978533896307719,0.08313460290849703,0.09308050000778112,0.09921574266355584,0.10128912096277995,0.09921574266355582,0.09308050000778122,0.08313460290849663,0.06978533896307768,0.053579610233586164,0.03518302374405393,0.015376620998058891]

    f_conf = zero(eltype(first(ρ)+lB))
    Q = zero(eltype(first(ρ)+lB))
    for i in 2:(N[1]-1)
        Hvals = [H(seq,i,lB,Γ,η,ω) for ω in ωGauss]
        Qi = sum(exp.(-Hvals).*wGauss)
        f_conf += -log.(Qi)*Qi
        Q += Qi
    end
    f_conf *= ρ[1]/Q/N[1]
    if verbose
        println("f0: ", f0)
        println("fhs: ", fhs)
        println("fel: ", fel)
        println("fch: ", fch)
        println("f_conf: ", f_conf)
    end

    return (f0+fhs+fel+fch+f_conf)
end

function Kcross(N,i,s1,s2,κ,r)
    α = 1/2*(s1+s2-s1^2/(i-1) - s2^2/(N-i))
    B = (s1/(i-1)+s2/(N-i))*r/2

    if κ == 0
        I_full = 1/B*(erfc(- B/(2sqrt(α))) - erfc(B/(2sqrt(α))))
        # Φ₊ = erf(1/(2sqrt(α))-B/(2sqrt(α))) - erf(- B/(2sqrt(α)))
        # Φ₋ = erf(1/(2sqrt(α))+B/(2sqrt(α))) - erf(+ B/(2sqrt(α)))
        # J_cut = 1/(4pi*B)*(Φ₊ - Φ₋)
        # println("I_full: ", I_full, ", J_cut: ", J_cut)
        return (I_full)/(4pi)
    else
        eακ2 = exp(α*κ^2)
        eBκ = exp(B*κ)
        I_full = eακ2/B*(erfc(sqrt(α)*κ - B/(2sqrt(α)))/eBκ - eBκ*erfc(sqrt(α)*κ + B/(2sqrt(α))))
        # Φ₊ = erf(1/(2sqrt(α))+sqrt(α)*κ-B/(2sqrt(α))) - erf(sqrt(α)*κ - B/(2sqrt(α)))
        # Φ₋ = erf(1/(2sqrt(α))+sqrt(α)*κ+B/(2sqrt(α))) - erf(sqrt(α)*κ + B/(2sqrt(α)))
        # J_cut = eακ2/(4pi*B)*(Φ₊/eBκ - eBκ*Φ₋)
        # println("I_full: ", I_full, ", J_cut: ", J_cut)
        return (I_full )/(4pi)
    end
end

function Kself(N,i,s1,s2,κ,r)
    α = abs(s1-s2)/6

    if κ == 0
        I_full = 1/sqrt(α*π)
        # Φ₊ = erf(1/(2sqrt(α))-B/(2sqrt(α))) - erf(- B/(2sqrt(α)))
        # Φ₋ = erf(1/(2sqrt(α))+B/(2sqrt(α))) - erf(+ B/(2sqrt(α)))
        # J_cut = 1/(4pi*B)*(Φ₊ - Φ₋)
        # println("I_full: ", I_full, ", J_cut: ", J_cut)
        return (I_full)/(4pi)
    else
        I_full = 1/sqrt(α*π)-κ*exp(α*κ^2)*erfc(sqrt(α)*κ)
        # Φ₊ = erf(1/(2sqrt(α))+sqrt(α)*κ-B/(2sqrt(α))) - erf(sqrt(α)*κ - B/(2sqrt(α)))
        # Φ₋ = erf(1/(2sqrt(α))+sqrt(α)*κ+B/(2sqrt(α))) - erf(sqrt(α)*κ + B/(2sqrt(α)))
        # J_cut = eακ2/(4pi*B)*(Φ₊/eBκ - eBκ*Φ₋)
        # println("I_full: ", I_full, ", J_cut: ", J_cut)
        return (I_full )/(4pi)
    end
end

function H(z,i,lB,κ,η,ω)
    N = length(z)
    r = sqrt(2*(1 - cos(ω)))
    ξ           = sin(ω/2)
    p1          = (-23.97695745797514,  195.4172655243969, -671.548442953145, 1210.8414381853227, -1263.8677912611697,  732.0871407433209, -181.9169471540193)
    p2          = (780.8573801610571, -6500.320617503649, 22545.43504610945, -41636.48654322127, 43139.52673720643, -23729.840872111734, 5403.381526901659)
    a, b    = 0.0, 0.0
    for i in 0:6
        a += p1[i+1] * ξ^i
        b += p2[i+1] * ξ^i
    end
    y3        = (1+η*a+η^2*b)/(1-η)^3#*exp(lB/(2*ξpp)*exp(-Γ*2*ξpp))
    S         = log(y3)

    A1 = z[i-1]*z[i+1]/r*exp(-κ*r)

    A2 = 0.0
    for j in 1:(i-1)
        for k in 1:(N-i)
            if j == (i-1) && k == (N-i)
                continue
            end
            A2 += z[j]*z[N-k+1]*Kcross(N,i,j,k,κ,r)
        end
    end

    for j in 1:(i-1)
        for k in (j+1):(i-1)
            A2 += z[j]*z[k]*Kself(N,i,j,k,κ,r)
        end
    end

    for j in 1:(N-i)
        for k in (j+1):(N-i)
            A2 += z[N-j+1]*z[N-k+1]*Kself(N,i,j,k,κ,r)
        end
    end
    return lB*(A1 + A2) - S
end