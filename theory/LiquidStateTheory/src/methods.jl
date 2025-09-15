function μ_Π(model,T,ρT,ρ)
    Z = model.params.Z.values

    fun(x)  = a(model,T,[x[1],x[1],ρT-2x[1]])
    df(x) = Clapeyron.ForwardDiff.derivative(fun,x)
    
    μ = df(ρ)
    f = fun(ρ)

    return μ, f-μ*ρ
end

function phase_eq_pure(model,T,ρT,x0=[-20.,-0.305])
    F = zeros(2)
    # x0 = log10.(x0)
    f!(F,x) = phaseqObj_chrg(F,model,T,ρT,exp10.(x[1]),exp10.(x[2]))
    results = Clapeyron.Solvers.nlsolve(f!,big.(x0))
    sol = Clapeyron.Solvers.x_sol(results)
    return exp10.(sol)
end

function phaseqObj_chrg(F,model,T,ρT,ρsup,ρco)
    μsup, Πsup = μ_Π(model,T,ρT,ρsup)
    μco, Πco = μ_Π(model,T,ρT,ρco)
    
    F[1] = μsup-μco
    
    F[2] = Πsup-Πco

    # Add mass balance
    return F
end

function phase_eq_mix(model,T,ρT,ρ_0,x0=[1e-10,1e-10,0.5,0.5,0.5])
    F = zeros(5)
    x0 = log10.(x0)
    f!(F,x) = phaseqObj_chrg_mix(F,model,T,ρT,exp10.(x[1:2]),exp10.(x[3:4]),exp10.(x[5]),ρ_0)
    results = Clapeyron.Solvers.nlsolve(f!,x0,TrustRegion(Newton(),NWI()))
    sol = Clapeyron.Solvers.x_sol(results)
    return exp10.(sol)
end

function phaseqObj_chrg_mix(F,model,T,ρT,ρsup,ρco,phi_sup,ρ_0)
    fun(x)  = a(model,T,[x[1],x[1]+x[2],x[2],ρT-2x[1]-2x[2]])
    df(x) = Clapeyron.ForwardDiff.gradient(fun,x)

    μsup = df(ρsup)
    μco  = df(ρco)

    F[1:2] .= μsup .- μco
    F[3] = fun(ρsup) - fun(ρco) - dot(μsup,ρsup) + dot(μco,ρco)
    F[4:5] .= phi_sup*ρsup+(1-phi_sup)*ρco - ρ_0 

    return F
end

function electrochemical_potential_pure(model,T,ρsup,ρco)
    fun(x)  = a(model,T,[x[1],x[2],x[3]])
    df(x) = Clapeyron.ForwardDiff.gradient(fun,x)
    
    μsup = df(ρsup)
    μco  = df(ρco)
    ψ = μsup .- μco
    return ψ
end

function phase_eq_pure_p(model,T,p,x0=[-20.,-0.305])
    F = zeros(4)
    # x0 = log10.(x0)
    f!(F,x) = phaseqObj_isobaric(F,model,T,p,exp10.(x[1:2]),exp10.(x[3:4]))
    results = Clapeyron.Solvers.nlsolve(f!,big.(x0))
    sol = Clapeyron.Solvers.x_sol(results)
    return exp10.(sol)
end

function phaseqObj_isobaric(F,model,T,p,ρsup,ρco)
    fun(x)  = a(model,T,[x[1],x[1],x[2]])
    df(x) = Clapeyron.ForwardDiff.gradient(fun,x)

    μsup = df(ρsup)
    μco  = df(ρco)
    Πsup = fun(ρsup) - dot(μsup,ρsup)
    Πco  = fun(ρco) - dot(μco,ρco)

    F[1:2] .= μsup .- μco
    F[3] = Πsup - Πco
    F[4] = Πsup - p

    return F
end

function crit_pure(model::Union{LSModel},ρT,x0=big.([0.004,1e-2]))
    # F = zeros(2)
    x0[2] = log10(x0[2])
    f!(F,x) = critObj_chrg(F,model,ρT,x[1],exp10(x[2]))
    results = Clapeyron.Solvers.nlsolve(f!,x0)
    sol = Clapeyron.Solvers.x_sol(results)
    return sol[1], exp10(sol[2])
end

function critObj_chrg(F,model,ρT,lBc,ρc)
    N = model.params.N.values
    Z = model.params.Z.values
    ν = abs(Z[1]/Z[2])

    fun(x)  = a(model,lBc,[x[1],x[1],ρT-2x[1]])
    df(x) = Clapeyron.ForwardDiff.derivative(fun,x)
    d2f(x) = Clapeyron.ForwardDiff.derivative(df,x)
    d3f(x) = Clapeyron.ForwardDiff.derivative(d2f,x)
        
    F[1] = d2f(ρc)
    
    F[2] = d3f(ρc)
    # println(F)
    return F
end

function crit_mix(model::Union{LSModel},ρT,T,x0=big.([0.004,1e-2]))
    # F = zeros(2)
    x0 = log10.(x0)
    f!(F,x) = critObj_mix(F,model,ρT,T,exp10.(x))
    results = Clapeyron.Solvers.nlsolve(f!,x0)
    sol = Clapeyron.Solvers.x_sol(results)
    return exp10.(sol)
end

function critObj_mix(F,model,ρT,lBc,ρc)
    N = model.params.N.values
    Z = model.params.Z.values

    fun(x)  = a(model,lBc,[x[1],x[1]+x[2],x[2],ρT-2x[1]-2x[2]])
    H(x) = ForwardDiff.hessian(fun,x) #∂A/∂zᵢ∂zⱼ == ∂A/∂zⱼ∂zᵢ
    L(x) = det(Symmetric(H(x)))
    dL(x) = ForwardDiff.gradient(L,x)
    HH = H(ρc)
    LL = det(HH)
    Mᵢ = @view(HH[end,:])
    Mᵢ .=  dL(ρc)
    MM = HH
    F[1] = LL
    F[2] = det(MM)
    # println(F)
    return F
end
