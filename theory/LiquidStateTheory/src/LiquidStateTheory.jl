module LiquidStateTheory
    using LinearAlgebra
    using Clapeyron, NLSolvers, ForwardDiff, GenericLinearAlgebra, Polynomials

    using Clapeyron, NLSolvers
    import Clapeyron: @registermodel, @comps, AssocOptions, has_sites, getparams
    import Clapeyron: N_A, k_B, e_c, ϵ_0
    import Clapeyron: sigma_LorentzBerthelot, epsilon_LorentzBerthelot, assoc_mix
    import Clapeyron: Solvers

    include("LS.jl")
    include("paLS.jl")
    include("sequences.jl")
    include("methods.jl")

    export paLS, LS, phase_eq_pure, phase_eq_pure_p, phase_eq_mix, electrochemical_potential_pure, μ_Π, dielectric_constant, crit_pure, crit_mix
    export electrostatic_F_and_grad, pairwise_distances_and_grads, coordinates
end