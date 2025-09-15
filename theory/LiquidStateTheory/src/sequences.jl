"""
pairwise_distances_and_grads(α; ℓ=1.0)

Planar chain with equal bond length ℓ and N-2 signed turns α (radians).
Atoms: N = length(α) + 2.

Returns:
  r  :: Matrix{Float64}      # N×N, pairwise distances r[i,j]
  dr :: Array{Float64,3}     # N×N×(N-2), dr[i,j,k] = ∂ r_ij / ∂ α_k

Angles α_k: signed in-plane turn from bond k to bond k+1.
Bond 1’s orientation is fixed to 0.
"""
function pairwise_distances_and_grads(α::AbstractVector{<:Real}; ℓ::Real=1.0)
    NA = length(α)         # number of angles = N-2
    N  = NA + 2            # atoms
    NB = N - 1             # bonds

    # Bond orientations φ_p, p=1..NB; φ_1=0
    φ = zeros(Float64, NB)
    @inbounds for p = 2:NB
        φ[p] = φ[p-1] + float(α[p-1])
    end

    # Unit bond vectors u_p = (cosφ, sinφ) and +90° rotations
    c  = cos.(φ);  s  = sin.(φ)
    rc = -s;       rs =  c

    # Prefix sums (with zero at index 1) for bonds and rotations
    Ux = zeros(Float64, NB+1); Uy = similar(Ux)
    Uy .= 0.0
    Rx = zeros(Float64, NB+1); Ry = similar(Rx)
    Ry .= 0.0
    @inbounds for p = 1:NB
        Ux[p+1] = Ux[p] + c[p]
        Uy[p+1] = Uy[p] + s[p]
        Rx[p+1] = Rx[p] + rc[p]
        Ry[p+1] = Ry[p] + rs[p]
    end

    # Outputs
    r  = zeros(Float64, N, N)
    dr = zeros(Float64, N, N, NA)

    @inbounds for i = 1:N
        r[i,i] = 0.0
        for j = i+1:N
            # Bonds used: p = i .. j-1
            Rix = ℓ * (Ux[j] - Ux[i])
            Riy = ℓ * (Uy[j] - Uy[i])
            rij = hypot(Rix, Riy)
            r[i,j] = rij
            r[j,i] = rij

            if rij > 0.0
                # Angles that influence r_ij: α_k with k ∈ [i, j-2]
                kmax = min(j-2, NA)
                p1   = j - 1
                for k = i:kmax
                    p0  = max(i, k+1)
                    dRx = ℓ * (Rx[p1+1] - Rx[p0])
                    dRy = ℓ * (Ry[p1+1] - Ry[p0])
                    drdα = (Rix*dRx + Riy*dRy) / rij
                    dr[i,j,k] = drdα
                    dr[j,i,k] = drdα
                end
            end
        end
    end
    return r, dr
end

function coordinates(α::AbstractVector{<:Real}; ℓ::Real=1.0)
    N = length(α) + 2
    coords = zeros(Float64, N, 2)
    coords[1, :] .= [0.0, ℓ]  # First atom at origin

    # Compute coordinates iteratively
    for i = 3:N
        coords[i, 1] = coords[i-1, 1] + ℓ * cos(sum(α[1:i-2]))
        coords[i, 2] = coords[i-1, 2] + ℓ * sin(sum(α[1:i-2]))
    end

    return coords
end


"""
electrostatic_U_and_grad(α, z; ℓ=1.0, Γ=1.0, ℓB=1.0, min_sep=1, eps_reg=1e-12)

Planar chain with equal bond length ℓ.
Angles α are the N-2 signed in-plane turns at joints 1..N-2 (radians).
Charges z are on atoms 1..N (N = length(α) + 2).

Computes
    U = ∑_{i<j} (ℓB*z[i]*z[j]/r_ij) * exp(-Γ*r_ij)
and g[k] = ∂U/∂α_k for k=1..N-2.

`min_sep` lets you exclude near neighbors (e.g., 2 excludes 1–2; 3 excludes 1–2 and 1–3).
"""
function electrostatic_U_and_grad(α::AbstractVector{<:Real}, z::AbstractVector{<:Real};
                                  ℓ::Real=1.0, Γ::Real=1.0, lB::Real=1.0,
                                  min_sep::Int=1, eps_reg::Real=1e-12)

    NA = length(α)                  # number of turns = N-2
    N  = NA + 2                     # number of atoms
    NB = N - 1                      # number of bonds
    @assert length(z) == N  "z must have length N = length(α)+2"
    @assert 1 ≤ min_sep ≤ NB "min_sep must be in 1..(N-1)"

    # Bond orientations φ_p for p=1..NB; φ_1 = 0
    φ = zeros(typeof(α[1]),NB)
    @inbounds for p = 2:NB
        φ[p] = φ[p-1] + α[p-1]   # α indices 1..NB-1 = 1..N-2
    end

    # Unit bond vectors and their +90° rotations
    c = cos.(φ);  s = sin.(φ)           # u_p = (c[p], s[p])
    rc = -s;        rs =  c             # Rot(u_p) = (rc[p], rs[p])

    # Prefix sums (with leading zero) for fast segment sums
    Ux = zeros(typeof(α[1]), NB+1); Uy = similar(Ux)
    Uy .= 0
    Rx = zeros(typeof(α[1]), NB+1); Ry = similar(Rx)
    Ry .= 0
    @inbounds for p = 1:NB
        Ux[p+1] = Ux[p] + c[p]
        Uy[p+1] = Uy[p] + s[p]
        Rx[p+1] = Rx[p] + rc[p]
        Ry[p+1] = Ry[p] + rs[p]
    end

    U = 0.0
    dU = zeros(typeof(α[1]), NA)

    @inbounds for i = 1:N-1
        jmin = i + min_sep
        if jmin > N; continue; end
        zi = z[i]
        for j = jmin:N
            zj = z[j]

            # R_ij and r_ij
            Rix = ℓ * (Ux[j] - Ux[i])
            Riy = ℓ * (Uy[j] - Uy[i])
            r   = hypot(Rix, Riy)
            r_eff = max(r, eps_reg)

            # Pair energy and dU/dr for A_ij = ℓB z_i z_j, B_ij = Γ
            Aij  = lB * zi * zj
            e    = exp(-Γ * r_eff)
            U   += Aij * e / r_eff
            dUdr = Aij * e * (-Γ / r_eff - 1.0 / (r_eff^2))

            # Pure 12-power repulsion
            sr12 = (ℓ / r_eff)^12
            U   += sr12
            dUdr += -12.0 * sr12 / r_eff    # d/dr (σ/r)^12 = -12 (σ/r)^12 / r

            # Angles that influence (i,j) are k = i .. j-2, but also ≤ N-2
            p1 = j - 1
            kmax = min(j-2, NA)
            for k = i:kmax
                p0  = max(i, k+1)                 # sum Rot(u_p) for p=p0..p1
                dRx = ℓ * (Rx[p1+1] - Rx[p0])
                dRy = ℓ * (Ry[p1+1] - Ry[p0])
                drdα = (Rix*dRx + Riy*dRy) / r_eff
                dU[k] += dUdr * drdα
            end
        end
    end
    return U, dU
end

# Coefficients (degree 6) as constants
const P1 = (-23.97695745797514,  195.4172655243969, -671.548442953145,
            1210.8414381853227, -1263.8677912611697,  732.0871407433209,
            -181.9169471540193)
const P2 = ( 780.8573801610571, -6500.320617503649, 22545.43504610945,
           -41636.48654322127,  43139.52673720643, -23729.840872111734,
             5403.381526901659)

# Simultaneous polynomial and derivative evaluation (Horner)
@inline function _evalpoly_and_deriv(x, c::NTuple{N,T}) where {N,T}
    S = promote_type(typeof(x), T)
    p  = S(c[N])
    dp = zero(S)
    @inbounds for k = N-1:-1:1
        dp = muladd(dp, x, p)        # dp = dp*x + p
        p  = muladd(p,  x, S(c[k]))  # p  = p*x  + c[k]
    end
    return p, dp                     # returns (P(x), P'(x))
end
"""
    electrostatic_S_and_grad(α; η=0.4)

Compute
    S = ∑_i log g_i,  with  g_i = (1 + η a(ξ_i) + η^2 b(ξ_i)) / (1-η)^3,
where ξ_i = cos(α_i/2), and a, b are degree-6 polynomials with coefficients P1, P2.

Returns:
    S::T,  grad::Vector{T}   where grad[i] = ∂S/∂α_i

Notes:
- Uses the correct chain rule: dξ/dα = -½ sin(α/2).
- Equivalent to your original function but clearer and faster.
"""
function electrostatic_S_and_grad(α::AbstractVector{<:Real}; η::Real=0.4)
    T  = promote_type(eltype(α), typeof(η))
    ηT = T(η)
    oneT = one(T); halfT = oneT/2
    # Precompute the constant factor (1-η)^(-3)
    δ     = oneT - ηT
    invδ3 = inv(δ*δ*δ)

    # Promote coefficients once to T
    p1T = ntuple(i->T(P1[i]), length(P1))
    p2T = ntuple(i->T(P2[i]), length(P2))

    S = zero(T)
    dS = similar(α, T)
    dS .= 0

    @inbounds for i in eachindex(α)
        a2 = T(α[i]) / 2
        s = sin(a2)          # sin(α/2)
        ξ = cos(a2)          # = sin((π-α)/2)

        a, da_dξ = _evalpoly_and_deriv(ξ, p1T)
        b, db_dξ = _evalpoly_and_deriv(ξ, p2T)

        num  = oneT + ηT*a + (ηT^2)*b       # numerator of g_i
        gi   = num * invδ3                  # g_i
        S   += log(gi)

        # Chain rule: dξ/dα = -½ sin(α/2)
        dnum_dα = (ηT*da_dξ + (ηT^2)*db_dξ) * (-halfT*s)
        dS[i] = dnum_dα / num             # since ∂ log(gi)/∂α = (1/gi)*(∂gi/∂α) = (1/num)*(∂num/∂α)
    end

    return S, dS
end

function electrostatic_F_and_grad(α::AbstractVector{<:Real}, z::AbstractVector{<:Real};
                                  ℓ::Real=1.0, Γ::Real=1.0, lB::Real=1.0, η=0.4,
                                  min_sep::Int=1, eps_reg::Real=1e-12)
    U, dU = electrostatic_U_and_grad(α, z; ℓ=ℓ, Γ=Γ, lB=lB,
                                       min_sep=min_sep, eps_reg=eps_reg)
    S, dS = electrostatic_S_and_grad(α; η=η)
    return U - S, dU .- dS
end

function obj_fun(F, α, z; ℓ=1.0, Γ=1.0, lB=1.0, η=0.4,
                 min_sep=1, eps_reg=1e-12)
    U, dU = electrostatic_F_and_grad(α, z; ℓ=ℓ, Γ=Γ, lB=lB, η=η,
              min_sep=min_sep, eps_reg=eps_reg)
    println(α)
    F = dU
    println(F)
    return F
end