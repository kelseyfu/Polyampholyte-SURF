# ---------- Pairwise geometry helpers (unchanged) ----------
"""
pairwise_distances_and_grads(α; ℓ=1.0)

Planar chain with equal bond length ℓ and N-2 signed turns α (radians).
Atoms: N = length(α) + 2.

Returns:
  r  :: Matrix{Float64}      # N×N, pairwise distances r[i,j]
  dr :: Array{Float64,3}     # N×N×(N-2), dr[i,j,k] = ∂ r_ij / ∂ α_k
"""
function pairwise_distances_and_grads(α::AbstractVector{<:Real}; ℓ::Real=1.0)
    NA = length(α); N = NA + 2; NB = N - 1
    φ = zeros(Float64, NB)
    @inbounds for p = 2:NB
        φ[p] = φ[p-1] + float(α[p-1])
    end
    c  = cos.(φ);  s  = sin.(φ)   # u_p
    rc = -s;       rs =  c        # Rot(u_p)

    Ux = zeros(Float64, NB+1);
    Uy = zeros(Float64, NB+1);
    Rx = zeros(Float64, NB+1);
    Ry = zeros(Float64, NB+1);
    @inbounds for p = 1:NB
        Ux[p+1] = Ux[p] + c[p]
        Uy[p+1] = Uy[p] + s[p]
        Rx[p+1] = Rx[p] + rc[p]
        Ry[p+1] = Ry[p] + rs[p]
    end

    r  = zeros(Float64, N, N)
    dr = zeros(Float64, N, N, NA)

    @inbounds for i = 1:N
        for j = i+1:N
            Rix = ℓ * (Ux[j] - Ux[i])
            Riy = ℓ * (Uy[j] - Uy[i])
            rij = hypot(Rix, Riy)
            r[i,j] = rij; r[j,i] = rij

            if rij > 0.0
                kmax = min(j-2, NA); p1 = j - 1
                for k = i:kmax
                    p0  = max(i, k+1)
                    dRx = ℓ * (Rx[p1+1] - Rx[p0])
                    dRy = ℓ * (Ry[p1+1] - Ry[p0])
                    drdα = (Rix*dRx + Riy*dRy) / rij
                    dr[i,j,k] = drdα; dr[j,i,k] = drdα
                end
            end
        end
    end
    return r, dr
end

# ---------- U(α): energy, gradient, Hessian ----------
"""
electrostatic_U_and_grad(α, z; ℓ=1.0, Γ=1.0, lB=1.0, min_sep=1, eps_reg=1e-12)

Now returns: U, gU, HU
  U   :: Float64
  gU  :: Vector{Float64}        (length N-2)
  HU  :: Matrix{Float64}        ((N-2)×(N-2), symmetric)

Potential per pair (i<j):
  f(r) = (lB*z[i]*z[j]/r) * exp(-Γ r)  +  (ℓ/r)^12
"""
function electrostatic_U_and_grad(α::AbstractVector{<:Real}, z::AbstractVector{<:Real};
                                  ℓ::Real=1.0, Γ::Real=1.0, lB::Real=1.0,
                                  min_sep::Int=1, eps_reg::Real=1e-12)

    NA = length(α); N = NA + 2; NB = N - 1
    @assert length(z) == N  "z must have length N = length(α)+2"
    @assert 1 ≤ min_sep ≤ NB

    # Bond orientations and prefix sums (as before)
    φ = zeros(Float64, NB)
    @inbounds for p = 2:NB
        φ[p] = φ[p-1] + float(α[p-1])
    end
    c  = cos.(φ);  s  = sin.(φ)     # u_p
    rc = -s;       rs =  c          # Rot(u_p)

    Ux = zeros(Float64, NB+1);
    Uy = zeros(Float64, NB+1);
    Rx = zeros(Float64, NB+1);
    Ry = zeros(Float64, NB+1);
    @inbounds for p = 1:NB
        Ux[p+1] = Ux[p] + c[p]
        Uy[p+1] = Uy[p] + s[p]
        Rx[p+1] = Rx[p] + rc[p]
        Ry[p+1] = Ry[p] + rs[p]
    end

    U  = 0.0
    gU = zeros(Float64, NA)
    HU = zeros(Float64, NA, NA)

    @inbounds for i = 1:N-1
        jmin = i + min_sep
        jmin > N && continue
        zi = float(z[i])
        for j = jmin:N
            zj = float(z[j])

            # Geometry for (i,j)
            p1 = j - 1
            Rix = ℓ * (Ux[j] - Ux[i]);  Riy = ℓ * (Uy[j] - Uy[i])
            r   = hypot(Rix, Riy)
            r̃   = max(r, eps_reg)
            r2  = r̃*r̃;  r3 = r2*r̃

            # Potential pieces and derivatives wrt r
            Aij = lB * zi * zj
            e   = exp(-Γ * r̃)

            # f'(r) = d/dr [A e^{-Γ r}/r + (ℓ/r)^12]
            fp   = Aij * e * (-Γ / r̃ - 1.0 / r2)  +  ( -8.0 * (ℓ / r̃)^8 / r̃ )
            # f''(r) = d^2/dr^2 [...]
            fpp  = Aij * e * (Γ^2 / r̃ + 2Γ / r2 + 2.0 / r3)  +  156.0 * (ℓ / r̃)^8 / (r2)

            # Accumulate energy
            U += Aij * e / r̃ + (ℓ / r̃)^8

            # Angles that affect (i,j): k ∈ [i, j-2] ∩ [1, NA]
            kmax = min(j-2, NA)
            if kmax >= i && r̃ > 0
                # Precompute dot with R for reuse
                for k = i:kmax
                    p0k = max(i, k+1)
                    dRxk = ℓ * (Rx[p1+1] - Rx[p0k])
                    dRyk = ℓ * (Ry[p1+1] - Ry[p0k])
                    RdotdRk = Rix*dRxk + Riy*dRyk
                    drk = RdotdRk / r̃
                    gU[k] += fp * drk

                    # Hessian block H_{k,t}
                    for t = k:kmax
                        p0t  = max(i, t+1)
                        dRxt = ℓ * (Rx[p1+1] - Rx[p0t])
                        dRyt = ℓ * (Ry[p1+1] - Ry[p0t])
                        RdotdRt = Rix*dRxt + Riy*dRyt
                        drt = RdotdRt / r̃

                        # Second derivative of r: r_{kt}
                        p02  = max(i, max(k,t)+1)
                        sumUx = ℓ * (Ux[p1+1] - Ux[p02])   # sum u over p=p02..p1
                        sumUy = ℓ * (Uy[p1+1] - Uy[p02])
                        Rdd   = - (Rix*sumUx + Riy*sumUy)

                        dRkdotdRt = dRxk*dRxt + dRyk*dRyt
                        r_kt = (dRkdotdRt + Rdd) / r̃ - (RdotdRk * RdotdRt) / r3

                        Hkt = fpp * drk * drt + fp * r_kt
                        HU[k,t] += Hkt
                        if t != k
                            HU[t,k] += Hkt
                        end
                    end
                end
            end
        end
    end
    return U, gU, HU
end

# ---------- S(α): value, gradient, Hessian (diagonal) ----------
const P1 = (-23.97695745797514,  195.4172655243969, -671.548442953145,
            1210.8414381853227, -1263.8677912611697,  732.0871407433209,
            -181.9169471540193)
const P2 = ( 780.8573801610571, -6500.320617503649, 22545.43504610945,
           -41636.48654322127,  43139.52673720643, -23729.840872111734,
             5403.381526901659)

@inline function _evalpoly_012(x, c::NTuple{N,T}) where {N,T}
    # Horner: returns (P, P', P'')
    S = promote_type(typeof(x), T)
    p  = S(c[N]); dp = zero(S); d2p = zero(S)
    @inbounds for k = N-1:-1:1
        d2p_new = muladd(d2p, x, 2*dp)
        dp_new  = muladd(dp,  x, p)
        p       = muladd(p,   x, S(c[k]))
        d2p = d2p_new; dp = dp_new
    end
    return p, dp, d2p
end

"""
electrostatic_S_and_grad(α; η=0.4)

Now returns: S, gS, HS
  S  :: T
  gS :: Vector{T}        (length N-2)
  HS :: Matrix{T}        (diagonal; off-diagonals are zero)
"""
function electrostatic_S_and_grad(α::AbstractVector{<:Real}; η::Real=0.4)
    T  = promote_type(eltype(α), typeof(η))
    ηT = T(η); oneT = one(T); halfT = oneT/2

    δ     = oneT - ηT
    invδ3 = inv(δ*δ*δ)

    p1T = ntuple(i->T(P1[i]), length(P1))
    p2T = ntuple(i->T(P2[i]), length(P2))

    S  = zero(T)
    gS = zeros(T, length(α))
    HS = zeros(T, length(α), length(α))

    @inbounds for i in eachindex(α)
        a2 = T(α[i]) / 2
        s  = sin(a2)           # sin(α/2)
        ξ  = cos(a2)           # cos(α/2)
        a, ap, app = _evalpoly_012(ξ, p1T)
        b, bp, bpp = _evalpoly_012(ξ, p2T)

        num = oneT + ηT*a + (ηT^2)*b      # numerator of g_i
        gi  = num * invδ3
        S  += log(gi)

        dξ_dα  = -halfT * s
        d2ξ_dα2 = - (oneT/4) * ξ

        num_p = (ηT*ap + (ηT^2)*bp) * dξ_dα
        num_pp = ηT*(app*(dξ_dα^2) + ap*d2ξ_dα2) +
                 (ηT^2)*(bpp*(dξ_dα^2) + bp*d2ξ_dα2)

        gS[i] = num_p / num
        HS[i,i] = (num_pp*num - num_p^2) / (num*num)
        # off-diagonals remain zero (angles independent in S)
    end
    return S, gS, HS
end

# ---------- Combined F = U - S ----------
function electrostatic_F_and_grad(α::AbstractVector{<:Real}, z::AbstractVector{<:Real};
                                  ℓ::Real=1.0, Γ::Real=1.0, lB::Real=1.0, η=0.4,
                                  min_sep::Int=1, eps_reg::Real=1e-12)
    U, gU, _ = electrostatic_U_and_grad(α, z; ℓ=ℓ, Γ=Γ, lB=lB,
                                        min_sep=min_sep, eps_reg=eps_reg)
    S, gS, _ = electrostatic_S_and_grad(α; η=η)
    return U - S, gU .- gS
end

# New: full (value, grad, Hessian) for F
function electrostatic_F_grad_hess(α::AbstractVector{<:Real}, z::AbstractVector{<:Real};
                                   ℓ::Real=1.0, Γ::Real=1.0, lB::Real=1.0, η=0.4,
                                   min_sep::Int=2, eps_reg::Real=1e-12)
    U, gU, HU = electrostatic_U_and_grad(α, z; ℓ=ℓ, Γ=Γ, lB=lB,
                                         min_sep=min_sep, eps_reg=eps_reg)
    S, gS, HS = electrostatic_S_and_grad(α; η=η)
    F  = U - S
    gF = gU .- gS
    HF = HU .- HS
    return F, gF, HF
end

function obj_fun(F, α, z; ℓ=1.0, Γ=1.0, lB=1.0, η=0.4,
                 min_sep=1, eps_reg=1e-12)
    U, dU, _ = electrostatic_F_and_grad(α, z; ℓ=ℓ, Γ=Γ, lB=lB, η=η,
              min_sep=min_sep, eps_reg=eps_reg)
    S, dS, _ = electrostatic_S_and_grad(α; η=η)
    
    F = dU.-dS
    # println(F)
    return F
end