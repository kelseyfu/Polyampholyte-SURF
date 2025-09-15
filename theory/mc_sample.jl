
using Pkg, Revise, CSV, DataFrames
using Base.Threads
Pkg.activate("LiquidStateTheory")
using LiquidStateTheory, Statistics


# -- Utilities --

# Wrap to (-π, π], removing -pi/3 to pi/3
@inline function wrap_angle(x::Real)
    y = mod(x + π, 2π)
    if y < π/3
        return -π/3
    elseif y > 2π/3
        return π/3
    else
        return y - π
    end
end

# -- Metropolis sampler over α (N-2 angles) --

"""
    sample_angle_distribution(z; kwargs...) -> (edges, pdfs, α_samples, meta)

Run Metropolis MC to sample P(α | z) ∝ exp(-β U(α, z)), returning per-angle histograms.

Inputs (keywords):
- ℓ=1.0, Γ=1.0, lB=1.0, σ=1.0 : potential parameters
- β=1.0                        : inverse temperature
- nsteps=200_000               : total Metropolis steps
- burnin=20_000, thin=20       : burn-in and thinning
- step=0.25                    : proposal width (radians) for local angle updates
- min_sep=1                    : exclude near neighbors (2→exclude 1–2; 3→exclude 1–3)
- init=:straight or init=α0    : initial angles (vector or symbol)
- nbins=72                     : histogram bins over (-π, π]
- seed::Int = nothing          : RNG seed (optional)

Outputs:
- edges :: Vector{Float64}       # bin edges length nbins+1 over (-π, π]
- pdfs  :: Matrix{Float64}       # (N-2)×nbins, normalized per angle
- α_samples :: Matrix{Float64}   # nsamp × (N-2), the kept samples (post burn-in & thinning)
- meta :: NamedTuple             # acceptance statistics and means
"""
function sample_angle_distribution(z::AbstractVector{<:Real};
        ℓ::Real=1.0, Γ::Real=1.0, lB::Real=1.0, σ::Real=1.0, β::Real=1.0,
        nsteps::Int=200_000, burnin::Int=20_000, thin::Int=20,
        step::Real=0.25, min_sep::Int=1, init=:straight,
        nbins::Int=72, seed::Union{Nothing,Int}=nothing, eps_reg::Real=1e-12)

    seed === nothing || Random.seed!(seed)

    N  = length(z)
    NA = N - 2
    @assert NA ≥ 1 "Need at least 3 atoms (N≥3) so there is at least one angle."
    @assert burnin < nsteps

    # Initialize α
     if init === :straight
        α = ones(NA)*pi/4
        α[1:2:end] .*= -1
    elseif init isa AbstractVector
        @assert length(init) == NA
        Float64.(init)
    else
        error("init must be :straight or a vector of length N-2")
    end

    # Initial energy
    F, _ = electrostatic_F_and_grad(α, z; ℓ=ℓ, Γ=Γ, lB=lB, min_sep=min_sep, eps_reg=eps_reg)
    # println(F)
    # Storage
    nsamp = (nsteps - burnin) ÷ thin
    α_samples = Array{Float64}(undef, nsamp, NA)
    Rg = Array{Float64}(undef, nsamp)
    Re2e = Array{Float64}(undef, nsamp)
    A = Array{Float64}(undef, nsamp)
    accept = zeros(Int, NA)
    tried  = zeros(Int, NA)
    Q = zeros(Float64, nsamp)  # angle variances
    F_samp = zeros(Float64, nsamp)

    # Histogram setup over (-π, π]
    edges = range(-π, stop=π, length=nbins+1) |> collect
    counts = zeros(Float64, NA, nbins)

    # MC loop: local single-angle random-walk proposals
    sidx = 0
    @inbounds for t = 1:nsteps
        k = rand(1:NA)                           # pick which angle to move
        tried[k] += 1

        δp  = ((2rand() - 1) * pi/3  + 2pi/3)*step               # uniform in [-step, step]
        δm  = -((2rand() - 1) * pi/3  + 2pi/3)*step               # uniform in [-step, step]
        old = α[k]
        α[k]   = rand() < 0.5 ? δp : δm
        # α[k] = wrap_angle(old + δ)

        Fnew, _ = electrostatic_F_and_grad(α, z; ℓ=ℓ, Γ=Γ, lB=lB, min_sep=min_sep, eps_reg=eps_reg)
        ΔF = Fnew - F
        # println(ΔF)

        if (ΔF ≤ 0) || (rand() < exp(-β * ΔF / N))
            # accept
            F = Fnew
            accept[k] += 1
            
        else
            # reject → revert
            α[k] = old
        end

        # Record (thinned) samples after burn-in
        if t > burnin && ((t - burnin) % thin == 0)
            sidx += 1
            α_samples[sidx, :] = α
            r, _ = pairwise_distances_and_grads(α_samples[sidx, :])
            Rg[sidx] = sqrt(sum(r.^2) / (2*N^2))
            Re2e[sidx] = r[1,end]
            coords = coordinates(α_samples[sidx, :]; ℓ=ℓ)
            # Calculate the shape anisotropy
            rij1 = coords[:,1] .- coords[:,1]'
            rij2 = coords[:,2] .- coords[:,2]'
            R1 = sqrt(sum(rij1.^2) / (2*N^2))
            R2 = sqrt(sum(rij2.^2) / (2*N^2))
            A[sidx]  = 1/2*((R1^2-R2^2)^2/(Rg[sidx]^4))
            F_samp[sidx] = F
            Q[sidx] = sum(exp.(-β*(F_samp[1:sidx])/N))/sidx
        end
        
    end

    # Normalize histograms to PDFs (each row sums to 1)
    # pdfs = counts ./ sum(counts, dims=2)

    meta = (
        acceptance = accept ./ max.(tried, 1),
        mean_angles = vec(mean(α_samples; dims=1)),
        step = step, β = β, nsamples = nsamp
    )

    return α_samples, F_samp, Q, meta, Rg, Re2e, A
end

function parallel_sample_angle_distribution(z::AbstractVector{<:Real}; nchains::Int=4, kwargs...)
    results = Vector{Any}(undef, nchains)
    @threads for i in 1:nchains
        # Use a different seed for each chain for independence
        # seed = isnothing(kwargs[:seed]) ? rand(1:10^6) : kwargs[:seed] + i
        α_samples, F_samp, Q, meta, Rg, Re2e, A = sample_angle_distribution(z; kwargs...)
        results[i] = (α_samples, F_samp, Q, meta, Rg, Re2e, A)
    end
    return results
end

# Example usage:
nsteps = 2_000_000

sequences = [
 "+++++++++++++++++++++++++-------------------------",# 50/50   
 "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-",
 "+++++-----+++++-----+++++-----+++++-----+++++-----",
 "++++++++++++++++++++-+-+-+-+-+--------------------",
 "++++++++++++++-+-+-+-+-+-+-+-+-+-+-+--------------",
 "++++++++++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+----------",
 "++++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+----",
 "-------------------++++++++++++++++++++-+-+-+-+-+-",
 "-------------++++++++++++++-+-+-+-+-+-+-+-+-+-+-+-",
 "---------++++++++++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-",
 "---++++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-",
 "++++++++++++++++++++++++--------------------------",
 "++++++++++++++++++++++----------------------------",
 "++++++++++++++++++++------------------------------",
 "++++++++++++++++++--------------------------------",
 "++++++++++++++++----------------------------------",
 "++++++++++++++------------------------------------",
 "++++++++++++--------------------------------------",
 "++++++++++----------------------------------------",
 "++++++++------------------------------------------",
 "++++++--------------------------------------------",
 "++++----------------------------------------------",
 "++------------------------------------------------",
 "--------------------------------------------------",
 "+++++++++++++++++++-+-+-+-+-+---------------------",
 "+++++++++++++++++-+-+-+-+-+-----------------------",
 "+++++++++++++++-+-+-+-+-+-------------------------",
 "+++++++++++++-+-+-+-+-+---------------------------",
 "+++++++++++-+-+-+-+-+-----------------------------",
 "+++++++++-+-+-+-+-+-------------------------------",
 "+++++++-+-+-+-+-+---------------------------------",
 "+++++-+-+-+-+-+-----------------------------------",
 "+++-+-+-+-+-+-------------------------------------",
 "+-+-+-+-+-+---------------------------------------",
 "++++++++++++-+-+-+-+-+-+-+-+-+-+-+----------------",
 "++++++++++-+-+-+-+-+-+-+-+-+-+-+------------------",
 "++++++++-+-+-+-+-+-+-+-+-+-+-+--------------------",
 "++++-+-+-+-+-+-+-+-+-+-+-+------------------------",
 "++-+-+-+-+-+-+-+-+-+-+-+--------------------------",
 "++++++++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+------------",
 "++++++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--------------",
 "++++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+----------------",
 "++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+------------------",
 "++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+------",
 "++---++++-+++-+-++++++-++----+----+---+-++--+++--+",
 "-+-+-+++-+++++++-+----+++++-++--+++--+-----+--+++-",
 "-+---++---+----+-+-+-+++-+++--+++-++-++-+--++++-++",
 "+++---++---+-+--+-+-++--+-+-+++---++-+-++-++-----+",
 "+-+-+--+--+-+-+-+--+-++-+--+++++-+---++++-+----+-+",
 "+++-+-++++++--++-+---+++-+--+--+-+-+++---+-+++-+++",
 "--++++++-+--+++++-+++++-++++-+-+-+--+-++----+--+-+",
 "--+-+---++-+-+++--++-----++++--++--+++----++++-+--",
 "--++++----++--------+-+-++++---++++--+++++++++-+++",
 "++-+-++++--+---++---+++---+---+----+---+-+---++-++",
 "-----+++--+++++++-+----++----+++--+-------+++-+++-",
 "-+++++-----++-++-+-++++++-+--+--+----+--+--++-+-++",
 "--+++++-+-++--+-+-+++-+----++----+++-++----++---+-",
 "++++++++-++++-++-+------++--+--++++++++-+--+++--+-",
 "--+--++--+-+--+++---++++---+-++++-+---+-+-+-+-+++-",
 "-+++-+-+---++-++---+---++----+-++--++-+----+-++++-",
]


mkpath("data")
z = ones(50)
z[2:2:end] .*= -1
# z[21:2:25] .*= -1
# z[26:2:30] .*= -1

lB = LinRange(1,10,10)
Γ = LinRange(0,2,5)
sequence_name = "ABAB"

nchains = Threads.nthreads()  # Use all available threads

for (i, c) in enumerate(sequences)
    mkpath("data/sequ_$i")
    println("Processing sequence $i / $(length(sequences))")
    for j in 1:length(c)
        if c[j] == '+'
            z[j] = 1.0
        else
            z[j] = -1.0
        end
    end
    for j in 1:length(lB)
        for k in 1:length(Γ)
            mkpath("data/sequ_$i/lB_$(round(lB[j]; digits=1))_Γ_$(round(Γ[k]; digits=1))")
            results = parallel_sample_angle_distribution(z; Γ=0.0, lB=lB_i, nsteps=nsteps, step=0.5, nchains=nchains)
            α_samples = zeros(Int((nsteps-20_000)/20*nchains), length(z)-2)
            Rg = zeros(Int((nsteps-20_000)/20*nchains))
            Re2e = zeros(Int((nsteps-20_000)/20*nchains))
            A = zeros(Int((nsteps-20_000)/20*nchains))
            Q = zeros(Int((nsteps-20_000)/20*nchains))
            for (i, result) in enumerate(results)
                α_samples[Int((i-1)*((nsteps-20_000)/20))+1:i*Int((nsteps-20_000)/20), :] = result[1]
                Rg[Int((i-1)*((nsteps-20_000)/20))+1:i*Int((nsteps-20_000)/20)] = result[5]
                Re2e[Int((i-1)*((nsteps-20_000)/20))+1:i*Int((nsteps-20_000)/20)] = result[6]
                A[Int((i-1)*((nsteps-20_000)/20))+1:i*Int((nsteps-20_000)/20)] = result[7]
                Q[Int((i-1)*((nsteps-20_000)/20))+1:i*Int((nsteps-20_000)/20)] = result[3]
            end

            # Sample results in CSV
            df = DataFrame(α_samples, :auto)
            CSV.write("data/sequ_$i/α_samples_lB_$(round(lB_i; digits=1)).csv", df)

            df = DataFrame(hcat(Rg, Re2e, A, Q), [:Rg, :Re2e, :A, :Q])
            CSV.write("data/sequ_$i/Rg_Re2e_A_lB_$(round(lB_i; digits=1)).csv", df)
        end
    end
end