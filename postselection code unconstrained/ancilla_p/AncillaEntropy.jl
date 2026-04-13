# This code computes the ancilla entropies at times = L,2L,3L,4L
# for a measurement rate p, no. of qubits: L

include("HaarRandomCircuitAncilla.jl")
using .HaarRandomCircuitAncilla
using LinearAlgebra
using CSV, DataFrames,Random


function main(p::Float64,L::Int64,r::Int64)
    Random.seed!(r)
    sL_2, sL, s3L_2, s2L  = run_periodic(p,L)
    dfL_2 = DataFrame(sL_2',:auto)
    dfL = DataFrame(sL',:auto)
    df3L_2 = DataFrame(s3L_2',:auto)
    df2L = DataFrame(s2L',:auto)
    CSV.write("ancilla_entropy(L_2)_$(p)_$(L)_$(r).csv",dfL_2,header=["s1", "s2", "s3","sinf"])
    CSV.write("ancilla_entropy(L)_$(p)_$(L)_$(r).csv",dfL,header=["s1", "s2", "s3","sinf"])
    CSV.write("ancilla_entropy(3L_2)_$(p)_$(L)_$(r).csv",df3L_2,header=["s1", "s2", "s3","sinf"])
    CSV.write("ancilla_entropy(2L)_$(p)_$(L)_$(r).csv",df2L,header=["s1", "s2", "s3","sinf"])
end

# Maximally entangles the ancilla with the whole system
function entangle_ancilla!(ψ::Vector{ComplexF64},L::Int64)

    ψ[1:2^L] .= normalize(randn(ComplexF64,2^L))./sqrt(2)
    random_vec = randn(ComplexF64, 2^L)
    Q, _ = qr([ψ[1:2^L] random_vec])  
    ψ[2^L+1:end] = normalize(Q[:,2])./sqrt(2)

end

# Circuit dynamics to compute the ancilla entropy
function run_periodic(p::Float64,L::Int64)
    
    sL_2 = zeros(Float64,4)
    s3L_2 = zeros(Float64,4)
    sL = zeros(Float64,4)
    s2L = zeros(Float64,4)
    T = Int(L/2)

    ψ = zeros(ComplexF64,2^(L+1))

    entangle_ancilla!(ψ,L)  # Initialize ψ 

    HaarRandomCircuitAncilla.full_evolution_ancilla!(ψ,T,p,L+1)
    sL_2 .= HaarRandomCircuitAncilla.ancilla_entropy(ψ)
    HaarRandomCircuitAncilla.full_evolution_ancilla!(ψ,T,p,L+1)
    sL .= HaarRandomCircuitAncilla.ancilla_entropy(ψ)
    HaarRandomCircuitAncilla.full_evolution_ancilla!(ψ,T,p,L+1)
    s3L_2 .= HaarRandomCircuitAncilla.ancilla_entropy(ψ)
    HaarRandomCircuitAncilla.full_evolution_ancilla!(ψ,T,p,L+1)
    s2L .= HaarRandomCircuitAncilla.ancilla_entropy(ψ)

    return sL_2, sL, s3L_2, s2L 
end


p = parse(Float64,ARGS[1])
L = parse(Int,ARGS[2])
r = parse(Int,ARGS[3])

main(p,L,r)

