include("HaarRandomCircuitAncilla.jl")
using .HaarRandomCircuitAncilla
using LinearAlgebra
using CSV, DataFrames, Random


function main(p::Float64,L::Int64,r::Int64)
    Random.seed!(r)
    s = run_periodic(p,L)
    df = DataFrame(s,:auto)
    CSV.write("ancilla_entropy_$(p)_$(L)_$(r).csv",df,header=["s1", "s2", "s3","sinf"]) 
end


function entangle_ancilla!(ψ::Vector{ComplexF64},L::Int64)

    ψ[1:2^L] .= normalize(randn(ComplexF64,2^L))./sqrt(2)
    random_vec = randn(ComplexF64, 2^L)
    Q, _ = qr([ψ[1:2^L] random_vec])  

    ψ[2^L+1:end] = normalize(Q[:,2])./sqrt(2)

end

function run_periodic(p::Float64,L::Int64)
    
    s = zeros(Float64,2L+1,4)
    T = 2L

    # Defining the initial state
    ψ = zeros(ComplexF64,2^(L+1))

    entangle_ancilla!(ψ,L)
    HaarRandomCircuitAncilla.full_evolution_ancilla!(ψ,T,p,L+1,s)
    
    return s
end

p = parse(Float64,ARGS[1])
L = parse(Int,ARGS[2])
r = parse(Int,ARGS[3])

main(p,L,r)

