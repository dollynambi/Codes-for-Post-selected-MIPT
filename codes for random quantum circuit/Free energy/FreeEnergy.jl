include("HaarRandomCircuit.jl")
using .HaarRandomCircuit
using Random, CSV, DataFrames


function main(p::Float64,L::Int64,r::Int64)
    T =2L
    Random.seed!(r)
    s = run_circuit_periodic_t(L,p,T)
    df_logpsum = DataFrame(logpsum = s)
    CSV.write("logp_$(p)_$(L)_$(r).csv", df_logpsum)
end


function run_circuit_periodic_t(L,p,T)

    ψ = zeros(ComplexF64,2^L)
    ψ[1] = 1

    if iseven(div(L,4))
        HaarRandomCircuit.fullcircuit!(ψ,L,p,T;BCs="periodic")
    else
        HaarRandomCircuit.fullcircuit_OE!(ψ,L,p,T;BCs="periodic")
    end

    if iseven(div(L,4))
        s = HaarRandomCircuit.fullcircuit!(ψ,L,p,6T;BCs="periodic")
    else
        s = HaarRandomCircuit.fullcircuit_OE!(ψ,L,p,6T;BCs="periodic")
    end

    return s
end


p = parse(Float64,ARGS[1])
L = parse(Int,ARGS[2])
r = parse(Int,ARGS[3])

main(p,L,r)

