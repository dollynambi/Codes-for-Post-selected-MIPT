include("HaarRandomCircuit.jl")
using .HaarRandomCircuit
using Random, CSV, DataFrames

function entropies(ψ::Vector{ComplexF64})
    MI = zeros(Float64,5)
    SAB = zeros(Float64,5)
    for jj=1:(div(L,2))
        A = collect(jj:(jj+div(L,4)-1))
        B = mod1.(collect((jj+div(L,4)):(jj+2*div(L,4)-1)),L)
        C = mod1.(collect((jj+2*div(L,4)):(jj+3*div(L,4)-1)),L)
        r = HaarRandomCircuit.top_mutual_information(ψ,A,B,C,L,[0,1,2,3,Inf])
        MI .+= r[1]
        SAB .+= r[5]
    end
    MI ./= div(L,2)
    SAB ./= div(L,2)
    return MI, SAB
end

function run_circuit_periodic(L,p,r)
    Random.seed!(r)
    ψ = HaarRandomCircuit.initialize_state(2^L)

    if iseven(div(L,4))
        HaarRandomCircuit.fullcircuit!(ψ,L,p,L;BCs="periodic")
    else
        HaarRandomCircuit.fullcircuit_OE!(ψ,L,p,L;BCs="periodic")
    end
    MI, SAB = entropies(ψ)
    df_tmi = DataFrame(transpose(MI),:auto)
    df_sab = DataFrame(transpose(SAB),:auto)
    CSV.write("tmi_$(p)_$(L)_$(r)_L.csv", df_tmi, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sab_$(p)_$(L)_$(r)_L.csv", df_sab, header=["S0","S1","S2","S3","Inf"])
    
    if iseven(div(L,4))
        HaarRandomCircuit.fullcircuit!(ψ,L,p,L;BCs="periodic")
    else
        HaarRandomCircuit.fullcircuit_OE!(ψ,L,p,L;BCs="periodic")
    end
    MI, SAB = entropies(ψ)
    df_tmi = DataFrame(transpose(MI),:auto)
    df_sab = DataFrame(transpose(SAB),:auto)
    CSV.write("tmi_$(p)_$(L)_$(r)_2L.csv", df_tmi, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sab_$(p)_$(L)_$(r)_2L.csv", df_sab, header=["S0","S1","S2","S3","Inf"])

end


p = parse(Float64,ARGS[1])
L = parse(Int,ARGS[2])
r = parse(Int,ARGS[3])
run_circuit_periodic(p,L,r)