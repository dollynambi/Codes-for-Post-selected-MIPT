include("HaarRandomCircuit.jl")
using .HaarRandomCircuit
using Random, CSV, DataFrames


function main(p::Float64,L::Int64,r::Int64)
    T =2L
    Random.seed!(r)
    MI,SA,SB,SC,SAB,SAC,SBC,SABC = run_circuit_periodic_t(L,p,T)
    df_tmi = DataFrame((MI),:auto)
    df_sa = DataFrame((SA),:auto)
    df_sb = DataFrame((SB),:auto)
    df_sc = DataFrame((SC),:auto)
    df_sab = DataFrame((SAB),:auto)
    df_sac = DataFrame((SAC),:auto)
    df_sbc = DataFrame((SBC),:auto)
    df_sabc = DataFrame((SABC),:auto)
    CSV.write("tmi_$(p)_$(L)_$(r).csv", df_tmi, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sa_$(p)_$(L)_$(r).csv", df_sa, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sb_$(p)_$(L)_$(r).csv", df_sb, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sc_$(p)_$(L)_$(r).csv", df_sc, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sab_$(p)_$(L)_$(r).csv", df_sab, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sac_$(p)_$(L)_$(r).csv", df_sac, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sbc_$(p)_$(L)_$(r).csv", df_sbc, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sabc_$(p)_$(L)_$(r).csv", df_sabc, header=["S0","S1","S2","S3","Inf"])
end


function run_circuit_periodic_t(L,p,T)

    MI = zeros(Float64,2T+1,5)
    SA = zeros(Float64,2T+1,5)
    SB = zeros(Float64,2T+1,5)
    SC = zeros(Float64,2T+1,5)
    SAB = zeros(Float64,2T+1,5)
    SAC = zeros(Float64,2T+1,5)
    SBC = zeros(Float64,2T+1,5)
    SABC = zeros(Float64,2T+1,5)

    ψ = zeros(ComplexF64,2^L)
    ψ[1] = 1

    if iseven(div(L,4))
        HaarRandomCircuit.fullcircuit!(ψ,L,p,T;BCs="periodic",MI,SA,SB,SC,SAB,SAC,SBC,SABC)
    else
        HaarRandomCircuit.fullcircuit_OE!(ψ,L,p,T;BCs="periodic",MI,SA,SB,SC,SAB,SAC,SBC,SABC)
    end

    return MI,SA,SB,SC,SAB,SAC,SBC,SABC
end



function run_circuit_periodic(L,p)

    ψ = HaarRandomCircuit.initialize_state(2^L)

    if iseven(div(L,4))
        HaarRandomCircuit.fullcircuit!(ψ,L,p,2L;BCs="periodic")
    else
        HaarRandomCircuit.fullcircuit_OE!(ψ,L,p,2L;BCs="periodic")
    end

    MI = zeros(Float64,5)
    SA = zeros(Float64,5)
    SB = zeros(Float64,5)
    SC = zeros(Float64,5)
    SAB = zeros(Float64,5)
    SAC = zeros(Float64,5)
    SBC = zeros(Float64,5)
    SABC = zeros(Float64,5)
    for jj=1:(div(L,2))
        A = collect(jj:(jj+div(L,4)-1))
        B = mod1.(collect((jj+div(L,4)):(jj+2*div(L,4)-1)),L)
        C = mod1.(collect((jj+2*div(L,4)):(jj+3*div(L,4)-1)),L)
        r = HaarRandomCircuit.top_mutual_information(ψ,A,B,C,L,[0,1,2,3,Inf])
        MI .+= r[1]
        SA .+= r[2]
        SB .+= r[3]
        SC .+= r[4]
        SAB .+= r[5]
        SAC .+= r[6]
        SBC .+= r[7]
        SABC .+= r[8]
    end
    MI ./= div(L,2)
    SA ./= div(L,2)
    SB ./= div(L,2)
    SC ./= div(L,2)
    SAB ./= div(L,2)
    SAC ./= div(L,2)
    SBC ./= div(L,2)
    SABC ./= div(L,2)
    return MI,SA,SB,SC,SAB,SAC,SBC,SABC
end


function multiple_runs(L::Int64,N::Int64,k::Int64)   # Here k is the division of job
    p = range(0,stop=0.3,length=61)
    np = length(p)
    MI = zeros(Float64,np,5)
    SA = zeros(Float64,np,5)
    SB = zeros(Float64,np,5)
    SC = zeros(Float64,np,5)
    SAB = zeros(Float64,np,5)
    SAC = zeros(Float64,np,5)
    SBC = zeros(Float64,np,5)
    SABC = zeros(Float64,np,5)
    for j in 1:np
        println("p = ",p[j])
        for i in 1:N
            Random.seed!(N*(k-1)+i)
            r = run_circuit_periodic(L,p[j])
            MI[j,:] += r[1]
            SA[j,:] += r[2]
            SB[j,:] += r[3]
            SC[j,:] += r[4]
            SAB[j,:] += r[5]
            SAC[j,:] += r[6]
            SBC[j,:] += r[7]
            SABC[j,:] += r[8]
            if any(isinf,r[2])
                print("p = ",p[j], "r = ",i)
                print(r[2])
            end
        end
    end
    MI ./= N
    SA ./= N
    SB ./= N
    SC ./= N
    SAB ./= N
    SAC ./= N
    SBC ./= N
    SABC ./= N
    df_tmi = DataFrame((MI),:auto)
    df_sa = DataFrame((SA),:auto)
    df_sb = DataFrame((SB),:auto)
    df_sc = DataFrame((SC),:auto)
    df_sab = DataFrame((SAB),:auto)
    df_sac = DataFrame((SAC),:auto)
    df_sbc = DataFrame((SBC),:auto)
    df_sabc = DataFrame((SABC),:auto)
    CSV.write("tmi_$(L)_$(k).csv", df_tmi, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sa_$(L)_$(k).csv", df_sa, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sb_$(L)_$(k).csv", df_sb, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sc_$(L)_$(k).csv", df_sc, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sab_$(L)_$(k).csv", df_sab, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sac_$(L)_$(k).csv", df_sac, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sbc_$(L)_$(k).csv", df_sbc, header=["S0","S1","S2","S3","Inf"])
    CSV.write("sabc_$(L)_$(k).csv", df_sabc, header=["S0","S1","S2","S3","Inf"])

end

p = parse(Float64,ARGS[1])
L = parse(Int,ARGS[2])
r = parse(Int,ARGS[3])

main(p,L,r)




