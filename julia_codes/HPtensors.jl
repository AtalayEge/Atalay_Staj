using ITensors
using JLD2

#Set parameters to run MPS computation:
N = 100
cutoff = 1E-8
tau = 0.1
ttotal = 5.0

#The function for obtaining n_c expectation values:
function One_dim_spin_wave(N,cutoff,tau,ttotal)

    Nsteps = Int(ttotal/tau)
    Jx = -1
    Jz = -1
    s = siteinds("Boson",N)
    S = 1/2
    ID = op("Id",s[1])


    gates = ITensor[]
    #Make gates:
    for j=1:N-1
        s1 = s[j]
        s2 = s[j+1]
        b1 = op("a",s1)
        b1dag = op("adag",s1)
        b2 = op("a",s2)
        b2dag = op("adag",s2)

        hj = Jx*S*(b1*b2dag+b1dag*b2)
        Gj = exp(-1.0im * tau/2 * hj)
        push!(gates,Gj)
    end
    #Add also the reverse gates:
    append!(gates,reverse(gates))

    #Initial state as |010001..> randomly
    psi = productMPS(s, n -> n==50 ? "0" : "1")

    c = div(N,2)

    t1 = zeros(Nsteps)
    y1 = zeros(Nsteps)
    y2 = zeros(Nsteps)
    y3 = zeros(Nsteps)
    t = 0
    #Calculate the expectation values for n_c at given time
    for step=1:Nsteps
        psi = apply(gates, psi; cutoff=cutoff)
        t1[step] = t
        t += tau
        n_c = S*expect(psi,"Id";site_range=c:c)-expect(psi,"n";site_range=c:c)
        n_c1 = S*expect(psi,"Id";site_range=c+1:c+1)-
                 expect(psi,"n";site_range=c+1:c+1)
        n_c2 = S*expect(psi,"Id";site_range=c:c)-
                 expect(psi,"n";site_range=c-1:c-1)
        y1[step] = n_c
        y2[step] = n_c1
        y3[step] = n_c2
    end
    return t1,y1,y2,y3
end

t1,y1,y2,y3 = One_dim_spin_wave(N,cutoff,tau,ttotal)
jldsave("../data/HPtensors.jld2";t1,y1,y2,y3)
