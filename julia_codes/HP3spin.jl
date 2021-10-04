using QuantumOptics
using JLD2

N = 3 #Chains size
Jx = -1
Jz = -1
S = 1//2 #Spin of a single site
spin_basis = SpinBasis(S)
f_basis = FockBasis(Int(2*S))

#Function to create 2^N sized basis
function N_basis(basis,N)
    N_basis = tensor(basis,basis)
    for i=1:N-2
        N_basis = tensor(N_basis,basis)
    end
    return N_basis
end

Nbasis = N_basis(f_basis,N) #Basis for our problem

#Function to find Hamiltonian of H-P transformation
function One_dim_chain(N_basis,S,N)
    largeID = identityoperator(N_basis)
    H_hp = 0*largeID
    b = destroy(f_basis)
    b_dag = dagger(b)
    for i=1:N-1
        j = i+1
        b_i = embed(N_basis,i,b)
        b_dagi = dagger(b_i)
        b_j = embed(N_basis,j,b)
        b_dagj = dagger(b_j)
        n_i = b_dagi*b_i
        n_j = b_dagj*b_j
        H_hp += Jx*S*(b_i*b_dagj+b_dagi*b_j)
        H_hp += Jx*(-1/4)*(b_i*b_dagj*n_j+n_i*b_i*b_dagj+b_dagi*(n_i+n_j)b_j)
        H_hp += Jx*(1/(16*S))*(n_i*b_i*b_dagj*n_j+b_dagi*n_i*n_j*b_j)
        H_hp += Jz*(S^2*largeID-S*n_i-S*n_j+n_i*n_j)
    end
    return H_hp,largeID
end

#Call the function
H_hp,largeID = One_dim_chain(Nbasis,S,N)

#S1z opertor for H-P transformation by definiton
S1z = S*largeID - dagger(b1)*b1

#Create the initial state |010>
f_state = fockstate(f_basis,0)
f_state2 = fockstate(f_basis,1)
f_state3 = fockstate(f_basis,0)
init_state = tensor(f_state,tensor(f_state2,f_state3))

#Time evolution:
T = 0:0.1:10;
tout, psi_t = timeevolution.schroedinger(T, init_state, H_hp);

#Find expectation values of S1z
y1 = zeros(size(T)[1])
for i=1:size(T)[1]
    psi = psi_t[i]
    rho1 = tensor(psi,dagger(psi))
    y1[i]=tr(rho1*S1z)
end

#Save the data:
jldsave("../data/HPmagne.jld2";T,y1)
