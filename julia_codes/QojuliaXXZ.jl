using QuantumOptics
using JLD2

N = 3 #Chains size
Jx = -1
Jz = -1
S = 1//2 #Spin of a single site
spin_basis = SpinBasis(S)


function N_basis(basis,N)
    N_basis = tensor(basis,basis)
    for i=1:N-2
        N_basis = tensor(N_basis,basis)
    end
    return N_basis
end

Nbasis = N_basis(spin_basis,N)

#1D spin chains XXZ hamiltonian:
function One_dim_chain(N_basis,N)
    largeID = identityoperator(N_basis)
    H_xxz = largeID*0
    sigma_x = 1/2*sigmax(spin_basis)
    sigma_y = 1/2*sigmay(spin_basis)
    sigma_z = 1/2*sigmaz(spin_basis)
    for i=1:N-1
        j = i+1
        S_ix = embed(N_basis,i,sigma_x)
        S_iy = embed(N_basis,i,sigma_y)
        S_iz = embed(N_basis,i,sigma_z)
        S_jx = embed(N_basis,j,sigma_x)
        S_jy = embed(N_basis,j,sigma_y)
        S_jz = embed(N_basis,j,sigma_z)
        H_xxz+=Jx/2*(S_ix*S_jx+S_iy*S_jy)+Jz*S_iz*S_jz
    end
    return H_xxz
end

H = One_dim_chain(Nbasis,N)
S1z = embed(Nbasis,1,1/2*sigmaz(spin_basis))
S2z = embed(Nbasis,2,1/2*sigmaz(spin_basis))
S3z = embed(Nbasis,3,1/2*sigmaz(spin_basis))

#Initial state of |010>
init_state = spinup(spin_basis)⊗spindown(spin_basis)⊗spinup(spin_basis)
T = 0:0.1:10;
tout, psi_t = timeevolution.schroedinger(T, init_state, H);

#Expectation values of S1z
y1 = zeros(size(T)[1])
y2 = zeros(size(T)[1])
y3 = zeros(size(T)[1])
for i=1:size(T)[1]
    psi = psi_t[i]
    rho1 = tensor(psi,dagger(psi))
    y1[i]=tr(rho1*S1z)
    y2[i]=tr(rho1*S2z)
    y3[i]=tr(rho1*S3z)
end

jldsave("../data/XXZmagne.jld2";T,y1,y2,y3)
