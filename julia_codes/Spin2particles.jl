#Data script

using QuantumOptics
using JLD2


J = -1 #Interaction strenght, +1 or -1 doesn't matter for this problem
spin_basis = SpinBasis(1//2) #Set the one particle basis
#Initialize the XX interaction Hamiltonian for two particles:
ID = identityoperator(spin_basis)
#S^x for particle number 1:
S_1 = tensor(sigmax(spin_basis),ID)
#S^x for particle number 2:
S_2 = tensor(ID,sigmax(spin_basis))
H_x = J*S_1*S_2 #Hamiltonian for the entire system
psi0 = tensor(spinup(spin_basis),spinup(spin_basis))#Initial state |00>

T = [0:0.1:3;] #Time evolution:
tout, psi_t = timeevolution.schroedinger(T, psi0, H_x)

#We can check how does the magnetization of first particles change over time
#Define S^z for the first particle:
S_1z = tensor(sigmaz(spin_basis),ID)
#Initialize the figure:

y = zeros(size(T)[1]) #Placeholder for the expectaton values
y2 = zeros(size(T)[1]) #Placeholder for density matrix method of
#Loop to compute expectation values:
for i=1:size(T)[1]
    psi = psi_t[i]
    rho1 = tensor(psi,dagger(psi)) #Density matrix method
    y2[i]=tr(rho1*S_1z)
    y[i] = expect(S_1z,psi) #Built in funtion
end

# Save the data in jld2 format:
# Use relative path
jldsave("../data/magnedata2.jld2";T,y,y2)
