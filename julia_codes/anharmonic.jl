
# write a description for the code and date it.

using QuantumOptics
using PyPlot

m = 1. #Mass value
omega = 1. #Traps strenght
lambda = 0.01 #anharmonic terms strenght

#Position basis setup:
xmin = -5
xmax = 5
points = 100
pos_basis = PositionBasis(xmin,xmax,points)

p = momentum(pos_basis)
x = position(pos_basis)
H_0 = p^2/(2*m) + 1/2*m*omega^2*dense(x^2) #Harmonic Hamiltonian
H_1 = H_0-lambda*dense(x^4) #Adding the anharmonic term

#Initial state:
x0 = 1.5
p0 = 0
sigma = 0.6
psi0 = gaussianstate(pos_basis, x0, p0, sigma);

T = [0:0.1:3;] #Time range

#Time evolotion for both Hamiltonians:
tout0, psi_t0 = timeevolution.schroedinger(T, psi0, H_0);
tout1, psi_t1 = timeevolution.schroedinger(T, psi0, H_1);

# we can separate plotting and do it in another file.

#Plotting points:
x_points = samplepoints(pos_basis)

#Initialize the figure and subplots:
figure(figsize=(10,10))
subplot(211)
title("Harmonic Trap",fontsize=15)
ylabel(L"| \Psi(t) |^2",fontsize=15)
#Plotting loop for first figure:
for i=1:size(T)[1]
    psi = psi_t0[i]
    n = abs.(psi.data).^2
    plot(x_points, n, "b", alpha=0.9*(float(i)/length(T))^8+0.1)
end
#Initialize second figure:
subplot(212)
title("Harmonic Trap with Weak Anharmonicity",fontsize=15)
ylabel(L"| \Psi(t) |^2",fontsize=15)
xlabel(L"x",fontsize=15)
#Plotting loop for the second figure
for i=1:size(T)[1]
    psi = psi_t1[i]
    n = abs.(psi.data).^2
    plot(x_points, n, "r", alpha=0.9*(float(i)/length(T))^8+0.1)
end

# save figure in a dedicated figures directory
# use an example figure and insert it in to the report.
#Save the subplots:
#Alternatively use command: gcf() to show
#savefig("anharmonic2.png")
