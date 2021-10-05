#Plotting script

using PyPlot
using JLD2

#Initialize the figure:
figure(figsize=(8,5))
title("Magnetization at Site 1 for XXZ Spin Hamiltonian")
ylabel(L"\langle S^z_1 \rangle",fontsize=15)
xlabel("t",fontsize=15)

#Open the data:
f = jldopen("../data/XXZmagne.jld2","r")
T = read(f,"T")
y1 = read(f,"y1")
close("../data/XXZmagne.jld2")

plot(T,y1,"-or")

#Save to latex templates pics directory
savefig("../report_latex/Pics/XXZmagne.png")
