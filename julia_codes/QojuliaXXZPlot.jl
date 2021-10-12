#Plotting script

using PyPlot
using JLD2

#Initialize the figure:
figure(figsize=(8,5))
title("Magnetization at Sites 1,2 and 3 for XXZ Spin Hamiltonian")
ylabel(L"\langle S^z \rangle",fontsize=15)
xlabel("t",fontsize=15)

#Open the data:
f = jldopen("../data/XXZmagne.jld2","r")
T = read(f,"T")
y1 = read(f,"y1")
y2 = read(f,"y2")
y3 = read(f,"y3")
close("../data/XXZmagne.jld2")

plot(T,y1,"-or",T,y2,"-ok",T,y3,"-*b")
legend([L"\langle S^z_1 \rangle",L"\langle S^z_2 \rangle",
        L"\langle S^z_3 \rangle"])

#Save to latex templates pics directory
savefig("../report_latex/Pics/XXZmagne.png")
