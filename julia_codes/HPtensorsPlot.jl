#Plotting script

using PyPlot
using JLD2

#Initialize the figure:
figure(figsize=(8,5))
title("Expectation Value Change for N=100 H-P Chain")
ylabel(L"\langle n_c \rangle",fontsize=15)
xlabel("t",fontsize=15)

#Open the data:
f = jldopen("../data/HPtensors.jld2","r")
T = read(f,"t1")
y1 = read(f,"y1")
close("../data/HPtensors.jld2")

plot(T,y1,"-or")

#Save to latex templates pics directory
savefig("../report_latex/Pics/HPtensors.png")
