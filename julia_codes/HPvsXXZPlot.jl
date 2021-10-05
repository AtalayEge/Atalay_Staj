#Plotting script

using PyPlot
using JLD2

#Initialize the figure:
figure(figsize=(8,5))
title("Magnetization Comparison for XXZ and H-P Models")
ylabel(L"\langle S^z_1 \rangle",fontsize=15)
xlabel("t",fontsize=15)

#Open the data:
f = jldopen("../data/XXZmagne.jld2","r")
f2 = jldopen("../data/HPmagne.jld2","r")
T = read(f,"T")
y1 = read(f,"y1")
y2 = read(f2,"y1")
close("../data/XXZmagne.jld2")
close("../data/HPmagne.jld2")

plot(T,y1,"-or",T,y2,"--k")
legend(["XXZ Spin","H-P Bosons"])
#Save to latex templates pics directory
savefig("../report_latex/Pics/XXZvsHP.png")
