#Plotting script

using PyPlot
using JLD2

#Initialize the figure:
figure(figsize=(8,5))
title("Expectation Value Change for N=100 H-P Chain")
ylabel(L"\langle S_c^z \rangle",fontsize=15)
xlabel("t",fontsize=15)

#Open the data:
f = jldopen("../data/HPtensors.jld2","r")
T = read(f,"t1")
y1 = read(f,"y1")
y2 = read(f,"y2")
y3 = read(f,"y3")
close("../data/HPtensors.jld2")

plot(T,y1,"-or",T,y2,"ob",T,y3,"-*k")
legend([L"\langle S_c^z \rangle",L"\langle S_{c+1}^z \rangle",
        L"\langle S_{c-1}^z \rangle"])

#Save to latex templates pics directory
savefig("../report_latex/Pics/HPtensors.png")
