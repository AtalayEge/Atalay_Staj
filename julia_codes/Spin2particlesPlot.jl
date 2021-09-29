#Plotting script

using PyPlot
using JLD2

#Initialize the figure:
figure(figsize=(8,5))
title("Magnetization of a Single Particle")
ylabel(L"\langle S^z_1 \rangle",fontsize=15)
xlabel("t",fontsize=15)

#Open the data:
f = jldopen("C:/Users/asus/Desktop/staj2/staj_dosya/julia_codes/magnedata2.jld2","r")
T = read(f,"T")
y = read(f,"y")
y2 = read(f,"y2")
close("C:/Users/asus/Desktop/staj2/staj_dosya/julia_codes/magnedata2.jld2")

plot(T,y2,"r",T,y,"--k")
legend([L"tr(\rho S_1^z)",L"\langle\psi|S_1^z|\psi\rangle"],fontsize=15)

#Save to latex templates pics directory
savefig("C:/Users/asus/Desktop/staj2/staj_dosya/report_latex/Pics/spinmagne.png")
