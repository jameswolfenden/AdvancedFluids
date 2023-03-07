using Plots, Plots.PlotMeasures

meshsize = [1000,100,200,300,400,34,67]
maxforce = [2.23,1.7678,2.02,2.12,2.17,0.905709343,1.513477389]

permsortmeshsize = sortperm(meshsize)
sortmeshsize = meshsize[permsortmeshsize]
sortmaxforce = maxforce[permsortmeshsize]

meshplot = sortmeshsize.^3

# plot with meshplot as a log scale
plot(meshplot, sortmaxforce, xlabel="Mesh Density (cells/m³)", ylabel="Max Force on Back Wall (N)", legend=:none, xticks=[10^5,10^6,10^7,10^8,10^9], yticks=0:0.25:2.5, title="Mesh Convergence", marker=:circle, markersize=4, markerstrokewidth=2, linewidth=2, color=:black, xscale=:log10, minorgrid=:true)