using DelimitedFiles, Plots, Plots.PlotMeasures
path = "3D/resultsdomain/"

p0 = readdlm(path*"domain0p.csv",',')
p1 = readdlm(path*"domain1p.csv",',')
p2 = readdlm(path*"domain2p.csv",',')
p3 = readdlm(path*"domain3p.csv",',')

config = readdlm(path*"iterations.csv",',')
iterations = Int(config[1]) + 1
n=iterations
cellspermeterx = Int(config[2])
cellspermetery = Int(config[3])
cellspermeterz = Int(config[4])

time = readdlm(path*"elapsedtime.csv",',')
time = time.*1000


println("Read in p")

size0 = Int(size(p0,1)/n)
size1 = Int(size(p1,1)/n)
size2 = Int(size(p2,1)/n)
size3 = Int(size(p3,1)/n)

i=120


p0i = p0[i*size0+1:(i+1)*size0,:]
p1i = p1[i*size1+1:(i+1)*size1,:]
p2i = p2[i*size2+1:(i+1)*size2,:]
p3i = p3[i*size3+1:(i+1)*size3,:]
pBig = [p2i[2:end-1, 2:end-1] p0i[2:end-1, 2:end-1] p3i[2:end-1, 2:end-1]]
clims=(101.3,101.525)
x = (0:size(pBig,1)-1) / cellspermeterx*100
y = (0:size(pBig,2)-1)/cellspermetery*100

heatmap(y,x,pBig[end:-1:1,:]./1e3, clim=clims, rightmargin=20px)
ylabel!("Distance up back wall (cm)")
xlabel!("Distance from back wall (cm)")
title!("Pressure (kPa) at t = "* first(string(time[i]),4)*" ms")

