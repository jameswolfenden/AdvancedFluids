p0 = readdlm("3D/resultsdomain/domain0p.csv",',')
p1 = readdlm("3D/resultsdomain/domain1p.csv",',')
p2 = readdlm("3D/resultsdomain/domain2p.csv",',')
p3 = readdlm("3D/resultsdomain/domain3p.csv",',')

config = readdlm("3D/resultsdomain/iterations.csv",',')
n = Int(readdlm("3D/resultsdomain/iterations.csv",',')[1] + 1)
iterations = Int(config[1]) + 1
cellspermeterx = Int(config[2])
cellspermetery = Int(config[3])
cellspermeterz = Int(config[4])

time = readdlm(path*"elapsedtime.csv",',')


println("Read in p")

size0 = Int(size(p0,1)/n)
size1 = Int(size(p1,1)/n)
size2 = Int(size(p2,1)/n)
size3 = Int(size(p3,1)/n)

i=55


p0i = p0[i*size0+1:(i+1)*size0,:]
p1i = p1[i*size1+1:(i+1)*size1,:]
p2i = p2[i*size2+1:(i+1)*size2,:]
p3i = p3[i*size3+1:(i+1)*size3,:]
pBig = [p2i[2:end-1, 2:end-1] p0i[2:end-1, 2:end-1] p3i[2:end-1, 2:end-1]]
clims=(101.3,101.4)
x = (0:size(pBig,1)-1) / cellspermeterx
y = (0:size(pBig,2)-1)/cellspermetery

heatmap(y,x,pBig[end:-1:1,:].*1e2, clim=clims, rightmargin=20px)
ylabel!("Distance up back wall (m)")
xlabel!("Distance from back wall (m)")
title!("Pressure (kPa) at t = "* first(string(time[i]),4)*" s")

