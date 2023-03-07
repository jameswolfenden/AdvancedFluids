using DelimitedFiles, Plots, LinearAlgebra

# read in the number of iterations
n = Int(readdlm("3D/resultsdomain/iterations.csv",',')[1] + 1)
# print the number of iterations(+1)
println("Number of iterations(+1): $(n)")

# read in the data
p0 = readdlm("3D/resultsdomain/domain0p.csv",',')
p1 = readdlm("3D/resultsdomain/domain1p.csv",',')
p2 = readdlm("3D/resultsdomain/domain2p.csv",',')
p3 = readdlm("3D/resultsdomain/domain3p.csv",',')

println("Read in p")

# find the highest pressure in the second column of p2
maxpressure_ = maximum(p2[:,2])
println("Max pressure: $(maxpressure_)")

toplay = 0:n-2

size0 = Int(size(p0,1)/n)
size1 = Int(size(p1,1)/n)
size2 = Int(size(p2,1)/n)
size3 = Int(size(p3,1)/n)

pAnim = Animation()
for i in toplay
    p0i = p0[i*size0+1:(i+1)*size0,:]
    p1i = p1[i*size1+1:(i+1)*size1,:]
    p2i = p2[i*size2+1:(i+1)*size2,:]
    p3i = p3[i*size3+1:(i+1)*size3,:]
    pBig = [p2i[2:end-1, 2:end-1] p0i[2:end-1, 2:end-1] p3i[2:end-1, 2:end-1]; zeros(size(p1i,1)-2, size(p2i,2)-2) p1i[2:end-1, 2:end-1] zeros(size(p1i,1)-2, size(p3i,2)-2)]
    clims = (101300, 101525)
    pPlt = heatmap(pBig, clim=clims, yflip=true)
    frame(pAnim, pPlt)
end
mp4(pAnim, "3D/resultsdomain/p.mp4", fps=30)

p0 = Nothing
p1 = Nothing
p2 = Nothing
p3 = Nothing

rho0 = readdlm("3D/resultsdomain/domain0rho.csv",',')
rho1 = readdlm("3D/resultsdomain/domain1rho.csv",',')
rho2 = readdlm("3D/resultsdomain/domain2rho.csv",',')
rho3 = readdlm("3D/resultsdomain/domain3rho.csv",',')

println("Read in rho")

rhoAnim = Animation()
for i in toplay
    rho0i = rho0[i*size0+1:(i+1)*size0,:]
    rho1i = rho1[i*size1+1:(i+1)*size1,:]
    rho2i = rho2[i*size2+1:(i+1)*size2,:]
    rho3i = rho3[i*size3+1:(i+1)*size3,:]
    rhoBig = [rho2i[2:end-1, 2:end-1] rho0i[2:end-1, 2:end-1] rho3i[2:end-1, 2:end-1]; zeros(size(rho1i,1)-2, size(rho2i,2)-2) rho1i[2:end-1, 2:end-1] zeros(size(rho1i,1)-2, size(rho3i,2)-2)]
    clims = (1.267, 1.269)
    rhoPlt = heatmap(rhoBig, clim=clims, yflip=true)
    frame(rhoAnim, rhoPlt)
end
mp4(rhoAnim, "3D/resultsdomain/rho.mp4", fps=30)

rho0 = Nothing
rho1 = Nothing
rho2 = Nothing
rho3 = Nothing

u0 = readdlm("3D/resultsdomain/domain0u.csv",',')
u1 = readdlm("3D/resultsdomain/domain1u.csv",',')
u2 = readdlm("3D/resultsdomain/domain2u.csv",',')
u3 = readdlm("3D/resultsdomain/domain3u.csv",',')

println("Read in u")

uAnim = Animation()
for i in toplay
    u0i = u0[i*size0+1:(i+1)*size0,:]
    u1i = u1[i*size1+1:(i+1)*size1,:]
    u2i = u2[i*size2+1:(i+1)*size2,:]
    u3i = u3[i*size3+1:(i+1)*size3,:]
    uBig = [u2i[2:end-1, 2:end-1] u0i[2:end-1, 2:end-1] u3i[2:end-1, 2:end-1]; zeros(size(u1i,1)-2, size(u2i,2)-2) u1i[2:end-1, 2:end-1] zeros(size(u1i,1)-2, size(u3i,2)-2)]
    clims = (-1, 1)
    uPlt = heatmap(uBig, clim=clims, yflip=true)
    frame(uAnim, uPlt)
end
mp4(uAnim, "3D/resultsdomain/u.mp4", fps=30)

u0 = Nothing
u1 = Nothing
u2 = Nothing
u3 = Nothing

v0 = readdlm("3D/resultsdomain/domain0v.csv",',')
v1 = readdlm("3D/resultsdomain/domain1v.csv",',')
v2 = readdlm("3D/resultsdomain/domain2v.csv",',')
v3 = readdlm("3D/resultsdomain/domain3v.csv",',')

println("Read in v")

vAnim = Animation()
for i in toplay
    v0i = v0[i*size0+1:(i+1)*size0,:]
    v1i = v1[i*size1+1:(i+1)*size1,:]
    v2i = v2[i*size2+1:(i+1)*size2,:]
    v3i = v3[i*size3+1:(i+1)*size3,:]
    vBig = [v2i[2:end-1, 2:end-1] v0i[2:end-1, 2:end-1] v3i[2:end-1, 2:end-1]; zeros(size(v1i,1)-2, size(v2i,2)-2) v1i[2:end-1, 2:end-1] zeros(size(v1i,1)-2, size(v3i,2)-2)]
    clims = (-1, 1)
    vPlt = heatmap(vBig, clim=clims, yflip=true)
    frame(vAnim, vPlt)
end
mp4(vAnim, "3D/resultsdomain/v.mp4", fps=30)

v0 = Nothing
v1 = Nothing
v2 = Nothing
v3 = Nothing

w0 = readdlm("3D/resultsdomain/domain0w.csv",',')
w1 = readdlm("3D/resultsdomain/domain1w.csv",',')
w2 = readdlm("3D/resultsdomain/domain2w.csv",',')
w3 = readdlm("3D/resultsdomain/domain3w.csv",',')

println("Read in w")

wAnim = Animation()
for i in toplay
    w0i = w0[i*size0+1:(i+1)*size0,:]
    w1i = w1[i*size1+1:(i+1)*size1,:]
    w2i = w2[i*size2+1:(i+1)*size2,:]
    w3i = w3[i*size3+1:(i+1)*size3,:]
    wBig = [w2i[2:end-1, 2:end-1] w0i[2:end-1, 2:end-1] w3i[2:end-1, 2:end-1]; zeros(size(w1i,1)-2, size(w2i,2)-2) w1i[2:end-1, 2:end-1] zeros(size(w1i,1)-2, size(w3i,2)-2)]
    clims = (-1, 1)
    wPlt = heatmap(wBig, clim=clims, yflip=true)
    frame(wAnim, wPlt)
end
mp4(wAnim, "3D/resultsdomain/w.mp4", fps=30)

w0 = Nothing
w1 = Nothing
w2 = Nothing
w3 = Nothing
