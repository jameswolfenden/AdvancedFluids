using DelimitedFiles, Plots, LinearAlgebra

n = 101

p0 = readdlm("2d/resultsdomain/domain0p.csv",',')
p1 = readdlm("2d/resultsdomain/domain1p.csv",',')
p2 = readdlm("2d/resultsdomain/domain2p.csv",',')
p3 = readdlm("2d/resultsdomain/domain3p.csv",',')
rho0 = readdlm("2d/resultsdomain/domain0rho.csv",',')
rho1 = readdlm("2d/resultsdomain/domain1rho.csv",',')
rho2 = readdlm("2d/resultsdomain/domain2rho.csv",',')
rho3 = readdlm("2d/resultsdomain/domain3rho.csv",',')
u0 = readdlm("2d/resultsdomain/domain0u.csv",',')
u1 = readdlm("2d/resultsdomain/domain1u.csv",',')
u2 = readdlm("2d/resultsdomain/domain2u.csv",',')
u3 = readdlm("2d/resultsdomain/domain3u.csv",',')
v0 = readdlm("2d/resultsdomain/domain0v.csv",',')
v1 = readdlm("2d/resultsdomain/domain1v.csv",',')
v2 = readdlm("2d/resultsdomain/domain2v.csv",',')
v3 = readdlm("2d/resultsdomain/domain3v.csv",',')

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
    pBig = [p2i[2:end-1, 2:end-2] p0i[2:end-1, 2:end-2] p3i[2:end-1, 2:end-2]; zeros(size(p1i,1)-2, size(p2i,2)-3) p1i[2:end-1, 2:end-2] zeros(size(p1i,1)-2, size(p3i,2)-3)]
    clims = (1, 1.06)
    pPlt = heatmap(pBig, clim=clims, yflip=true)
    frame(pAnim, pPlt)
end
mp4(pAnim, "2d/resultsdomain/p.mp4", fps=30)

rhoAnim = Animation()
for i in toplay
    rho0i = rho0[i*size0+1:(i+1)*size0,:]
    rho1i = rho1[i*size1+1:(i+1)*size1,:]
    rho2i = rho2[i*size2+1:(i+1)*size2,:]
    rho3i = rho3[i*size3+1:(i+1)*size3,:]
    rhoBig = [rho2i[2:end-1, 2:end-2] rho0i[2:end-1, 2:end-2] rho3i[2:end-1, 2:end-2]; zeros(size(rho1i,1)-2, size(rho2i,2)-3) rho1i[2:end-1, 2:end-2] zeros(size(rho1i,1)-2, size(rho3i,2)-3)]
    clims = (1.3, 1.35)
    rhoPlt = heatmap(rhoBig, clim=clims, yflip=true)
    frame(rhoAnim, rhoPlt)
end
mp4(rhoAnim, "2d/resultsdomain/rho.mp4", fps=30)

uAnim = Animation()
for i in toplay
    u0i = u0[i*size0+1:(i+1)*size0,:]
    u1i = u1[i*size1+1:(i+1)*size1,:]
    u2i = u2[i*size2+1:(i+1)*size2,:]
    u3i = u3[i*size3+1:(i+1)*size3,:]
    uBig = [u2i[2:end-1, 2:end-2] u0i[2:end-1, 2:end-2] u3i[2:end-1, 2:end-2]; zeros(size(u1i,1)-2, size(u2i,2)-3) u1i[2:end-1, 2:end-2] zeros(size(u1i,1)-2, size(u3i,2)-3)]
    clims = (-0.05, 0.05)
    uPlt = heatmap(uBig, clim=clims, yflip=true)
    frame(uAnim, uPlt)
end
mp4(uAnim, "2d/resultsdomain/u.mp4", fps=30)

vAnim = Animation()
for i in toplay
    v0i = v0[i*size0+1:(i+1)*size0,:]
    v1i = v1[i*size1+1:(i+1)*size1,:]
    v2i = v2[i*size2+1:(i+1)*size2,:]
    v3i = v3[i*size3+1:(i+1)*size3,:]
    vBig = [v2i[2:end-1, 2:end-2] v0i[2:end-1, 2:end-2] v3i[2:end-1, 2:end-2]; zeros(size(v1i,1)-2, size(v2i,2)-3) v1i[2:end-1, 2:end-2] zeros(size(v1i,1)-2, size(v3i,2)-3)]
    clims = (-0.05, 0.05)
    vPlt = heatmap(vBig, clim=clims, yflip=true)
    frame(vAnim, vPlt)
end
mp4(vAnim, "2d/resultsdomain/v.mp4", fps=30)

