using DelimitedFiles, Plots, VideoIO, LinearAlgebra

n = Int(length(readdir("2d/results"))/16)

p0 = zeros(eltype(Float64), size(readdlm("2d/results/0pCells0.csv",','))..., n)
rho0 = zeros(eltype(Float64), size(p0))
u0 = zeros(eltype(Float64), size(p0))
v0 = zeros(eltype(Float64), size(p0))
p1 = zeros(eltype(Float64), size(readdlm("2d/results/1pCells0.csv",','))..., n)
rho1 = zeros(eltype(Float64), size(p1))
u1 = zeros(eltype(Float64), size(p1))
v1 = zeros(eltype(Float64), size(p1))
p2 = zeros(eltype(Float64), size(readdlm("2d/results/2pCells0.csv",','))..., n)
rho2 = zeros(eltype(Float64), size(p2))
u2 = zeros(eltype(Float64), size(p2))
v2 = zeros(eltype(Float64), size(p2))
p3 = zeros(eltype(Float64), size(readdlm("2d/results/3pCells0.csv",','))..., n)
rho3 = zeros(eltype(Float64), size(p3))
u3 = zeros(eltype(Float64), size(p3))
v3 = zeros(eltype(Float64), size(p3))

toplay = 1:n
# toplay = 1:10:1001

for i in toplay
    println("Loading frame $(i)")
    p0[:, :, i] = readdlm("2d/results/0pCells" * string(i - 1) * ".csv",',')
    rho0[:, :, i] = readdlm("2d/results/0rhoCells" * string(i - 1) * ".csv",',')
    u0[:, :, i] = readdlm("2d/results/0uCells" * string(i - 1) * ".csv",',')
    v0[:, :, i] = readdlm("2d/results/0vCells" * string(i - 1) * ".csv",',')
    p1[:, :, i] = readdlm("2d/results/1pCells" * string(i - 1) * ".csv",',')
    rho1[:, :, i] = readdlm("2d/results/1rhoCells" * string(i - 1) * ".csv",',')
    u1[:, :, i] = readdlm("2d/results/1uCells" * string(i - 1) * ".csv",',')
    v1[:, :, i] = readdlm("2d/results/1vCells" * string(i - 1) * ".csv",',')
    p2[:, :, i] = readdlm("2d/results/2pCells" * string(i - 1) * ".csv",',')
    rho2[:, :, i] = readdlm("2d/results/2rhoCells" * string(i - 1) * ".csv",',')
    u2[:, :, i] = readdlm("2d/results/2uCells" * string(i - 1) * ".csv",',')
    v2[:, :, i] = readdlm("2d/results/2vCells" * string(i - 1) * ".csv",',')
    p3[:, :, i] = readdlm("2d/results/3pCells" * string(i - 1) * ".csv",',')
    rho3[:, :, i] = readdlm("2d/results/3rhoCells" * string(i - 1) * ".csv",',')
    u3[:, :, i] = readdlm("2d/results/3uCells" * string(i - 1) * ".csv",',')
    v3[:, :, i] = readdlm("2d/results/3vCells" * string(i - 1) * ".csv",',')

end

toplay = 1:n

pAnim = Animation()
for i in toplay
    pBig = [p2[2:end-1, 2:end-2, i] p0[2:end-1, 2:end-2, i] p3[2:end-1, 2:end-2, i]; zeros(size(p1,1)-2, size(p2,2)-3) p1[2:end-1, 2:end-2, i] zeros(size(p1,1)-2, size(p3,2)-3)]
    clims = (1, 1.06)
    pPlt = heatmap(pBig, clim=clims, yflip=true, colorbar=:none, showaxis=false, size=(1000,1000))
    frame(pAnim, pPlt)
end
mp4(pAnim, "2d/p.mp4", fps = 30)

rhoAnim = Animation()
for i in toplay
    rhoBig = [rho2[2:end-1, 2:end-2, i] rho0[2:end-1, 2:end-2, i] rho3[2:end-1, 2:end-2, i]; zeros(size(rho1,1)-2, size(rho2,2)-3) rho1[2:end-1, 2:end-2, i] zeros(size(rho1,1)-2, size(rho3,2)-3)]
    clims = (1.3, 1.35)
    rhoPlt = heatmap(rhoBig, clim=clims, yflip=true, colorbar=:none, showaxis=false, size=(1000,1000))
    frame(rhoAnim, rhoPlt)
end
mp4(rhoAnim, "2d/rho.mp4", fps = 30)

uAnim = Animation()
for i in toplay
    uBig = [u2[2:end-1, 2:end-2, i] u0[2:end-1, 2:end-2, i] u3[2:end-1, 2:end-2, i]; zeros(size(u1,1)-2, size(u2,2)-3) u1[2:end-1, 2:end-2, i] zeros(size(u1,1)-2, size(u3,2)-3)]
    clims = (-0.02, 0.02)
    uPlt = heatmap(uBig, clim=clims, yflip=true, colorbar=:none, showaxis=false, size=(1000,1000))
    frame(uAnim, uPlt)
end
mp4(uAnim, "2d/u.mp4", fps = 30)

vAnim = Animation()
for i in toplay
    vBig = [v2[2:end-1, 2:end-2, i] v0[2:end-1, 2:end-2, i] v3[2:end-1, 2:end-2, i]; zeros(size(v1,1)-2, size(v2,2)-3) v1[2:end-1, 2:end-2, i] zeros(size(v1,1)-2, size(v3,2)-3)]
    clims = (-0.02, 0.02)
    vPlt = heatmap(vBig, clim=clims, yflip=true, colorbar=:none, showaxis=false, size=(1000,1000))
    frame(vAnim, vPlt)
end
mp4(vAnim, "2d/v.mp4", fps = 30)


