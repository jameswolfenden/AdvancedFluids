using DelimitedFiles, Plots

p = readdlm("3D/resultsdomain/wallPressures.csv",',')
config = readdlm("3D/resultsdomain/iterations.csv",',')
iterations = Int(config[1]) + 1
cellspermetery = Int(config[3])
cellspermeterz = Int(config[4])
p = reshape(p, (Int(size(p,1)/iterations), iterations, size(p,2)))


toplay = 0:iterations-2
pWallAnim = Animation()
for i in toplay
    pWallPlt = heatmap(p[:,i+1,:], clim=(1.013,1.015))
    frame(pWallAnim, pWallPlt)
end
mp4(pWallAnim, "3D/resultsdomain/pWall.mp4", fps=30)

ATM = 1.01325
pAbs = p .- ATM

cellarea = 1/cellspermetery*1/cellspermeterz

forcecell = pAbs*1e5 .* cellarea

totalforce = dropdims(sum(sum(forcecell, dims=1),dims=3), dims=(1,3))

maxlocation = dropdims(findmax(p, dims=(1,3))[2], dims=(1,3))

maxrow = dropdims(findmax(sum(forcecell, dims=3), dims=1)[2], dims=(1,3))
maxrowarray = zeros(size(p,2))
for i in 1:size(p,2)
    maxrowarray[i] = Int(maxrow[i][1])
end

# replace any 1s with 200 in maxrowarray
maxrowarray[maxrowarray .== 1] .= 200

x = (size(p,1).-maxrowarray.+0.5)./cellspermetery

H = 2.0
E = 2.0e11
t = 0.25e-3
W = 1.0
I = t^3 * W / 12

# find the centre of pressure at each time step
rowforce = sum(forcecell, dims=3)

# array of the physical x positions of each row
xs = size(p,1)/cellspermetery.-((1:size(p,1)).-0.5)./cellspermetery

sxsf = dropdims(sum(rowforce.*xs, dims=1), dims=(1,3))

sf = dropdims(sum(rowforce, dims=1), dims=(1,3))

centrepressure = sxsf./ sf
centrepressure[sf.==0].=0

y = totalforce .* centrepressure .* (H^2 .- centrepressure.^2).^(3 / 2) / (9 * sqrt(3) * H * E * I)
