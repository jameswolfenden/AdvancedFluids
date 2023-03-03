

p = readdlm("3D/resultsdomain/wallPressures.csv",',')
iterations = Int(readdlm("3D/resultsdomain/iterations.csv",',')[1] + 1)
p = reshape(p, (Int(size(p,1)/iterations), iterations, size(p,2)))


toplay = 0:n-2
pWallAnim = Animation()
for i in toplay
    pWallPlt = heatmap(p[:,i+1,:], clim=(1.01325,1.017))
    frame(pWallAnim, pWallPlt)
end
mp4(pWallAnim, "3D/resultsdomain/pWall.mp4", fps=30)

ATM = 1.01325
pAbs = p .- ATM

cellspermeterz = 500
cellspermetery = 200
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
y = totalforce .* x .* (H^2 .- x.^2).^(3 / 2) / (9 * sqrt(3) * H * E * I)

