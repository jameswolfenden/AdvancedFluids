using DelimitedFiles, Plots

path = "3D/resultsdomain/"

p = readdlm(path*"wallPressures.csv",',')
config = readdlm(path*"iterations.csv",',')
iterations = Int(config[1]) + 1
cellspermetery = Int(config[3])
cellspermeterz = Int(config[4])
p = reshape(p, (Int(size(p,1)/iterations), iterations, size(p,2)))

time = readdlm(path*"elapsedtime.csv",',')
time = time.*1000


toplay = 0:iterations-2
pWallAnim = Animation()
for i in toplay
    pWallPlt = heatmap(p[:,i+1,:] ./1000, clim=(101.30,101.50))
    frame(pWallAnim, pWallPlt)
end
mp4(pWallAnim, path*"pWall.mp4", fps=30)

ATM = 101325
pAbs = p .- ATM

cellarea = 1/cellspermetery*1/cellspermeterz

forcecell = pAbs .* cellarea

totalforce = dropdims(sum(sum(forcecell, dims=1),dims=3), dims=(1,3))

maxlocation = dropdims(findmax(p, dims=(1,3))[2], dims=(1,3))

maxrow = dropdims(findmax(sum(forcecell, dims=3), dims=1)[2], dims=(1,3))
maxrowarray = zeros(size(p,2))
for i in 1:size(p,2)
    maxrowarray[i] = Int(maxrow[i][1])
end

H = 2.0
E = 2.0e11
t = 0.25e-3
W = 1.0
I = t^3 * W / 12

# find the centre of pressure at each time step
rowpressure = sum(forcecell, dims=3)
# make rowpressures below zero zero
rowpressure[rowpressure .< 0] .= 0

# array of the physical x positions of each row
xs = ((size(p,1):-1:1).-0.5)./cellspermetery

sxsf = dropdims(sum(rowpressure.*xs, dims=1), dims=(1,3))

sf = dropdims(sum(rowpressure, dims=1), dims=(1,3))

centrepressure = sxsf./ sf
centrepressure[sf.==0].=0
centrepressure = abs.(centrepressure)
# ensure centre pressure is between 0 and H
centrepressure[centrepressure .< 0] .= 0
centrepressure[centrepressure .> H] .= H
centrepressure = H .- centrepressure

maxpressure = maxrowarray./cellspermetery
maxpressure[maxpressure .== 0] .= H

y = zeros(size(centrepressure))
# equation only valid for centre pressure < H/2
y[centrepressure .< H/2] = totalforce[centrepressure .< H/2] .* centrepressure[centrepressure .< H/2] .* (H^2 .- centrepressure[centrepressure .< H/2].^2).^(3 / 2) / (9 * sqrt(3) * H * E * I)

# equation only valid for centre pressure > H/2
y[centrepressure .> H/2] = totalforce[centrepressure .> H/2] .* (H .- centrepressure[centrepressure .> H/2]) .* (H^2 .- (H .- centrepressure[centrepressure .> H/2]).^2).^(3 / 2) / (9 * sqrt(3) * H * E * I)

#y[maxpressure .> H/2] = totalforce[maxpressure .> H/2] .* (H .- maxpressure[maxpressure .> H/2]) .* (H^2 .- (H .- maxpressure[maxpressure .> H/2]).^2).^(3 / 2) / (9 * sqrt(3) * H * E * I)

# approximate as triangle
ylocation = sqrt.((H.^2 .- (H.-centrepressure).^2)/3)
a = sqrt.(ylocation.^2 .+ y.^2)
b = sqrt.((H .- ylocation).^2 .+ y.^2)
delta = a .+ b .- H
elongation = delta ./ H
# get the sign of y
signy = ones(size(y))
signy[y .!=0] = y[y .!=0] ./ abs.(y[y .!=0])
elongation = elongation .* signy

plotnum = 1
if plotnum == 1
    plot(time, totalforce, label=:none, xlims=(0,1.5), yticks=-2:0.5:3, linewidth=2, color=:black, xticks=0:0.25:2, right_margin=10px, minorgrid=:true)
    xlabel!("Time (ms)")
    ylabel!("Force (N)")
    title!("Total Force on Back Wall")
    # line at y=0
    plot!([0], linewidth=0.5, color=:black, seriestype=:hline, label=:none)
elseif plotnum == 2
    i = 60
    centrepressure[centrepressure .== H] .= NaN
    plot(time[1:i], (H.-centrepressure[1:i]) .*100, label=:none, ylims=(-0.5,25), xticks=0:0.2:1, xlims=(0,1), right_margin=10px, linewidth=2, color=:black, minorgrid=:true)
    xlabel!("Time (ms)")
    ylabel!("Height up Back Wall (cm)")
    title!("Position of Centre of Pressure on Back Wall")
elseif plotnum ==3
    i = 120
    zs = (-size(p,3)/2+0.5:size(p,3)/2-0.5)./cellspermeterz
    ys = (0:size(p,1)-1)./cellspermetery
    heatmap(zs.*100,ys.*100,p[end:-1:1,i,:]/1000, clim=(101.30,101.525), right_margin=20px)
    xlabel!("Distance along bottom of Back Wall (cm)")
    ylabel!("Distance up Back Wall (cm)")
    title!("Pressure on Back Wall (kPa) at t = "* first(string(time[i]),4)*" ms")
elseif plotnum == 4
    i = 60
    # plot elongation against time
    plot(time[1:i], elongation[1:i]*100, label=:none, ylims=(-0.1,4), xticks=0:0.2:1, xlims=(0,1), right_margin=10px, linewidth=2, color=:black, yticks = 0:1:10, minorgrid=:true)
    xlabel!("Time (ms)")
    ylabel!("Strain (%)")
    title!("Strain in Gold Leaf")
end
