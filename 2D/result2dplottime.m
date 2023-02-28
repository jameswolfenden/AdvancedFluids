n = numel(dir("results\*.csv"))/20;

p0 = zeros([size(readmatrix("results\0pCells0.csv")) n]);
rho0 = zeros(size(p0));
u0 = zeros(size(p0));
v0 = zeros(size(p0));
dye0 = zeros(size(p0));
p1 = zeros([size(readmatrix("results\1pCells0.csv")) n]);
rho1 = zeros(size(p1));
u1 = zeros(size(p1));
v1 = zeros(size(p1));
dye1 = zeros(size(p1));
p2 = zeros([size(readmatrix("results\2pCells0.csv")) n]);
rho2 = zeros(size(p2));
u2 = zeros(size(p2));
v2 = zeros(size(p2));
dye2 = zeros(size(p2));
p3 = zeros([size(readmatrix("results\3pCells0.csv")) n]);
rho3 = zeros(size(p3));
u3 = zeros(size(p3));
v3 = zeros(size(p3));
dye3 = zeros(size(p3));

toplay=1:n;
%toplay = 1:10:1001;

for i = toplay
    
    p0(:,:,i) = readmatrix("results\0pCells"+(i-1)+".csv");
    rho0(:,:,i) = readmatrix("results\0rhoCells"+(i-1)+".csv");
    u0(:,:,i) = readmatrix("results\0uCells"+(i-1)+".csv");
    v0(:,:,i) = readmatrix("results\0vCells"+(i-1)+".csv");
    dye0(:,:,i) = readmatrix("results\0dyeCells"+(i-1)+".csv");
    p1(:,:,i) = readmatrix("results\1pCells"+(i-1)+".csv");
    rho1(:,:,i) = readmatrix("results\1rhoCells"+(i-1)+".csv");
    u1(:,:,i) = readmatrix("results\1uCells"+(i-1)+".csv");
    v1(:,:,i) = readmatrix("results\1vCells"+(i-1)+".csv");
    dye1(:,:,i) = readmatrix("results\1dyeCells"+(i-1)+".csv");
    p2(:,:,i) = readmatrix("results\2pCells"+(i-1)+".csv");
    rho2(:,:,i) = readmatrix("results\2rhoCells"+(i-1)+".csv");
    u2(:,:,i) = readmatrix("results\2uCells"+(i-1)+".csv");
    v2(:,:,i) = readmatrix("results\2vCells"+(i-1)+".csv");
    dye2(:,:,i) = readmatrix("results\2dyeCells"+(i-1)+".csv");
    p3(:,:,i) = readmatrix("results\3pCells"+(i-1)+".csv");
    rho3(:,:,i) = readmatrix("results\3rhoCells"+(i-1)+".csv");
    u3(:,:,i) = readmatrix("results\3uCells"+(i-1)+".csv");
    v3(:,:,i) = readmatrix("results\3vCells"+(i-1)+".csv");
    dye3(:,:,i) = readmatrix("results\3dyeCells"+(i-1)+".csv");
    
end

figure
v = VideoWriter('p.mp4', 'MPEG-4');
v.Quality = 100;
open(v);
for i=toplay
    pBig = [p2(2:size(p2,1)-1,2:size(p2,2)-2,i), p0(2:size(p0,1)-1,2:size(p0,2)-2,i), p3(2:size(p3,1)-1,2:size(p3,2)-2,i); zeros(size(p1,1)-2,size(p2,2)-3), p1(2:size(p1,1)-1,2:size(p1,2)-2,i), zeros(size(p1,1)-2,size(p3,2)-3)];
    clims = [1 1.06];
    imagesc(pBig, clims);
    hold on
    fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[size(p2,1)-1.5 size(pBig,1)+0.5 size(pBig,1)+0.5 size(p2,1)-1.5], "w");
    fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p3,2)+3.5 size(pBig,2)-size(p3,2)+3.5],[size(p2,1)-1.5 size(pBig,1)+0.5 size(pBig,1)+0.5 size(p2,1)-1.5], "w");
    title("p");
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

figure
v = VideoWriter('rho.mp4', 'MPEG-4');
v.Quality = 100;
open(v);
for i=toplay
    rhoBig = [rho2(2:size(rho2,1)-1,2:size(rho2,2)-2,i), rho0(2:size(rho0,1)-1,2:size(rho0,2)-2,i), rho3(2:size(rho3,1)-1,2:size(rho3,2)-2,i); zeros(size(rho1,1)-2,size(rho2,2)-3), rho1(2:size(rho1,1)-1,2:size(rho1,2)-2,i), zeros(size(rho1,1)-2,size(rho3,2)-3)];
    clims = [1.3 1.35];
    imagesc(rhoBig, clims);
    hold on
    fill([0.5 0.5 size(rho2,2)-2.5 size(rho2,2)-2.5],[size(rho2,1)-1.5 size(rhoBig,1)+0.5 size(rhoBig,1)+0.5 size(rho2,1)-1.5], "w");
    fill([size(rhoBig,2)+0.5 size(rhoBig,2)+0.5 size(rhoBig,2)-size(rho3,2)+3.5 size(rhoBig,2)-size(rho3,2)+3.5],[size(rho2,1)-1.5 size(rhoBig,1)+0.5 size(rhoBig,1)+0.5 size(rho2,1)-1.5], "w");
    title("rho");
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

figure
v = VideoWriter('u.mp4', 'MPEG-4');
v.Quality = 100;
open(v);
for i=toplay
    uBig = [u2(2:size(u2,1)-1,2:size(u2,2)-2,i), u0(2:size(u0,1)-1,2:size(u0,2)-2,i), u3(2:size(u3,1)-1,2:size(u3,2)-2,i); zeros(size(u1,1)-2,size(u2,2)-3), u1(2:size(u1,1)-1,2:size(u1,2)-2,i), zeros(size(u1,1)-2,size(u3,2)-3)];
    clims = [-0.01 0.01];
    imagesc(uBig, clims);
    hold on
    fill([0.5 0.5 size(u2,2)-2.5 size(u2,2)-2.5],[size(u2,1)-1.5 size(uBig,1)+0.5 size(uBig,1)+0.5 size(u2,1)-1.5], "w");
    fill([size(uBig,2)+0.5 size(uBig,2)+0.5 size(uBig,2)-size(u3,2)+3.5 size(uBig,2)-size(u3,2)+3.5],[size(u2,1)-1.5 size(uBig,1)+0.5 size(uBig,1)+0.5 size(u2,1)-1.5], "w");
    title("u");
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

figure
v = VideoWriter('v.mp4', 'MPEG-4');
v.Quality = 100;
open(v);
for i=toplay
    vBig = [v2(2:size(v2,1)-1,2:size(v2,2)-2,i), v0(2:size(v0,1)-1,2:size(v0,2)-2,i), v3(2:size(v3,1)-1,2:size(v3,2)-2,i); zeros(size(v1,1)-2,size(v2,2)-3), v1(2:size(v1,1)-1,2:size(v1,2)-2,i), zeros(size(v1,1)-2,size(v3,2)-3)];
    clims = [-0.01 0.01];
    imagesc(vBig, clims);
    hold on
    fill([0.5 0.5 size(v2,2)-2.5 size(v2,2)-2.5],[size(v2,1)-1.5 size(vBig,1)+0.5 size(vBig,1)+0.5 size(v2,1)-1.5], "w");
    fill([size(vBig,2)+0.5 size(vBig,2)+0.5 size(vBig,2)-size(v3,2)+3.5 size(vBig,2)-size(v3,2)+3.5],[size(v2,1)-1.5 size(vBig,1)+0.5 size(vBig,1)+0.5 size(v2,1)-1.5], "w");
    title("v");
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

figure
v = VideoWriter('dye.mp4', 'MPEG-4');
v.Quality = 100;
open(v);
for i=toplay
    dyeBig = [dye2(2:size(dye2,1)-1,2:size(dye2,2)-2,i), dye0(2:size(dye0,1)-1,2:size(dye0,2)-2,i), dye3(2:size(dye3,1)-1,2:size(dye3,2)-2,i); zeros(size(dye1,1)-2,size(dye2,2)-3), dye1(2:size(dye1,1)-1,2:size(dye1,2)-2,i), zeros(size(dye1,1)-2,size(dye3,2)-3)];
    clims = [0 1];
    imagesc(dyeBig, clims);
    hold on
    fill([0.5 0.5 size(dye2,2)-2.5 size(dye2,2)-2.5],[size(dye2,1)-1.5 size(dyeBig,1)+0.5 size(dyeBig,1)+0.5 size(dye2,1)-1.5], "w");
    fill([size(dyeBig,2)+0.5 size(dyeBig,2)+0.5 size(dyeBig,2)-size(dye3,2)+3.5 size(dyeBig,2)-size(dye3,2)+3.5],[size(dye2,1)-1.5 size(dyeBig,1)+0.5 size(dyeBig,1)+0.5 size(dye2,1)-1.5], "w");
    title("dye");
    frame = getframe(gcf);
    writeVideo(v,frame);
end
