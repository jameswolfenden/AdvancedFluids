p0 = readmatrix('2d\pCells0.csv');
rho0 = readmatrix('2d\rhoCells0.csv');
u0 = readmatrix('2d\uCells0.csv');
v0 = readmatrix('2d\vCells0.csv');
p1 = readmatrix('2d\pCells1.csv');
rho1 = readmatrix('2d\rhoCells1.csv');
u1 = readmatrix('2d\uCells1.csv');
v1 = readmatrix('2d\vCells1.csv');
p2 = readmatrix('2d\pCells2.csv');
rho2 = readmatrix('2d\rhoCells2.csv');
u2 = readmatrix('2d\uCells2.csv');
v2 = readmatrix('2d\vCells2.csv');
p3 = readmatrix('2d\pCells3.csv');
rho3 = readmatrix('2d\rhoCells3.csv');
u3 = readmatrix('2d\uCells3.csv');
v3 = readmatrix('2d\vCells3.csv');
% pStart = readmatrix('2d\pStartCells.csv');
% rhoStart = readmatrix('2d\rhoStartCells.csv');
% uStart = readmatrix('2d\uStartCells.csv');
% vStart = readmatrix('2d\vStartCells.csv');

pBig = [zeros(size(p1,1)-2,size(p2,2)-3), p1(2:size(p1,1)-1,2:size(p1,2)-2), zeros(size(p1,1)-2,size(p3,2)-3); p2(2:size(p2,1)-1,2:size(p2,2)-2), p0(2:size(p0,1)-1,2:size(p0,2)-2), p3(2:size(p3,1)-1,2:size(p3,2)-2)];
figure
clims = [0 1];
imagesc(pBig, clims);
hold on
fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p2,2)+3.5 size(pBig,2)-size(p2,2)+3.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
title("p");

rhoBig = [zeros(size(rho1,1)-2,size(rho2,2)-3), rho1(2:size(rho1,1)-1,2:size(rho1,2)-2), zeros(size(rho1,1)-2,size(rho3,2)-3); rho2(2:size(rho2,1)-1,2:size(rho2,2)-2), rho0(2:size(rho0,1)-1,2:size(rho0,2)-2), rho3(2:size(rho3,1)-1,2:size(rho3,2)-2)];
figure
clims = [0 1];
imagesc(rhoBig, clims);
hold on
fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p2,2)+3.5 size(pBig,2)-size(p2,2)+3.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
title("rho");

uBig = [zeros(size(u1,1)-2,size(u2,2)-3), u1(2:size(u1,1)-1,2:size(u1,2)-2), zeros(size(u1,1)-2,size(u3,2)-3); u2(2:size(u2,1)-1,2:size(u2,2)-2), u0(2:size(u0,1)-1,2:size(u0,2)-2), u3(2:size(u3,1)-1,2:size(u3,2)-2)];
figure
clims = [-1 1];
imagesc(uBig, clims);
hold on
fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p2,2)+3.5 size(pBig,2)-size(p2,2)+3.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
title("u");

vBig = [zeros(size(v1,1)-2,size(v2,2)-3), v1(2:size(v1,1)-1,2:size(v1,2)-2), zeros(size(v1,1)-2,size(v3,2)-3); v2(2:size(v2,1)-1,2:size(v2,2)-2), v0(2:size(v0,1)-1,2:size(v0,2)-2), v3(2:size(v3,1)-1,2:size(v3,2)-2)];
figure
clims = [-1 1];
imagesc(vBig, clims);
hold on
fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p2,2)+3.5 size(pBig,2)-size(p2,2)+3.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
title("v");



% clims = [0 1];
% figure
% subplot(2,3,5);
% imagesc(p0(2:size(p0,1)-1,2:size(p0,2)-2),clims)
% title("p0")
% subplot(2,3,2);
% imagesc(p1(2:size(p1,1)-1,2:size(p1,2)-2),clims)
% title("p1")
% subplot(2,3,4);
% imagesc(p2(2:size(p2,1)-1,2:size(p2,2)-2),clims)
% title("p2")
% subplot(2,3,6);
% imagesc(p3(2:size(p3,1)-1,2:size(p3,2)-2),clims)
% title("p3")
% 
% clims = [-1 1];
% figure
% subplot(2,3,5);
% imagesc(u0(2:size(u0,1)-1,2:size(u0,2)-2),clims)
% title("u0")
% subplot(2,3,2);
% imagesc(u1(2:size(u1,1)-1,2:size(u1,2)-2),clims)
% title("u1")
% subplot(2,3,4);
% imagesc(u2(2:size(u2,1)-1,2:size(u2,2)-2),clims)
% title("u2")
% subplot(2,3,6);
% imagesc(u3(2:size(u3,1)-1,2:size(u3,2)-2),clims)
% title("u3")
% 
% clims = [-1 1];
% figure
% subplot(2,3,5);
% imagesc(v0(2:size(v0,1)-1,2:size(v0,2)-2),clims)
% title("v0")
% subplot(2,3,2);
% imagesc(v1(2:size(v1,1)-1,2:size(v1,2)-2),clims)
% title("v1")
% subplot(2,3,4);
% imagesc(v2(2:size(v2,1)-1,2:size(v2,2)-2),clims)
% title("v2")
% subplot(2,3,6);
% imagesc(v3(2:size(v3,1)-1,2:size(v3,2)-2),clims)
% title("v3")
% 
% clims = [0 1];
% figure
% subplot(2,3,5);
% imagesc(rho0(2:size(rho0,1)-1,2:size(rho0,2)-2),clims)
% title("rho0")
% subplot(2,3,2);
% imagesc(rho1(2:size(rho1,1)-1,2:size(rho1,2)-2),clims)
% title("rho1")
% subplot(2,3,4);
% imagesc(rho2(2:size(rho2,1)-1,2:size(rho2,2)-2),clims)
% title("rho2")
% subplot(2,3,6);
% imagesc(rho3(2:size(rho3,1)-1,2:size(rho3,2)-2),clims)
% title("rho3")