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
p4 = readmatrix('2d\pCells4.csv');
rho4 = readmatrix('2d\rhoCells4.csv');
u4 = readmatrix('2d\uCells4.csv');
v4 = readmatrix('2d\vCells4.csv');
% pStart = readmatrix('2d\pStartCells.csv');
% rhoStart = readmatrix('2d\rhoStartCells.csv');
% uStart = readmatrix('2d\uStartCells.csv');
% vStart = readmatrix('2d\vStartCells.csv');

clims = [0 1];
figure
subplot(2,3,5);
imagesc(p1(2:size(p1,1)-1,2:size(p1,2)-2),clims)
title("p1")
subplot(2,3,2);
imagesc(p2(2:size(p2,1)-1,2:size(p2,2)-2),clims)
title("p2")
subplot(2,3,4);
imagesc(p3(2:size(p3,1)-1,2:size(p3,2)-2),clims)
title("p3")
subplot(2,3,6);
imagesc(p4(2:size(p4,1)-1,2:size(p4,2)-2),clims)
title("p4")

clims = [-1 1];
figure
subplot(2,3,5);
imagesc(u1(2:size(u1,1)-1,2:size(u1,2)-2),clims)
title("u1")
subplot(2,3,2);
imagesc(u2(2:size(u2,1)-1,2:size(u2,2)-2),clims)
title("u2")
subplot(2,3,4);
imagesc(u3(2:size(u3,1)-1,2:size(u3,2)-2),clims)
title("u3")
subplot(2,3,6);
imagesc(u4(2:size(u4,1)-1,2:size(u4,2)-2),clims)
title("u4")

clims = [-1 1];
figure
subplot(2,3,5);
imagesc(v1(2:size(v1,1)-1,2:size(v1,2)-2),clims)
title("v1")
subplot(2,3,2);
imagesc(v2(2:size(v2,1)-1,2:size(v2,2)-2),clims)
title("v2")
subplot(2,3,4);
imagesc(v3(2:size(v3,1)-1,2:size(v3,2)-2),clims)
title("v3")
subplot(2,3,6);
imagesc(v4(2:size(v4,1)-1,2:size(v4,2)-2),clims)
title("v4")

clims = [0 1];
figure
subplot(2,3,5);
imagesc(rho1(2:size(rho1,1)-1,2:size(rho1,2)-2),clims)
title("rho1")
subplot(2,3,2);
imagesc(rho2(2:size(rho2,1)-1,2:size(rho2,2)-2),clims)
title("rho2")
subplot(2,3,4);
imagesc(rho3(2:size(rho3,1)-1,2:size(rho3,2)-2),clims)
title("rho3")
subplot(2,3,6);
imagesc(rho4(2:size(rho4,1)-1,2:size(rho4,2)-2),clims)
title("rho4")