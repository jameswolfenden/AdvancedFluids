p = readmatrix('2d\pEndCells.csv');
rho = readmatrix('2d\rhoEndCells.csv');
u = readmatrix('2d\uEndCells.csv');
v = readmatrix('2d\vEndCells.csv');
pStart = readmatrix('2d\pStartCells.csv');
rhoStart = readmatrix('2d\rhoStartCells.csv');
uStart = readmatrix('2d\uStartCells.csv');
vStart = readmatrix('2d\vStartCells.csv');

figure
imagesc(p(2:size(p,1)-1,2:size(p,2)-2))
title("p")
figure
imagesc(rho(2:size(rho,1)-1,2:size(rho,2)-2))
title("rho")
figure
imagesc(u(2:size(u,1)-1,2:size(u,2)-2))
title("u")
figure
imagesc(v(2:size(v,1)-1,2:size(v,2)-2))
title("v")
