i=39

p0 = readmatrix("results\0pCells"+i+".csv");
rho0 = readmatrix("results\0rhoCells"+i+".csv");
u0 = readmatrix("results\0uCells"+i+".csv");
v0 = readmatrix("results\0vCells"+i+".csv");
p1 = readmatrix("results\1pCells"+i+".csv");
rho1 = readmatrix("results\1rhoCells"+i+".csv");
u1 = readmatrix("results\1uCells"+i+".csv");
v1 = readmatrix("results\1vCells"+i+".csv");
p2 = readmatrix("results\2pCells"+i+".csv");
rho2 = readmatrix("results\2rhoCells"+i+".csv");
u2 = readmatrix("results\2uCells"+i+".csv");
v2 = readmatrix("results\2vCells"+i+".csv");
p3 = readmatrix("results\3pCells"+i+".csv");
rho3 = readmatrix("results\3rhoCells"+i+".csv");
u3 = readmatrix("results\3uCells"+i+".csv");
v3 = readmatrix("results\3vCells"+i+".csv");

pBig = [zeros(size(p1,1),size(p2,2)-1), p1(:,1:end-1), zeros(size(p1,1),size(p3,2)-1); p2(:,1:end-1), p0(:,1:end-1), p3(:,1:end-1)];
figure
clims = [0 1];
%imagesc(pBig, clims);
imagesc(pBig);
hold on
%fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
%fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p2,2)+3.5 size(pBig,2)-size(p2,2)+3.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
title("p");

rhoBig = [zeros(size(rho1,1),size(rho2,2)-1), rho1(:,1:end-1), zeros(size(rho1,1),size(rho3,2)-1); rho2(:,1:end-1), rho0(:,1:end-1), rho3(:,1:end-1)];
figure
clims = [0 1];
%imagesc(rhoBig, clims);
imagesc(rhoBig);
hold on
%fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
%fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p2,2)+3.5 size(pBig,2)-size(p2,2)+3.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
title("rho");

uBig = [zeros(size(u1,1),size(u2,2)-1), u1(:,1:end-1), zeros(size(u1,1),size(u3,2)-1); u2(:,1:end-1), u0(:,1:end-1), u3(:,1:end-1)];
figure
clims = [-1 1];
%imagesc(uBig, clims);
imagesc(uBig);
hold on
%fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
%fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p2,2)+3.5 size(pBig,2)-size(p2,2)+3.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
title("u");

vBig = [zeros(size(v1,1),size(v2,2)-1), v1(:,1:end-1), zeros(size(v1,1),size(v3,2)-1); v2(:,1:end-1), v0(:,1:end-1), v3(:,1:end-1)];
figure
clims = [-2 2];
%imagesc(vBig, clims);
imagesc(vBig);
hold on
%fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
%fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p2,2)+3.5 size(pBig,2)-size(p2,2)+3.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
title("v");
