i=100

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

pBig = [p2(2:size(p2,1)-1,2:size(p2,2)-2), p0(2:size(p0,1)-1,2:size(p0,2)-2), p3(2:size(p3,1)-1,2:size(p3,2)-2); zeros(size(p1,1)-2,size(p2,2)-3), p1(2:size(p1,1)-1,2:size(p1,2)-2), zeros(size(p1,1)-2,size(p3,2)-3)];
figure
clims = [0 1];
%imagesc(pBig, clims);
imagesc(pBig);
hold on
fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[size(p2,1)-1.5 size(pBig,1)+0.5 size(pBig,1)+0.5 size(p2,1)-1.5], "w");
fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p3,2)+3.5 size(pBig,2)-size(p3,2)+3.5],[size(p2,1)-1.5 size(pBig,1)+0.5 size(pBig,1)+0.5 size(p2,1)-1.5], "w");
title("p");

rhoBig = [rho2(2:size(rho2,1)-1,2:size(rho2,2)-2), rho0(2:size(rho0,1)-1,2:size(rho0,2)-2), rho3(2:size(rho3,1)-1,2:size(rho3,2)-2); zeros(size(rho1,1)-2,size(rho2,2)-3), rho1(2:size(rho1,1)-1,2:size(rho1,2)-2), zeros(size(rho1,1)-2,size(rho3,2)-3)];
figure
clims = [0 1];
%imagesc(rhoBig, clims);
imagesc(rhoBig);
hold on
fill([0.5 0.5 size(rho2,2)-2.5 size(rho2,2)-2.5],[size(rho2,1)-1.5 size(rhoBig,1)+0.5 size(rhoBig,1)+0.5 size(rho2,1)-1.5], "w");
fill([size(rhoBig,2)+0.5 size(rhoBig,2)+0.5 size(rhoBig,2)-size(rho3,2)+3.5 size(rhoBig,2)-size(rho3,2)+3.5],[size(rho2,1)-1.5 size(rhoBig,1)+0.5 size(rhoBig,1)+0.5 size(rho2,1)-1.5], "w");
title("rho");

uBig = [u2(2:size(u2,1)-1,2:size(u2,2)-2), u0(2:size(u0,1)-1,2:size(u0,2)-2), u3(2:size(u3,1)-1,2:size(u3,2)-2); zeros(size(u1,1)-2,size(u2,2)-3), u1(2:size(u1,1)-1,2:size(u1,2)-2), zeros(size(u1,1)-2,size(u3,2)-3)];
figure
clims = [-1 1];
%imagesc(uBig, clims);
imagesc(uBig);
hold on
fill([0.5 0.5 size(u2,2)-2.5 size(u2,2)-2.5],[size(u2,1)-1.5 size(uBig,1)+0.5 size(uBig,1)+0.5 size(u2,1)-1.5], "w");
fill([size(uBig,2)+0.5 size(uBig,2)+0.5 size(uBig,2)-size(u3,2)+3.5 size(uBig,2)-size(u3,2)+3.5],[size(u2,1)-1.5 size(uBig,1)+0.5 size(uBig,1)+0.5 size(u2,1)-1.5], "w");
title("u");

vBig = [v2(2:size(v2,1)-1,2:size(v2,2)-2), v0(2:size(v0,1)-1,2:size(v0,2)-2), v3(2:size(v3,1)-1,2:size(v3,2)-2); zeros(size(v1,1)-2,size(v2,2)-3), v1(2:size(v1,1)-1,2:size(v1,2)-2), zeros(size(v1,1)-2,size(v3,2)-3)];
figure
clims = [-2 2];
%imagesc(vBig, clims);
imagesc(vBig);
hold on
fill([0.5 0.5 size(v2,2)-2.5 size(v2,2)-2.5],[size(v2,1)-1.5 size(vBig,1)+0.5 size(vBig,1)+0.5 size(v2,1)-1.5], "w");
fill([size(vBig,2)+0.5 size(vBig,2)+0.5 size(vBig,2)-size(v3,2)+3.5 size(vBig,2)-size(v3,2)+3.5],[size(v2,1)-1.5 size(vBig,1)+0.5 size(vBig,1)+0.5 size(v2,1)-1.5], "w");
title("v");
