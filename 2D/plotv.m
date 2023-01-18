i=39


v0 = readmatrix("results\0vCells"+i+".csv");

v1 = readmatrix("results\1vCells"+i+".csv");

v2 = readmatrix("results\2vCells"+i+".csv");

v3 = readmatrix("results\3vCells"+i+".csv");


vBig = [zeros(size(v1,1)-2,size(v2,2)-3), v1(2:size(v1,1)-1,2:size(v1,2)-2), zeros(size(v1,1)-2,size(v3,2)-3); v2(2:size(v2,1)-1,2:size(v2,2)-2), v0(2:size(v0,1)-1,2:size(v0,2)-2), v3(2:size(v3,1)-1,2:size(v3,2)-2)];
figure
clims = [-2 2];
%imagesc(vBig, clims);
imagesc(vBig);
hold on
fill([0.5 0.5 size(p2,2)-2.5 size(p2,2)-2.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
fill([size(pBig,2)+0.5 size(pBig,2)+0.5 size(pBig,2)-size(p2,2)+3.5 size(pBig,2)-size(p2,2)+3.5],[0.5 size(p1,1)-1.5 size(p1,1)-1.5 0.5], 'w');
title("v");
