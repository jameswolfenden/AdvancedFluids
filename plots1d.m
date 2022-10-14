rho = readmatrix('rho.csv');
u = readmatrix('u.csv');
p = readmatrix('p.csv');
t = readmatrix('t.csv');


figure
imagesc(1:size(rho,2),t,rho)
title("rho")
figure
imagesc(1:size(u,2),t,u)
title("u")
figure
imagesc(1:size(p,2),t,p)
title("p")
