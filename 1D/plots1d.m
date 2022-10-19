rho = readmatrix('rho.csv');
u = readmatrix('u.csv');
p = readmatrix('p.csv');
t = readmatrix('t.csv');

rho = flip(rho,1);
u = flip(u,1);
p = flip(p,1);




figure
imagesc(1:size(rho,2),t,rho)
title("rho")
colormap(gray)
figure
imagesc(1:size(u,2),t,u)
title("u")
colormap(gray)
figure
imagesc(1:size(p,2),t,p)
title("p")
colormap(gray)
