points = readmatrix('graph1.csv');
pStar = points(1);
uStar = points(2);
rhoStarL = points(3);
rhoStarR = points(4);
uShock = points(5);
uHead = points(6);
uTail = points(7);

xtRange = [-2 2 0 1];

pL = points(8);
rhoL = points(9);
uL = points(10);
pR = points(11);
rhoR = points(12);
uR = points(13);

fan = readmatrix('graph2.csv');
xtFan = fan(:,1);
rhoFan = fan(:,2);
uFan = fan(:,3);
pFan = fan(:,4);

figure
hold on;
title("rho");
axis(xtRange);
plot(xtFan,rhoFan);
plot([-2 min(xtFan)],[rhoL rhoFan(xtFan==min(xtFan))]);
plot([max(xtFan) uStar uStar],[rhoStarL rhoStarL rhoStarR]);
plot([uStar uShock uShock],[rhoStarR rhoStarR rhoR]);
plot([uShock 2],[rhoR rhoR]);

figure
hold on;
title("u");
axis(xtRange);
plot(xtFan,uFan);
plot([-2 min(xtFan)],[uL uFan(xtFan==min(xtFan))]);
plot([max(xtFan) uShock uShock],[uStar uStar uR]);
plot([uShock 2],[uR uR]);

figure
hold on;
title("p");
axis(xtRange);
plot(xtFan,pFan);
plot([-2 min(xtFan)],[pL pFan(xtFan==min(xtFan))]);
plot([max(xtFan) uShock uShock],[pStar pStar pR]);
plot([uShock 2],[pR pR]);

