gamma = 1.4;

%p = -1:0.001:2;


pL=1;
rhoL=1;
uL=0;
pR=0.1;
rhoR=0.125;
uR=0;

pL=0.4;
rhoL=1;
uL=-2;
pR=0.4;
rhoR=1;
uR=2;

p=[0.4,-0.198665]

AL = 2/((gamma+1)*rhoL);
BL = (gamma-1)/(gamma+1)*pL;
aL = sqrt(gamma*pL/rhoL);
fL1 = real((p-pL).*(AL./(p+BL)).^0.5);
fL2 = real(2.*aL./(gamma-1).*((p./pL).^((gamma-1)./(2.*gamma))-1));
fL1_ = real((AL./(BL+p)).^0.5.*(1-(p-pL)./(2.*(BL+p))));
fL2_ = real(1./(pL.*aL).*(p./pL).^(-(gamma+1)./(2.*gamma)));


AR = 2/((gamma+1)*rhoR);
BR = (gamma-1)/(gamma+1)*pR;
aR = sqrt(gamma*pR/rhoR);
fR1 = real((p-pR).*(AR./(p+BR)).^0.5);
fR2 = real(2.*aR./(gamma-1).*((p./pR).^((gamma-1)./(2.*gamma))-1));
fR1_ = real((AR./(BR+p)).^0.5.*(1-(p-pR)./(2.*(BR+p))));
fR2_ = real(1./(pR.*aR).*(p./pR).^(-(gamma+1)./(2.*gamma)));




% figure
% hold on
% plot(p,fL1+fR1+uR-uL);
% plot(p,fL1+fR2+uR-uL);
% plot(p,fL2+fR1+uR-uL);
% plot(p,fL2+fR2+uR-uL);
% yline(0);
% axis([-1 2 -15 15]);
% 
% figure
% hold on
% plot(p,fL1_+fR1_);
% plot(p,fL1_+fR2_);
% plot(p,fL2_+fR1_);
% plot(p,fL2_+fR2_);
% axis([-1 2 -15 15]);
% 
% figure
% plot(p,[fL1_; fL2_; fR1_; fR2_])


