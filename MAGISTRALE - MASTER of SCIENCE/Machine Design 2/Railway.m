clear all; close all;

mtot = 17000;
mj = 15578;
b=1015;
s=746.5;
h1=1280;
R=430;
Rb1=282;
Rb2=272;
Ff=45000;
gamma=0.35;
y2=396.5;
m1=115;
m2=138;
g=9.81;
coeff=0.175;

%CALCOLO DI FORZE E MOMENTI
P1=(0.625+0.0875*h1/b)*mj*g;
P2=(0.625-0.0875*h1/b)*mj*g;
H=coeff*mj*g;

Y2=H;
Y1=2*Y2;

F1=m1*g;
F2=m2*g;

Q2=(P2*(b+s)-P1*(b-s)-H*R-F2*y2-F1*2*s)/(2*s);
Q1=P1+P2-F1*2-F2-Q2;

M1=@(x) P1*x;
M2=@(x) P1*(b-s)+Y1*R+(P1-Q1-F1)*x;
M3=@(x) P2*x;
M4=@(x) P2*(b-s)+Y2*R+(P2-Q2-F1)*x;

if (M2(y2)- M4(2*s-y2))<10^(-4)
    disp('Momenti ok');
end

%DIAGRAMMA
% figure(1);hold on;
% title('Axial load');
% xlabel('Distance (mm)');
% ylabel('Force (N)');
num=1000;
dx=(b-s)/num;

% xdisp=b-s:dx:b+s;
% plot(xdisp,Y1);
% 
% xdisp=b+s:dx:2*b;
% plot(xdisp,H);

figure(2); hold on;
title('Momentum');
xlabel('Distance (mm)');
ylabel('Momentum (N*mm)');

xdisp=0:dx:(b-s);
plot(xdisp,M1(xdisp));

xdisp=0:dx:y2;
plot(xdisp+(b-s),M2(xdisp));

xdisp=0:dx:(b-s);
plot(2*b-xdisp,M3(xdisp));

xdisp=0:dx:2*s-y2;
plot(b+s-xdisp,M4(xdisp));

line([0 2*b],[0 0],'linewidth',4);

%Braking system
FfS=15750;
A1=10590;
B1=5159.5;
A2=FfS;
B2=FfS;

% figure(3);hold on;
% xlabel('Distance (mm)');
% ylabel('Momentum (N*mm)');
% title('Momentum YZ');
Mb1=@(x) A1*x;
Mb2=@(x) B1*x;
Mb3=@(x) A2*x;
Mb4=@(x) B2*x;

% xdisp=0:dx:(b-s);
% plot(xdisp,Mb1(xdisp)+Mb3(xdisp));
% 
% xdisp=(b-s):dx:(b-s+y2);
% plot(xdisp,Mb1(xdisp)+Mb3(b-s));
% 
% xdisp=0:dx:(b-s);
% plot(2*b-xdisp,Mb2(xdisp)+Mb4(xdisp));
% 
% xdisp=(b-s):dx:(b+s-y2);
% plot(2*b-xdisp,Mb2(xdisp)+Mb4(b-s));

%
Ffx=15310;
Mb5=@(x) Ffx*x;
% 
% figure(4);hold on;
% xlabel('Distance (mm)');
% ylabel('Momentum (N*mm)');
% title('Momentum YX');
% 
% xdisp=0:dx:b-s;
% plot(xdisp,Mb5(xdisp));
% plot(2*b-xdisp,Mb5(xdisp));
% 
% xdisp=b-s:dx:b+s;
% plot(xdisp,Mb5(b-s));
% 
% figure(5);hold on;
% xlabel('Distance (mm)');
% ylabel('Momentum (N*mm)');
% title('Torque');
% 
 Mt=@(x) FfS*Rb1-Ffx*R;
% 
% xdisp=0:dx:y2;
% plot(b-s+xdisp,-Mt(xdisp));
% 
% xdisp=0:dx:2*s-y2;
% plot(y2+b-s+xdisp,Mt(xdisp));

%TOTALE
figure(6); hold on;
xlabel('Distance (mm)');
ylabel('Momentum (N*mm)');
title('Total Moment YZ');

xdisp=0:dx:b-s;
plot(xdisp,M1(xdisp)+Mb1(xdisp)+Mb3(xdisp));

xdisp=0:dx:y2;
plot(xdisp+b-s,M2(xdisp)+Mb1(xdisp+b-s)+Mb3(b-s));

xdisp=0:dx:b-s;
plot(2*b-xdisp,M3(xdisp)+Mb2(xdisp)+Mb4(xdisp));

xdisp=0:dx:2*s-y2;
plot(b+s-xdisp,M4(xdisp)+Mb2(xdisp+b-s)+Mb4(b-s));

%% most stressed section
D=180;
d=60.5;
dist=128.5;

Mx= M2(0+dist)+Mb1(b-s+dist)+Mb3(b-s+dist);
Mz= Mb5(b-s+dist);
Mtot=sqrt(Mx^2+Mz^2);

sigma=32*Mtot*D/(pi*(D^4-d^4))+Y1/(pi/4*(D^2-d^2))

tau=16*Mt(0)*D/(pi*(D^4-d^4))

vm=sqrt(sigma^2+3*tau^2)

%% Stressed section solo peso
sigma_w=32*M2(dist)*D/(pi*(D^4-d^4))+Y1/(pi/4*(D^2-d^2))

