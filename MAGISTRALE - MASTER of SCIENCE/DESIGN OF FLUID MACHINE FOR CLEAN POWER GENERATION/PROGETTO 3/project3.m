%% Kaplan Design
clear;clc;
Hm=12; %m
Q=10; %m^3/s
g=9.81; %m/s^2
p=8; %poles couple number
rpm=60*50/p;
omega=rpm*2*pi/60; %rad/s
omega_s=omega*sqrt(Q)/((g*Hm)^(3/4));
D_s=1.9; %from graph
eta_stat=0.87;
D=D_s*sqrt(Q)/(g*Hm)^(1/4);
P_kW=g*Hm*Q*1000*eta_stat/1000; % kW
P_HP=g*Hm*Q*1000*eta_stat/735.5;%HP
leu=g*Hm*eta_stat;
ns=rpm*sqrt(P_HP)/(Hm^(5/4));
n1=177; %from graph
Dt=D; %[m]
BonDt=0.395; %from graph
B=BonDt*Dt; %[m] wicket gate blade heigth
DronDt=0.42; %from graph
Dr=DronDt*Dt; %[m] % Root diameter
Dm=(Dr+Dt)/2;%[m] %mean diameter
Dh=Dr; %diameter of root
hblade=(Dt-Dr)/2; %blade heigth [m]
load NACA0012.txt;load GOE303.txt;load GOE517.txt;load GOE210.txt;
wg=NACA0012;
hub=GOE303;
mid=GOE517;
tip=GOE210;
alfa=-7:0.25:15;
hubcl=interp1(GOE303(:,1),GOE303(:,2),alfa);hubcd=interp1(GOE303(:,1),GOE303(:,3),alfa);
midcl=interp1(GOE517(:,1),GOE517(:,2),alfa);midcd=interp1(GOE517(:,1),GOE517(:,3),alfa);
tipcl=interp1(GOE210(:,1),GOE210(:,2),alfa);tipcd=interp1(GOE210(:,1),GOE210(:,3),alfa);

%% thomann correlation
Db=Dt;
Dold=0;
while abs(Db-Dold)>1e-5
    Dold=Db;
    Zd=sqrt(Db*1000)/4+7;
    Rb=(125/(80+Zd))*Dt/2;
    Db=Rb*2;
end
Zd=round(Zd);Rb=(125/(80+Zd))*Dt/2;Db=Rb*2; %Zd is the number of blades of wicket gate
c=2*pi*Rb/(0.72*Zd); %wicket blade chord 
%% velocity triangles
f=0.9; %relaxing factor for convergence issues
v1t=leu/(omega*Dm/2); %m/s
v2m=Q/((Dt^2-Dh^2)*pi/4); %the component of mass flow movement
v1=sqrt(v1t^2+v2m^2); %m/s
vbt=(Dm*v1t/2)/Rb; %m/s
vbm=Q/(Db*pi*B); %m/s conserved along the wicket blade
vb=sqrt(vbm^2+vbt^2); %m/s
%k=0.07;
%vb=v1/(1-k)^(1/2);
%vbm=sqrt(vb^2-vbt^2);
stag=atan(vbm/vbt);stagD=stag*180/pi; % in degree
vam=Q/(pi*(Db+2*sin(stag)*c)*B); %keeping also the chord part to increase diameter length
vat=vbt*0.5; % [deg]as first value, to start the evaluation of a real one
vold=0;
inc=[];
inc(1)=5;
n=0;
while abs(vat-vold)>1e-4
    n=n+1;
    vold=vat;
    vinft=vat+vbt;
    vinfm=vam+vbm;
    vinf=sqrt(vinft^2/4+vinfm^2/4);
    %a=atan(vinfm/vinft);
    %inc(n)=(a-stag)*180/pi;
    Ft=Q*1000*(vbt-vat);
    Clnew(n)=(Ft/(1000*vinf^2*c*B*Zd/2)+interp1(wg(:,1),wg(:,3),inc(n))*cos(inc(n)*pi/180+stag))/sin(inc(n)*pi/180+stag);
    inc(n+1)=(interp1(wg(:,2),wg(:,1),Clnew(n)));
    a=stag+inc(n+1)*pi/180;
    %vinft=vinfm/tan(a);
    %vinf=sqrt(vinft^2/4+vinfm^2/4);
    Faero=1000*vinf^2*c*B*Zd*(-interp1(wg(:,1),wg(:,3),inc(n+1))*cos(a)+Clnew(n)*sin(a))/2;
    vat=vbt-Faero/(1000*Q);
    vat=vat*f+(1-f)*vold;
end
%% Weining correction for cascade
s=2*pi*(Rb+sin(stag)*c/2)/Zd;solinv=s/c;
stagnew=57; %degrees from graph and from axyal direction
stagcor=90-stagnew; %put in my correct direction
correction=(stagcor-stagD); %degree i need to change in wicket gate
%% Rotor or Runner design
r_h=0;r_m=0.5;r_t=1;
r=r_h:0.01:r_t; 
Zr=5; %from graph
ah=0.1566;am=0.1049;at=0.1308; %rad, angle beetween mean line of profiles and chord
v2=v2m; %because it's in optimal condition design
v1m=v2m; %for conservation
R=r.*(Dt-Dh)/2+Dh/2;
U=omega*R;
W2=sqrt(U.^2+v2m^2);
gammar=atan(v2m./U);
v1t=leu./U;
w1t=v1t-U;
w1m=v1m;
W1=sqrt(w1t.^2+w1m^2);
for i=1:length(r)
    if w1t(i)>=0
        Winf=sqrt(((w1m+v2m)/2)^2+((abs(w1t)-U)/2).^2);
        A(i)=pi-atan((w1m+v2m)/(w1t(i)-U(i)));
    else
        Winf=sqrt(((w1m+v2m)/2)^2+((abs(w1t)+U)/2).^2);
        A(i)=atan((w1m+v2m)/(abs(w1t(i))+U(i)));
    end
end
ai=blade(ah,ah,am,am,at,at,r_h,r_m,r_t,r); %[rad] the angle beetween chord and midline of profiles along the blade
[CL,CD]=blade(hubcl,hubcd,midcl,midcd,tipcl,tipcd,r_h,r_m,r_t,r); %estruion of cl and cd along blade
stagbD=(gammar+ai)*180/pi;
incb=gammar+ai-A;incbD=incb*180/pi;
dr=R(2)-R(1);
for i=1:length(r)-1
    pwi(i)=(R(i)+dr/2)*2*pi*dr*v2m*leu;
    chord(i)=pwi(i)*1000/(omega*(R(i)+dr/2)*Zr*1000*dr*Winf(i)^2*(interp1(alfa,CL(:,i),incbD(i))*sin(A(i))-interp1(alfa,CD(:,i),incbD(i))*cos(A(i)))/2);
end
sb=2*pi*R(1:end-1)/Zr;
solinvb=sb./chord;%with weining correction, all the angles will start from 50 degree to 40 in the end, from axyal direction
solb=1./solinvb;
%% Losses evaluation 
v3=1; %m/s
Ydiff=0.8*(v2^2-v3^2)/(2*g); %[m]
Ydistr=0.07*vb^2/(2*g); %[m]
for i=1:length(r)-1
    y(i)=interp1(alfa,CD(:,i),incb(i))*cos((pi/2-gammar(i)))^2/((cos(incb(i))^3)*solinvb(i));
    Y(i)=y(i)*Winf(i)^2/(2*g);
    dp(i)=1000*leu+1000*Winf(i)^2*y(i)/2; %[Pa] of total pressure drop for each section
end
dP=sum(dp); %Total pressure drop in the rotor
Yrun=sum(Y);
eta=1-(Ydiff+Ydistr+Yrun)/Hm;
%% cavitation control
thoma=(ns^1.67)/exp(10.74);
NPSHr=thoma*Hm*eta_stat;
S=omega*sqrt(Q)/(g*NPSHr)^(3/4);
Patm=101325; %[Pa]
Pv_gas=10000; %[Pa] as a reference value of a vapour water pressure at 45°C
Dz=(Patm-Pv_gas)/(1000*g)+Ydiff+v3^2/(2*g)-NPSHr; % i need to stay under this value to avoid cavitation
%% conclusion





