%% PROJECT 1 - HAWT @ DESIGN CONDITIONS
clear all
close all
clc
%% DATA
D = 26.5 +(15-1)*3; R=D/2; %m
V0 = 9; % m/s Annual mean velocity at Pantelleria
Nb = 3; 
lambda_d = 7; % from statistical diagr 
Cp_stat = 0.48; % from statistical diagr 

omega_d = lambda_d*V0/R; % rad/s, first HP
rpm = omega_d*60/2/pi;
ro_aria=1.225; % kg/mc
visc_aria=1.81e-5; % Pa*s
ni=visc_aria/ro_aria;

r_H = 0.2; % *D/2
r_P = 0.7;
r_T = 0.96;
r = r_H:0.02:r_T;
% r_cut = (0.7:0.05:0.9); 
% settore di riferimento per linearizzazione intorno a 0.8

tol=1e-5;

% from wind.nrel.gov/airfoils/AirfoilFamilies for rotor R up to 50m HUB NREL S818, PRIM NREL S827, TIP NREL S828
load hub_1e6.txt;load pri_1e6.txt;load tip_1e6.txt;load hub_5e5.txt;load pri_5e5.txt;load tip_5e5.txt;

alfa=[-4:0.5:12.5,13.5:0.5:16];
hub1=[posvet(alfa,hub_5e5(:,1),hub_5e5(:,2));posvet(alfa,hub_5e5(:,1),hub_5e5(:,3))];
hub2=[posvet(alfa,hub_1e6(:,1),hub_1e6(:,2));posvet(alfa,hub_1e6(:,1),hub_1e6(:,3))];
pri1=[posvet(alfa,pri_5e5(:,1),pri_5e5(:,2));posvet(alfa,pri_5e5(:,1),pri_5e5(:,3))];
pri2=[posvet(alfa,pri_1e6(:,1),pri_1e6(:,2));posvet(alfa,pri_1e6(:,1),pri_1e6(:,3))];
tip1=[posvet(alfa,tip_5e5(:,1),tip_5e5(:,2));posvet(alfa,tip_5e5(:,1),tip_5e5(:,3))];
tip2=[posvet(alfa,tip_1e6(:,1),tip_1e6(:,2));posvet(alfa,tip_1e6(:,1),tip_1e6(:,3))];

%% OPTIMAL CHORD
x = (omega_d*(r*R)/V0)';
[a, a_p] = glauert(x);

%% Opt chord - H,P,T
%hub
Re=900000;
Re_H=800000;
nh=0;
x_H=r_H*R*omega_d/V0;
a_ph=a_p(x==x_H);
while abs(Re_H-Re)>tol
    nh=nh+1;
    Re=Re_H;
    a_pold=1;
    while abs(a_ph-a_pold)>tol
        a_pold=a_ph;
        phi=atan((1-a(x==x_H))/((1+a_ph)*x_H)); 
        %!!!non posso usare tan(phi)=a'/a perchè il triangolo delle velocità indotte 
        %non è parallelo alla direzione del lift essendoci drag, uso sempre il triangolo di Woo.
        W=(1-a(x==x_H))*V0/sin(phi);
        f=Nb*(R-r_H*R)/(2*r_H*R*sin(phi));
        F=2/pi*acos(exp(-f));
        if Re<=1e6 
                RE=Re;
        else
                RE=1e6;
        end
        medH=Re_medio(hub1(1,:),hub2(1,:),hub1(2,:),hub2(2,:),RE,1e6,5e5);
        eff=medH(1,:)./medH(2,:); % eff=CL/CD
        beta_inf=alfa(eff==max(eff));
        calettamento_H=phi*180/pi-beta_inf; %SU SLIDES BETA_C
        Cl=medH(1,eff==max(eff));Cd=medH(2,eff==max(eff));
        sigmaH=((4*F*a(x==x_H)*sin(phi)*tan(phi)/(1-a(x==x_H)))/(1+Cd*tan(phi)/Cl))/Cl;
        corda_H=sigmaH*2*pi*r_H*R/Nb;
        a_ph=sigmaH*W*(Cl-Cd/tan(phi))/(4*omega_d*r_H*R*F);
    end
    Re_H=W*corda_H/ni;
    
end 

%PRIMARY

Re=900000;
Re_P=800000;
np=0;
x_P=r_P*R*omega_d/V0;
a_pp=a_p(x==x_P);
while abs(Re_P-Re)>tol
    np=np+1;
    Re=Re_P;
    a_pold=1;
    while abs(a_pp-a_pold)>tol
        a_pold=a_pp;
        phi=atan((1-a(x==x_P))/((1+a_pp)*x_P)); %!!!non posso usare il semplice rapporto a'/a perchè il triangolo delle velocità indotte non è parallelo alla direzione del lift con il drag, quindi uso sempre il triangolo effettivo.
        W=(1-a(x==x_P))*V0/sin(phi);
        f=Nb*(R-r_P*R)/(2*r_P*R*sin(phi));
        F=2/pi*acos(exp(-f));
        if Re<=1e6 
                RE=Re;
        else
                RE=1e6;
        end
        medP=Re_medio(pri1(1,:),pri2(1,:),pri1(2,:),pri2(2,:),RE,1e6,5e5);
        eff=medP(1,:)./medP(2,:);
        beta_inf=alfa(eff==max(eff));
        calettamento_P=phi*180/pi-beta_inf;
        Cl=medP(1,eff==max(eff));Cd=medP(2,eff==max(eff));
        sigmaP=((4*F*a(x==x_P)*sin(phi)*tan(phi)/(1-a(x==x_P)))/(1+Cd*tan(phi)/Cl))/Cl;
        corda_P=sigmaP*2*pi*r_P*R/Nb;
        a_pp=sigmaP*W*(Cl-Cd/tan(phi))/(4*omega_d*r_P*R*F);
    end
    Re_P=W*corda_P/ni;
    
end

% TIP

Re=900000;
Re_T=800000;
nt=0;
x_T=r_T*R*omega_d/V0;
a_pt=a_p(x==x_T);
while abs(Re_T-Re)>tol
    nt=nt+1;
    Re=Re_T;
    a_pold=1;
    while abs(a_pt-a_pold)>tol
        a_pold=a_pt;
        phi=atan((1-a(x==x_T))/((1+a_pt)*x_T)); %!!!non posso usare il semplice rapporto a'/a perchè il triangolo delle velocità indotte non è parallelo alla direzione del lift con il drag, quindi uso sempre il triangolo effettivo.
        W=(1-a(x==x_T))*V0/sin(phi);
        f=Nb*(R-r_T*R)/(2*r_T*R*sin(phi));
        F=2/pi*acos(exp(-f));
        if Re<=1e6 
                RE=Re;
        else
                RE=1e6;
        end
        medT=Re_medio(tip1(1,:),tip2(1,:),tip1(2,:),tip2(2,:),RE,1e6,5e5);
        eff=medT(1,:)./medT(2,:);
        beta_inf=alfa(eff==max(eff));
        calettamento_T=phi*180/pi-beta_inf;
        Cl=medT(1,eff==max(eff));Cd=medT(2,eff==max(eff));
        sigmaT=((4*F*a(x==x_T)*sin(phi)*tan(phi)/(1-a(x==x_T)))/(1+Cd*tan(phi)/Cl))/Cl;
        corda_T=sigmaT*2*pi*r_T*R/Nb;
        a_pt=sigmaT*W*(Cl-Cd/tan(phi))/(4*omega_d*r_T*R*F);
    end
    Re_T=W*corda_T/ni;
    
end 
%% Opt chord - BLADESPAN
[CL, CD]=blade(medH(1,:),medH(2,:),medP(1,:),medP(2,:),medT(1,:),medT(2,:),r_H,r_P,r_T,r);
EFF=CL./CD;
BETAinf=ones(1,length(r));
sigma=BETAinf;
corda=BETAinf;
for i=1:length(r)
    BETAinf(i)=alfa(EFF(:,i)==max(EFF(:,i)));
end
a_pi=a_p;
for i=1:length(r)
    a_pold=1;
    while abs(a_pi(i)-a_pold)>tol
        a_pold=a_pi(i);
        phi(i)=atan((1-a(i))/((1+a_pi(i))*x(i)));
        W(i)=(1-a(i))*V0/sin(phi(i));
        f(i)=Nb*(R-r(i)*R)/(2*r(i)*R*sin(phi(i)));
        F(i)=2/pi*acos(exp(-f(i)));
        Cl(i)=CL(EFF(:,i)==max(EFF(:,i)),i);Cd(i)=CD(EFF(:,i)==max(EFF(:,i)),i);
        sigma(i)=((4*F(i)*a(i)*sin(phi(i))*tan(phi(i))/(1-a(i)))/(1+Cd(i)*tan(phi(i))/Cl(i)))/Cl(i);
        corda(i)=sigma(i)*2*pi*r(i)*R/Nb;
        a_pi(i)=sigma(i)*W(i)*(Cl(i)-Cd(i)/tan(phi(i)))/(4*omega_d*r(i)*R*F(i));
    end
end
calettamento=phi*180/pi-BETAinf;

%% Opt chord - PLOTS
% figure(1)
%     plot(r,corda,'LineWidth',4);hold on;
%     p1 = plot(r_H,corda_H,'o','MarkerSize',10,'MarkerFaceColor','r'); 
%     p2 = plot(r_P,corda_P,'o','MarkerSize',10,'MarkerFaceColor','y'); 
%     p3 = plot(r_T,corda_T,'o','MarkerSize',10,'MarkerFaceColor','g');
%     hold off;
%     legend([p1 p2 p3],'hub','primary','tip');
%     xlabel('r/R'); ylabel('chord [m]');
%     title('Optimized chord'); saveas(gcf,'OptChordPlot','jpg')
% figure(2)
%     plot(r,calettamento,'LineWidth',4);hold on;
%     p1 = plot(r_H,calettamento_H,'o','MarkerSize',10,'MarkerFaceColor','r');
%     p2 = plot(r_P,calettamento_P,'o','MarkerSize',10,'MarkerFaceColor','y');
%     p3 = plot(r_T,calettamento_T,'o','MarkerSize',10,'MarkerFaceColor','g');
%     hold off;
%     legend([p1 p2 p3],'hub','primary','tip');
%     xlabel('r/R'); title('Stagger angle [deg]'); saveas(gcf,'OptPitchPlot','jpg')

%% LINEARIZED CHORD

% TBN: "it's kitchen!" 
% Massima corda ha c/R<=0.12
c1=corda(1);
if c1/R<=0.12
    c1=corda(1);
else
    c1=0.12*R;
end
c2=corda(x==x_P);
mc=(c1-c2)/(r_H-r_P);
cordalin=(r-r_H)*mc+c1;


sol=(cordalin*Nb)./(2*pi*R*r); % sol stays for solidity
e=0.5; % variable for convergence of cycle
%% lin chord - H,P,T

%hub

aH=a(1);
a_pH=a_p(1);
aold=1;
a_pold=1;
nh=0;
while abs(a_pH-a_pold)>tol 
    nh=nh+1;
    a_pold=a_pH;
    while abs(aH-aold)>tol
        aold=aH;
        PHI=atan((1-aH)/((1+a_pH)*x_H));
        Winf=(1-aH)*V0/sin(PHI);
        f=Nb*(R-r_H*R)/(2*r_H*R*sin(PHI));
        F=2/pi*acos(exp(-f));
        ReyH=Winf*cordalin(x==x_H)/ni;
            if ReyH<=1e6
                RE=ReyH;
            else
                RE=1e6;
            end
        medH=Re_medio(hub1(1,:),hub2(1,:),hub1(2,:),hub2(2,:),RE,1e6,5e5);
        eff=medH(1,:)./medH(2,:);
        beta_inf=alfa(eff==max(eff));
        beta_cH=PHI*180/pi-beta_inf;
        Cl=medH(1,eff==max(eff));Cd=medH(2,eff==max(eff));
        aH=sol(x==x_H)*Winf*(Cl/tan(PHI)+Cd)/(4*V0*F);
        aH=e*aH+(1-e)*aold;
    end
    a_pH=sol(x==x_H)*Winf*(Cl-Cd/tan(PHI))/(4*omega_d*r_H*R*F);
    a_pH=e*a_pH+(1-e)*a_pold;
end

%primary

aP=a(x==x_P);
a_pP=a_p(x==x_P);
aold=1;
a_pold=1;
np=0;
while abs(a_pP-a_pold)>tol 
    np=np+1;
    a_pold=a_pP;
    while abs(aP-aold)>tol
        aold=aP;
        PHI=atan((1-aP)/((1+a_pP)*x_P));
        Winf=(1-aP)*V0/sin(PHI);
        f=Nb*(R-r_P*R)/(2*r_P*R*sin(PHI));
        F=2/pi*acos(exp(-f));
        ReyP=Winf*cordalin(x==x_P)/ni;
            if ReyP<=1e6
                RE=ReyP;
            else
                RE=1e6;
            end
        medP=Re_medio(pri1(1,:),pri2(1,:),pri1(2,:),pri2(2,:),RE,1e6,5e5);
        eff=medP(1,:)./medP(2,:);
        beta_inf=alfa(eff==max(eff));
        beta_cP=PHI*180/pi-beta_inf;
        Cl=medP(1,eff==max(eff));Cd=medP(2,eff==max(eff));
        aP=sol(x==x_P)*Winf*(Cl/tan(PHI)+Cd)/(4*V0*F);
        aP=e*aP+(1-e)*aold;
    end
    a_pP=sol(x==x_P)*Winf*(Cl-Cd/tan(PHI))/(4*omega_d*r_P*R*F);
    a_pP=e*a_pP+(1-e)*a_pold;
end

%tip

aT=a(x==x_T);
a_pT=a_p(x==x_T);
aold=1;
a_pold=1;
nt=0;
while abs(a_pT-a_pold)>tol 
    nt=nt+1;
    a_pold=a_pT;
    while abs(aT-aold)>tol
        aold=aT;
        PHI=atan((1-aT)/((1+a_pT)*x_T));
        Winf=(1-aT)*V0/sin(PHI);
        f=Nb*(R-r_T*R)/(2*r_T*R*sin(PHI));
        F=2/pi*acos(exp(-f));
        ReyT=Winf*cordalin(x==x_T)/ni;
            if ReyT<=1e6
                RE=ReyT;
            else
                RE=1e6;
            end
        medT=Re_medio(tip1(1,:),tip2(1,:),tip1(2,:),tip2(2,:),RE,1e6,5e5);
        eff=medT(1,:)./medT(2,:);
        beta_inf=alfa(eff==max(eff));
        beta_cT=PHI*180/pi-beta_inf;
        Cl=medT(1,eff==max(eff));Cd=medT(2,eff==max(eff));
        aT=sol(x==x_T)*Winf*(Cl/tan(PHI)+Cd)/(4*V0*F);
        aT=e*aT+(1-e)*aold;
    end
    a_pT=sol(x==x_T)*Winf*(Cl-Cd/tan(PHI))/(4*omega_d*r_T*R*F);
    a_pT=e*a_pT+(1-e)*a_pold;
end

%% lin chord - BLADESPAN
[CL, CD]=blade(medH(1,:),medH(2,:),medP(1,:),medP(2,:),medT(1,:),medT(2,:),r_H,r_P,r_T,r);
EFF=CL./CD;
BETAinf=ones(1,length(r));
for i=1:length(r)
    BETAinf(i)=alfa(EFF(:,i)==max(EFF(:,i)));
    Cl(i)=CL(EFF(:,i)==max(EFF(:,i)),i);Cd(i)=CD(EFF(:,i)==max(EFF(:,i)),i);
end
ai=ones(1,length(r));
n=ai;a_pi=ai;n2=ai;
for i=1:length(r)
    ai(i)=a(i);
    a_pi(i)=a_p(i);
    aold=1;
    a_pold=1;
    n(i)=0;
    n2(i)=0;
    while abs(a_pi(i)-a_pold)>tol 
        n(i)=n(i)+1;
        a_pold=a_pi(i);
        while abs(ai(i)-aold)>tol
        aold=ai(i);
        n2(i)=n2(i)+1;
        PHI(i)=atan((1-ai(i))/((1+a_pi(i))*x(i)));
        Winf(i)=(1-ai(i))*V0/sin(PHI(i));
        f=Nb*(R-r(i)*R)/(2*r(i)*R*sin(PHI(i)));
        F(i)=2/pi*acos(exp(-f));
        beta_ci(i)=PHI(i)*180/pi-BETAinf(i);
        ai(i)=sol(i)*Winf(i)*(Cl(i)/tan(PHI(i))+Cd(i))/(4*V0*F(i));
        ai(i)=e*ai(i)+(1-e)*aold;
        end
        a_pi(i)=sol(i)*Winf(i)*(Cl(i)-Cd(i)/tan(PHI(i)))/(4*omega_d*r(i)*R*F(i));
        a_pi(i)=e*a_pi(i)+(1-e)*a_pold;
    end
end

%% Lin chord - PLOTS
% figure (3)
%     plot(r,cordalin,'LineWidth',4);hold on;
%     plot(r,corda);
%     hold off;
%     xlabel('r/R'); ylabel('corda [m]');
%     title('Linearization of blade cord'); saveas(gcf,'LinChordPlot','jpg')
% 
% figure(4)
%     plot(r,calettamento);hold on;
%     plot(r,beta_ci);hold on; 
%     plot(r_H,beta_cH,'o');
%     plot(r_P,beta_cP,'o');
%     plot(r_T,beta_cT,'o');
%     hold off;
%     legend('opt Cp','linearized','hub','primary','tip');
%     xlabel('r/R'); title('stagger angle [deg]');
%% POWER & Cp
Cx=ones(length(r),1);
% dr=r(i+1)-r(i); for us, dr is a constant difference for construction
% because r=r_H:0.05:r_T and so dr=0.05
dr=r(2)-r(1);
for i=1:length(r)
    Cx(i)=cordalin(i)*(Winf(i)^2)*(Cl(i)*sin(PHI(i))-Cd(i)*cos(PHI(i)))*r(i)*dr*(R^2);
end
Cp=Nb*omega_d*sum(Cx)/(pi*R^2*V0^3);
TotPower=ro_aria*Nb*omega_d*sum(Cx)/2;


% Calculation of maximun power and Cp using angular momentum balance with
% F correction
Cx2=ones(length(r),1);
% Taking into account also the dx effect, also x depends on r.
% It is possibile to define dx=dr*omega*R/V0, using the previous value for dr
dx=dr*omega_d*R/V0;
for i=1:length(r)
    Cx2(i)=F(i)*a_pi(i)*(1-ai(i))*(x(i)^3)*dx;
end
Cp2=8*sum(Cx2)/(lambda_d^2);
TotPower2=4*ro_aria*pi*(R^2)*(V0^3)*sum(Cx2)/(lambda_d^2); % W


% The two calculation must give the same results
Cpdev=abs(Cp-Cp2)*100/Cp; % Differenza percentuale tra i Cp con i due metodi
Powdev=abs(TotPower-TotPower2)*100/TotPower; % Differenza percentuale tra la potenze con i due metodi
Cp = (Cp+Cp2)/2; % Adimensional vs Cp=0.48 as maximum value for statistics 
                 % (evaluated as the mean between the two Cp fuond with
                 % power and momentum balance.
TotPower=(TotPower+TotPower2)/2/1000; % kW

% SHOW RESULTS
fprintf('\nFinal results: \nDiameter = %f m \nCoefficient of performance = %f \nPower = %f kW\n',D,Cp,TotPower);
%% evaluation of axial thrust
for i=1:length(r)
    fx(i)=(Nb*ro_aria*(Winf(i)^2)*cordalin(i)*(Cl(i)*cos(PHI(i))+Cd(i)*sin(PHI(i)))*dr*R)/2;
    fx2(i)=4*pi*r(i)*R*ro_aria*V0^2*ai(i)*(1-ai(i))*F(i)*dr*R;
end
Fa=sum(fx)/1000;
Fa2=sum(fx2)/1000;
Ct=Fa*1000/(ro_aria*pi*(R^2)*(V0^2)/2);
fprintf('\nAxial force = %f kN \nThrust Coefficient = %f \n',Fa,Ct);



