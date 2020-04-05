%PROJECT 2 VAWT DESIGN
clc
clear all
close all

%% DATA
D=3+0.17*(19-15);
Np=2;
c=0.05*D/Np; % chord
H=2*D;
AR=H/c;
Cx0=2;
tol=1e-4;
e=0.7; % for numerical convergence
rho=1.225; %kg/mc
mu=1.81e-5; %Pa*s
at=1-sqrt(Cx0)/2;
Cxt=Cx0-4*(sqrt(Cx0)-1)*(1-at); % Glauert correction
load NACA1e4.txt;load NACA16e4.txt;load NACA36e4.txt;load NACA7e5.txt;load NACA2e6.txt;load NACA5e6.txt;

Nstream_vect = 10:10:40; % for discretization sensibility analysis

for qq=1:length(Nstream_vect)
Nstream=Nstream_vect(qq); %10; %20; %30; %40; % choice
R=D/2; 
dang=180/Nstream;

%% Variation of v0 and lambda
vmin=1;vmax=30; %m/s
lambdamin=1;lambdamax=10; %speed ratio
V0=vmin:1:vmax;
Lambda=lambdamin:1:lambdamax;
for k=1:length(Lambda) 
    lambda=Lambda(k);
    for j=1:length(V0)
        v0=V0(j);
        omega=lambda*v0/(D/2); %rad/s
        %% Upwind part %%
        angle=0:dang:180;
        teta1=angle(1:end-1)'*(pi/180);
        dteta=dang*pi/180;
        dA=H*(D/2)*dteta*sin(teta1); % forntal area of the stream tube
        teta=teta1(:)+dteta/2;
        a=zeros(length(teta),1); a(1:end)=0.25; % first guess
        aold=ones(length(teta),1);
        n=0;
        while max(abs(a-aold))>tol
            aold=a;
            n=n+1;
            vD=v0*(1-a);
            betainf=atan(vD.*sin(teta)./(vD.*cos(teta)+omega*R))*180/pi;
            winf=sqrt((vD.*sin(teta)).^2+(vD.*cos(teta)+omega*R).^2);
            Re=winf*c*rho/mu;
            naca1=real_aero(NACA1e4(:,1),NACA1e4(:,2),NACA1e4(:,3),betainf,AR);naca2=real_aero(NACA16e4(:,1),NACA16e4(:,2),NACA16e4(:,3),betainf,AR);
            naca3=real_aero(NACA36e4(:,1),NACA36e4(:,2),NACA36e4(:,3),betainf,AR);naca4=real_aero(NACA7e5(:,1),NACA7e5(:,2),NACA7e5(:,3),betainf,AR);
            naca5=real_aero(NACA2e6(:,1),NACA2e6(:,2),NACA2e6(:,3),betainf,AR);naca6=real_aero(NACA5e6(:,1),NACA5e6(:,2),NACA5e6(:,3),betainf,AR);
            for i=1:length(Re)
                if Re(i)>=1e4 && Re(i)<16e4
                    mat=Re_medio(naca1(i,2),naca2(i,2),naca1(i,3),naca2(i,3),Re(i),16e4,1e4);
                    betar(i,1)=(naca2(i,1)-naca1(i,1))/(16e4-1e4)*(Re(i)-1e4)+naca1(i,1);
                end
                if Re(i)>=16e4 && Re(i)<36e4
                    mat=Re_medio(naca2(i,2),naca3(i,2),naca2(i,3),naca3(i,3),Re(i),36e4,16e4);
                    betar(i,1)=(naca3(i,1)-naca2(i,1))/(36e4-16e4)*(Re(i)-16e4)+naca2(i,1);
                end 
                if Re(i)>=36e4 && Re(i)<7e5
                    mat=Re_medio(naca3(i,2),naca4(i,2),naca3(i,3),naca4(i,3),Re(i),7e5,36e4);
                    betar(i,1)=(naca4(i,1)-naca3(i,1))/(7e5-36e4)*(Re(i)-36e4)+naca3(i,1);
                end
                if Re(i)>=7e5 && Re(i)<2e6
                    mat=Re_medio(naca4(i,2),naca5(i,2),naca4(i,3),naca5(i,3),Re(i),2e6,7e5);
                    betar(i,1)=(naca5(i,1)-naca4(i,1))/(2e6-7e5)*(Re(i)-7e5)+naca4(i,1);
                end
                if Re(i)>=2e6 && Re(i)<=5e6
                    mat=Re_medio(naca5(i,2),naca6(i,2),naca5(i,3),naca6(i,3),Re(i),5e6,2e6);
                    betar(i,1)=(naca6(i,1)-naca5(i,1))/(5e6-2e6)*(Re(i)-2e6)+naca5(i,1);
                end
                if Re(i)<1e4
                    mat=Re_medio(naca1(i,2),naca2(i,2),naca1(i,3),naca2(i,3),1e4,16e4,1e4);
                    betar(i,1)=naca1(i,1);
                end
                if Re(i)>5e6
                    mat=Re_medio(naca5(i,2),naca6(i,2),naca5(i,3),naca6(i,3),5e6,5e6,2e6);
                    betar(i,1)=naca6(i,1);
                end
                cl(i,1)=mat(1);cd(i,1)=mat(2);
            end
            betar=betar*(pi/180);
            cn=cl.*cos(betar)+cd.*sin(betar);
            ct1=cl.*sin(betar)-cd.*cos(betar);
            Cxi=(((Np*c/(2*pi*R))*(winf.^2)/(v0^2))./abs(sin(teta))).*(cn.*sin(teta)-ct1.*cos(teta));
            delta=16.-16*Cxi;
            a=((4-sqrt(delta))/8)*e+(1-e)*aold;
            a(a>at)=(1-((Cx0-Cxi(a>at))./(4*(sqrt(Cx0)-1))))*e+(1-e)*aold(a>at);
        end
        %% Downwind part %%
        angle2=-180:dang:0;
        teta21=angle2(1:end-1)'*(pi/180);
        dA2=H*(D/2)*dteta*abs(sin(teta21));
        teta2=teta21+dteta/2;
        a2=zeros(length(teta),1);a2(1:end)=0.25;
        aold=ones(length(teta),1);
        n2=0;
        v02=(1-2*a)*v0;
        while max(abs(a2-aold))>tol
            aold=a2;
            n2=n2+1;
            vD=v02.*(1-a2);
            betainf=atan(vD.*sin(teta2)./(vD.*cos(teta2)+omega*R))*180/pi;
            winf2=sqrt((vD.*sin(teta2)).^2+(vD.*cos(teta2)+omega*R).^2);
            Re2=winf2*c*rho/mu;
            naca1=real_aero(NACA1e4(:,1),NACA1e4(:,2),NACA1e4(:,3),betainf,AR);naca2=real_aero(NACA16e4(:,1),NACA16e4(:,2),NACA16e4(:,3),betainf,AR);
            naca3=real_aero(NACA36e4(:,1),NACA36e4(:,2),NACA36e4(:,3),betainf,AR);naca4=real_aero(NACA7e5(:,1),NACA7e5(:,2),NACA7e5(:,3),betainf,AR);
            naca5=real_aero(NACA2e6(:,1),NACA2e6(:,2),NACA2e6(:,3),betainf,AR);naca6=real_aero(NACA5e6(:,1),NACA5e6(:,2),NACA5e6(:,3),betainf,AR);
            for i=1:length(Re2)
                if Re2(i)>=1e4 && Re2(i)<16e4
                    mat=Re_medio(naca1(i,2),naca2(i,2),naca1(i,3),naca2(i,3),Re2(i),16e4,1e4);
                    betar2(i,1)=(naca2(i,1)-naca1(i,1))/(16e4-1e4)*(Re2(i)-1e4)+naca1(i,1);
                end
                if Re2(i)>=16e4 && Re2(i)<36e4
                    mat=Re_medio(naca2(i,2),naca3(i,2),naca2(i,3),naca3(i,3),Re2(i),36e4,16e4);
                    betar2(i,1)=(naca3(i,1)-naca2(i,1))/(36e4-16e4)*(Re2(i)-16e4)+naca2(i,1);
                end 
                if Re2(i)>=36e4 && Re2(i)<7e5
                    mat=Re_medio(naca3(i,2),naca4(i,2),naca3(i,3),naca4(i,3),Re2(i),7e5,36e4);
                    betar2(i,1)=(naca4(i,1)-naca3(i,1))/(7e5-36e4)*(Re2(i)-36e4)+naca3(i,1);
                end
                if Re2(i)>=7e5 && Re2(i)<2e6
                    mat=Re_medio(naca4(i,2),naca5(i,2),naca4(i,3),naca5(i,3),Re2(i),2e6,7e5);
                    betar2(i,1)=(naca5(i,1)-naca4(i,1))/(2e6-7e5)*(Re2(i)-7e5)+naca4(i,1);
                end
                if Re2(i)>=2e6 && Re2(i)<=5e6
                    mat=Re_medio(naca5(i,2),naca6(i,2),naca5(i,3),naca6(i,3),Re2(i),5e6,2e6);
                    betar2(i,1)=(naca6(i,1)-naca5(i,1))/(5e6-2e6)*(Re2(i)-2e6)+naca5(i,1);
                end
                if Re2(i)<1e4
                    mat=Re_medio(naca1(i,2),naca2(i,2),naca1(i,3),naca2(i,3),1e4,16e4,1e4);
                    betar2(i,1)=naca1(i,1);
                end
                if Re2(i)>5e6
                    mat=Re_medio(naca5(i,2),naca6(i,2),naca5(i,3),naca6(i,3),5e6,5e6,2e6);
                    betar2(i,1)=naca6(i,1);
                end
                cl(i,1)=mat(1);cd(i,1)=mat(2);
            end
            betar2=betar2*(pi/180);
            cn=cl.*cos(betar2)+cd.*sin(betar2);
            ct2=cl.*sin(betar2)-cd.*cos(betar2);
            Cxi2=(((Np*c/(2*pi*R))*(winf2.^2)./(v02.^2))./abs(sin(teta2))).*(cn.*sin(teta2)-ct2.*cos(teta2));
            delta=(4*(1-2*a).^2).^2.-16*((1-2*a).^2).*Cxi2;
            a2=((4*(1-2*a).^2-sqrt(delta))./(8*(1-2*a).^2))*e+(1-e)*aold;
            a2(a2>at)=(1-((Cx0-Cxi2(a2>at))./(4*(sqrt(Cx0)-1))))*e+(1-e)*aold(a2>at);
        end
        %% Torque, Power and Cp evalation %%
        T1=((rho*c*H*(winf.^2)/2).*ct1)*R;
        T2=((rho*c*H*(winf2.^2)/2).*ct2)*R;
        fuffa=[T1;T2];
        power=((sum(omega*T1.*dA)+sum(omega*T2.*dA2))/(2*D*H))*Np; %W
        Cp(k,j)=power/(rho*D*H*v0^3/2);Power(k,j)=power/1000;%kW
        for b=1:length(fuffa)
            Torque(k,j,b)=fuffa(b);
        end
    end
end
mesh(V0,Lambda,Cp);hold on
title('Cp variation');xlabel('V_0 [m/s]');ylabel('\lambda - Speed Ratio');hold off;
 for i=1:length(Lambda)
     if Lambda(i)==round(Lambda(i))
         figure()
         hold on;
         for j=1:length(V0)
            ciccio(:) =  Torque(i,j,:);
            plot([teta*180/pi;teta2*180/pi+360],ciccio,'DisplayName',['v_0 =',num2str(V0(j)),'m/s']);
            hold on;    
         end
         legend('show');
         plot([teta*180/pi;teta2*180/pi+360],zeros(1,length(ciccio)),'--','DisplayName','Torque=0');
         title(['Torque @ \lambda =',num2str(Lambda(i))]);
         xlabel('angle [deg]');
         ylabel('Torque [Nm]');


         hold off;
     end
 end
maxT(qq) = max(max(max(Torque))); % taken as a reference for discretization sensibility analysis
    for kk =1:size(Torque(10,3,:),3)
        torque_ref(kk,qq)=Torque(10,3,kk);% taken as a reference for discretization sensibility analysis
    end
end

% figure()
% plot(Nstream_vect,maxT)
% title('Discretization sensibility analysis'); xlabel('Number of stream tubes'); ylabel('Maximum torque acieved [Nm]')

% figure()
% for ii=1:size(torque_ref,2)
%     z=sum(sum(torque_ref(:,ii)==0));
%     torque_ref_plot=torque_ref(1:(end-z),ii);
%     plot(linspace(0,360,length(torque_ref_plot)),torque_ref_plot); hold on
% end
% hold off;
% title('Torque dependancy from discretization - lambda=2, V0=3m/s'); xlabel('angle [deg]'); ylabel('Torque [Nm]');
% legend('10 str. tub.','20 str. tub.','30 str. tub.','40 str. tub.')

