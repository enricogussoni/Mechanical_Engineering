%% PROJECT 1.2 - HAWT @ OFF DESIGN CONDITIONS
%i dati importati da porject 1_1 in clearvars -EXCEPT
%clearvars -EXCEPT V0 lambda_d omega_d ai a_pi R r r_H r_P r_T ...
     %tol cordalin  Nb ni...
     %beta_cH beta_cP beta_cT e sol
close all
clc
%% DATI
load hub_1e6.txt;load pri_1e6.txt;load tip_1e6.txt;load hub_5e5.txt;load pri_5e5.txt;load tip_5e5.txt;

angle=-9:0.1:12;ALFA=angle(end):0.1:90; %we are composing th vectors of cl and cd using also Viterna-Corrigan relation to extend the range over the stall angle
CLhub1 = interp1(hub_5e5(:,1),hub_5e5(:,2),angle);% we have chosen
CLoldold=CLhub1; angleoldold=angle;
CDhub1 = interp1(hub_5e5(:,1),hub_5e5(:,3),angle);fuffa=vit_corr(angle(end),CLhub1(angle==angle(end)),CDhub1(angle==angle(end)),1000,ALFA);
CLhub1=[CLhub1(1:end-1),fuffa(:,1)'];CDhub1=[CDhub1(1:end-1),fuffa(:,2)'];
CLpri1 = interp1(pri_5e5(:,1),pri_5e5(:,2),angle);
CDpri1 = interp1(pri_5e5(:,1),pri_5e5(:,3),angle);fuffa=vit_corr(angle(end),CLpri1(angle==angle(end)),CDpri1(angle==angle(end)),1000,ALFA);
CLpri1=[CLpri1(1:end-1),fuffa(:,1)'];CDpri1=[CDpri1(1:end-1),fuffa(:,2)'];
CLtip1 = interp1(tip_5e5(:,1),tip_5e5(:,2),angle);
CDtip1 = interp1(tip_5e5(:,1),tip_5e5(:,3),angle);fuffa=vit_corr(angle(end),CLtip1(angle==angle(end)),CDtip1(angle==angle(end)),1000,ALFA);
CLtip1=[CLtip1(1:end-1),fuffa(:,1)'];CDtip1=[CDtip1(1:end-1),fuffa(:,2)'];
CLhub2 = interp1(hub_1e6(:,1),hub_1e6(:,2),angle);
CDhub2 = interp1(hub_1e6(:,1),hub_1e6(:,3),angle);fuffa=vit_corr(angle(end),CLhub2(angle==angle(end)),CDhub2(angle==angle(end)),1000,ALFA);
CLhub2=[CLhub2(1:end-1),fuffa(:,1)'];CDhub2=[CDhub2(1:end-1),fuffa(:,2)'];
CLpri2 = interp1(pri_1e6(:,1),pri_1e6(:,2),angle);
CDpri2 = interp1(pri_1e6(:,1),pri_1e6(:,3),angle);fuffa=vit_corr(angle(end),CLpri2(angle==angle(end)),CDpri2(angle==angle(end)),1000,ALFA);
CLpri2=[CLpri2(1:end-1),fuffa(:,1)'];CDpri2=[CDpri2(1:end-1),fuffa(:,2)'];
CLtip2 = interp1(tip_1e6(:,1),tip_1e6(:,2),angle);
CDtip2 = interp1(tip_1e6(:,1),tip_1e6(:,3),angle);fuffa=vit_corr(angle(end),CLtip2(angle==angle(end)),CDtip2(angle==angle(end)),1000,ALFA);
CLtip2=[CLtip2(1:end-1),fuffa(:,1)'];CDtip2=[CDtip2(1:end-1),fuffa(:,2)'];
angle=[angle,ALFA(2:end)];
CLold=CLhub1; angleold=angle;

    

% Extension of data (for Pitch-to-Stall)
[CLhub1, CDhub1, angle1] = allungangle(CLhub1, CDhub1, angle);
[CLhub2, CDhub2, angle2] = allungangle(CLhub2, CDhub2, angle);
[CLtip1, CDtip1, angle3] = allungangle(CLtip1, CDtip1, angle);
[CLtip2, CDtip2, angle4] = allungangle(CLtip2, CDtip2, angle);
[CLpri1, CDpri1, angle5] = allungangle(CLpri1, CDpri1, angle);
[CLpri2, CDpri2, angle6] = allungangle(CLpri2, CDpri2, angle);

amin = max([angle1(1),angle2(1),angle3(1),angle4(1),angle5(1),angle6(1)]);

angle = amin:0.1:angle(end);
CLhub1= CLhub1(find(abs(angle1-amin)<1e-5):end); CLhub2= CLhub2(find(abs(angle2-amin)<1e-5):end);
CDhub1= CDhub1(find(abs(angle1-amin)<1e-5):end); CDhub2= CDhub2(find(abs(angle2-amin)<1e-5):end);
CLpri1= CLpri1(find(abs(angle5-amin)<1e-5):end); CLpri2= CLpri2(find(abs(angle6-amin)<1e-5):end);
CDpri1= CDpri1(find(abs(angle5-amin)<1e-5):end); CDpri2= CDpri2(find(abs(angle6-amin)<1e-5):end);
CLtip1= CLtip1(find(abs(angle3-amin)<1e-5):end); CLtip2= CLtip2(find(abs(angle4-amin)<1e-5):end);
CDtip1= CDtip1(find(abs(angle3-amin)<1e-5):end); CDtip2= CDtip2(find(abs(angle4-amin)<1e-5):end);

% figure()
%     plot(angle,CLhub1,'LineWidth',2); hold on;
%     plot(angleold,CLold,'LineWidth',2); hold on;
%     plot(angleoldold,CLoldold,'LineWidth',2); hold off;
%     title('C_L hub @ Re=0.5e6');
%     legend('mirroring','V-C','Original Data');
%     ylabel('C_L');
%     xlabel('attack angle [deg]');

Vcutin=3.5; %m/s %alta perchè macchina grossa
Vcutout = 28;%m/s, alta perchè su isola
Vc = V0;%.txt
Vr = 15;% from wind energy handbook jw and sons, 2001
step_v = 1;
vAC = (Vcutin+0.5):step_v:Vc; % partiamo 0.5 dopo cutin per questione vista in classe a slide 4
vCD = (Vc+step_v):step_v:Vr;
vDE = (Vr+step_v):step_v:Vcutout;

ad = ai;
apd = a_pi;
betacd=beta_ci;

e=0.5;
%% RANGE A-C
omega=ones(1,length(vAC));
CpAC=omega;PowerAC=omega;CpAC2=CpAC;PowerAC2=PowerAC;
for j=1:length(vAC)
    omega(j) = lambda_d*vAC(j)/R;
    xACj = (omega(j)*(r*R)/vAC(j));% non faccio xAC(j) perchè non serve salvarsi matrice, sovrascrive j-esimo vettore
    x_h = xACj(r==r_H);
    x_p = xACj(r==r_P);
    x_t = xACj(r==r_T);
    %% a-c: H,P,T
    %hub
    ah=ad(r==r_H);
    aph=apd(r==r_H);
    aold=1;
    apold=1;
    nh=0;
    while abs(aph-apold)>tol 
        nh=nh+1
        apold=aph;
        while abs(ah-aold)>tol
            aold=ah;
            phi=atan((1-ah)/((1+aph)*x_h));
            w=(1-ah)*vAC(j)/sin(phi);
            f=Nb*(R-r_H*R)/(2*r_H*R*sin(phi));
            F=2/pi*acos(exp(-f));
            Reh=w*cordalin(xACj==x_h)/ni;
                if Reh<=1e6
                    RE=Reh;
                else
                    RE=1e6;
                end
                if RE>=5e5
                    RE=RE;
                else
                    RE=5e5;
                end
            medh=Re_medio(hub1(1,:),hub2(1,:),hub1(2,:),hub2(2,:),RE,1e6,5e5);
            ah=sol(xACj==x_h)*w*(medh(1,alfa==BETAinf(r==r_H))/tan(phi)+medh(2,alfa==BETAinf(r==r_H)))/(4*vAC(j)*F);
            ah=e*ah+(1-e)*aold;
        end
        aph=sol(xACj==x_h)*w*(medh(1,alfa==BETAinf(r==r_H))-medh(2,alfa==BETAinf(r==r_H))/tan(phi))/(4*omega(j)*r_H*R*F);
        aph=e*aph+(1-e)*apold;
    end
    %primary
    ap=ad(r==r_P);
    app=apd(r==r_P);
    aold=1;
    apold=1;
    np=0;
    while abs(app-apold)>tol 
        np=np+1
        apold=app;
        while abs(ap-aold)>tol
            aold=ap;
            phi=atan((1-ap)/((1+app)*x_p));
            w=(1-ap)*vAC(j)/sin(phi);
            f=Nb*(R-r_P*R)/(2*r_P*R*sin(phi));
            F=2/pi*acos(exp(-f));
            Rep=w*cordalin(xACj==x_p)/ni;
                if Rep<=1e6
                    RE=Rep;
                else
                    RE=1e6;
                end
                if RE>=5e5
                    RE=RE;
                else
                    RE=5e5;
                end
            medp=Re_medio(pri1(1,:),pri2(1,:),pri1(2,:),pri2(2,:),RE,1e6,5e5);
            ap=sol(xACj==x_p)*w*(medp(1,alfa==BETAinf(r==r_P))/tan(phi)+medp(2,alfa==BETAinf(r==r_P)))/(4*vAC(j)*F);
            ap=e*ap+(1-e)*aold;
        end
        app=sol(xACj==x_p)*w*(medp(1,alfa==BETAinf(r==r_P))-medp(2,alfa==BETAinf(r==r_P))/tan(phi))/(4*omega(j)*r_P*R*F);
        app=e*app+(1-e)*apold;
    end
    %tip
    at=ad(r==r_T);
    apt=apd(r==r_T);
    aold=1;
    apold=1;
    nt=0;
    while abs(apt-apold)>tol 
        nt=nt+1
        apold=apt;
        while abs(at-aold)>tol
            aold=at;
            phi=atan((1-at)/((1+apt)*x_t));
            w=(1-at)*vAC(j)/sin(phi);
            f=Nb*(R-r_T*R)/(2*r_T*R*sin(phi));
            F=2/pi*acos(exp(-f));
            Ret=w*cordalin(xACj==x_t)/ni;
                if Ret<=1e6
                    RE=Ret;
                else
                    RE=1e6;
                end
                if RE>=5e5
                    RE=RE;
                else
                    RE=5e5;
                end
            medt=Re_medio(tip1(1,:),tip2(1,:),tip1(2,:),tip2(2,:),RE,1e6,5e5);
            at=sol(xACj==x_t)*w*(medt(1,alfa==BETAinf(r==r_T))/tan(phi)+medt(2,alfa==BETAinf(r==r_T)))/(4*vAC(j)*F);
            at=e*at+(1-e)*aold;
        end
        apt=sol(xACj==x_t)*w*(medt(1,alfa==BETAinf(r==r_T))-medt(2,alfa==BETAinf(r==r_T))/tan(phi))/(4*omega(j)*r_T*R*F);
        apt=e*apt+(1-e)*apold;
    end

    
    %% a-c: BLADESPAN
    % blade distribution of cl and cd
    [cl, cd]=blade(medh(1,:),medh(2,:),medp(1,:),medp(2,:),medt(1,:),medt(2,:),r_H,r_P,r_T,r);
    % calculation of a and a' for all blade's parts
    aAC=ones(1,length(r));
    n=aAC;apAC=aAC;n2=aAC;
    for i=1:length(r)
        aAC(i)=ad(i);
        apAC(i)=apd(i);
        aold=1;
        apold=1;
        n(i)=0;
        n2(i)=0;
        while abs(apAC(i)-apold)>tol 
            n(i)=n(i)+1;
            apold=apAC(i);
            while abs(aAC(i)-aold)>tol
                aold=aAC(i);
                n2(i)=n2(i)+1;
                phi(i)=atan((1-aAC(i))/((1+apAC(i))*xACj(i)));
                w(i)=(1-aAC(i))*vAC(j)/sin(phi(i));
                f=Nb*(R-r(i)*R)/(2*r(i)*R*sin(phi(i)));
                F_AC(i)=2/pi*acos(exp(-f)); 
                aAC(i)=sol(i)*w(i)*(cl(alfa==BETAinf(r==r(i)),i)/tan(phi(i))+cd(alfa==BETAinf(r==r(i)),i))/(4*vAC(j)*F_AC(i));
                aAC(i)=e*aAC(i)+(1-e)*aold;
            end
            apAC(i)=sol(i)*w(i)*(cl(alfa==BETAinf(r==r(i)),i)-cd(alfa==BETAinf(r==r(i)),i)/tan(phi(i)))/(4*omega(j)*r(i)*R*F_AC(i));
            apAC(i)=e*apAC(i)+(1-e)*apold;
        end
    end
    % power and Cp calculation 
    cc=ones(length(r),1);
    for i=1:length(r)
        cc(i)=cordalin(i)*(w(i)^2)*(cl(alfa==BETAinf(r==r(i)),i)*sin(phi(i))-cd(alfa==BETAinf(r==r(i)),i)*cos(phi(i)))*r(i)*dr*(R^2);
    end
    CpAC(j)=Nb*omega(j)*sum(cc)/(pi*R^2*vAC(j)^3);
    PowerAC(j)=ro_aria*Nb*omega(j)*sum(cc)/2;
    cc2=ones(length(r),1);
    for i=1:length(r)
        cc2(i)=F_AC(i)*apAC(i)*(1-aAC(i))*(xACj(i)^3)*dx;
    end
    CpAC2(j)=8*sum(cc2)/(lambda_d^2);
    PowerAC2(j)=4*ro_aria*pi*(R^2)*(vAC(j)^3)*sum(cc2)/(lambda_d^2);
end
%% Range C-D
betac=betacd;
daB=0.1;
for j=1:length(vCD) 
    bool=ones(1,length(r));
    xCDj = (omega_d*(r*R)/vCD(j));% non faccio xAC(j) perchè non serve salvarsi matrice, sovrascrive j-esimo vettore
    x_h = xCDj(r==r_H);
    x_p = xCDj(r==r_P);
    x_t = xCDj(r==r_T);
    while sum(bool)>=2/10*length(r)
        %% c-d: H,P,T
        %hub
        ah=ad(r==r_H);
        aph=apd(r==r_H);
        aold=1;
        apold=1;
        nh=0;
        while abs(aph-apold)>tol 
            nh=nh+1;
            apold=aph;
            while abs(ah-aold)>tol
                aold=ah;
                phi=atan((1-ah)/((1+aph)*x_h));
                w=(1-ah)*vCD(j)/sin(phi);
                f=Nb*(R-r_H*R)/(2*r_H*R*sin(phi));
                F=2/pi*acos(exp(-f));
                Reh=w*cordalin(xCDj==x_h)/ni;
                    if Reh<=1e6
                        RE=Reh;
                    else
                        RE=1e6;
                    end
                    if RE>=5e5
                        RE=RE;
                    else
                        RE=5e5;
                    end
                medh=Re_medio(CLhub1,CLhub2,CDhub1,CDhub2,RE,1e6,5e5);
                betainf=phi*180/pi-betac(r==r_H);
                ah=sol(xCDj==x_h)*w*(interp1(angle,medh(1,:),betainf)/tan(phi)+interp1(angle,medh(2,:),betainf))/(4*vCD(j)*F);
                ah=e*ah+(1-e)*aold;
            end
            aph=sol(xCDj==x_h)*w*(interp1(angle,medh(1,:),betainf)-interp1(angle,medh(2,:),betainf)/tan(phi))/(4*omega_d*r_H*R*F);
            aph=e*aph+(1-e)*apold;
        end 
      
        %primary
        ap=ad(r==r_P);
        app=apd(r==r_P);
        aold=1;
        apold=1;
        np=0;
        while abs(app-apold)>tol 
            np=np+1;
            apold=app;
            while abs(ap-aold)>tol
                aold=ap;
                phi=atan((1-ap)/((1+app)*x_p));
                w=(1-ap)*vCD(j)/sin(phi);
                f=Nb*(R-r_P*R)/(2*r_P*R*sin(phi));
                F=2/pi*acos(exp(-f));
                Rep=w*cordalin(xCDj==x_p)/ni;
                    if Rep<=1e6
                        RE=Rep;
                    else
                        RE=1e6;
                    end
                    if RE>=5e5
                        RE=RE;
                    else
                        RE=5e5;
                    end
                medp=Re_medio(CLpri1,CLpri2,CDpri1,CDpri2,RE,1e6,5e5);
                betainf=phi*180/pi-betac(r==r_P);
                ap=sol(xCDj==x_p)*w*(interp1(angle,medp(1,:),betainf)/tan(phi)+interp1(angle,medp(2,:),betainf))/(4*vCD(j)*F);
                ap=e*ap+(1-e)*aold;
            end
            app=sol(xCDj==x_p)*w*(interp1(angle,medp(1,:),betainf)-interp1(angle,medp(2,:),betainf)/tan(phi))/(4*omega_d*r_P*R*F);
            app=e*app+(1-e)*apold;
        end
        %tip
        at=ad(r==r_T);
        apt=apd(r==r_T);
        aold=1;
        apold=1;
        nt=0;
        while abs(apt-apold)>tol 
            nt=nt+1;
            apold=apt;
            while abs(at-aold)>tol
                aold=at;
                phi=atan((1-at)/((1+apt)*x_t));
                w=(1-at)*vCD(j)/sin(phi);
                f=Nb*(R-r_T*R)/(2*r_T*R*sin(phi));
                F=2/pi*acos(exp(-f));
                Ret=w*cordalin(xCDj==x_t)/ni;
                    if Ret<=1e6
                        RE=Ret;
                    else
                        RE=1e6;
                    end
                    if RE>=5e5
                        RE=RE;
                    else
                        RE=5e5;
                    end
                medt=Re_medio(CLtip1,CLtip2,CDtip1,CDtip2,RE,1e6,5e5);
                betainf=phi*180/pi-betac(r==r_T);
                at=sol(xCDj==x_t)*w*(interp1(angle,medt(1,:),betainf)/tan(phi)+interp1(angle,medt(2,:),betainf))/(4*vCD(j)*F);
                at=e*at+(1-e)*aold;
            end
            apt=sol(xCDj==x_t)*w*(interp1(angle,medt(1,:),betainf)-interp1(angle,medt(2,:),betainf)/tan(phi))/(4*omega_d*r_T*R*F);
            apt=e*apt+(1-e)*apold;
        end
        
        %% c-d: BLADESPAN
        % blade distribution of cl and cd
        [cl, cd]=blade(medh(1,:),medh(2,:),medp(1,:),medp(2,:),medt(1,:),medt(2,:),r_H,r_P,r_T,r);
        eff=cl./cd;
        % calculation of a and a' for all blade's parts
        aCD=ones(1,length(r));
        n=aCD;apCD=aCD;n2=aCD;
        for i=1:length(r)
            aCD(i)=ad(i);
            apCD(i)=apd(i);
            aold=1;
            apold=1;
            n(i)=0;
            n2(i)=0;
            while abs(apCD(i)-apold)>tol 
                n(i)=n(i)+1;
                apold=apCD(i);
                while abs(aCD(i)-aold)>tol
                    aold=aCD(i);
                    n2(i)=n2(i)+1;
                    phi(i)=atan((1-aCD(i))/((1+apCD(i))*xCDj(i)));
                    w(i)=(1-aCD(i))*vCD(j)/sin(phi(i));
                    f=Nb*(R-r(i)*R)/(2*r(i)*R*sin(phi(i)));
                    F_CD(i)=2/pi*acos(exp(-f));
                    betainf(i)=phi(i)*180/pi-betac(i);
                    aCD(i)=sol(i)*w(i)*(interp1(angle,cl(:,i),betainf(i))/tan(phi(i))+interp1(angle,cd(:,i),betainf(i)))/(4*vCD(j)*F_CD(i));
                    aCD(i)=e*aCD(i)+(1-e)*aold;
                end
                apCD(i)=sol(i)*w(i)*(interp1(angle,cl(:,i),betainf(i))-interp1(angle,cd(:,i),betainf(i))/tan(phi(i)))/(4*omega_d*r(i)*R*F_CD(i));
                apCD(i)=e*apCD(i)+(1-e)*apold;
            end
        end
        % bool definition for stall condition
        bool=zeros(1,length(r));
        for k=1:length(r)
            bool(k)=(interp1(angle,eff(:,k),betainf(k))<0.95*max(eff(:,k))&&betainf(k)>ALFA(1));
        end
        if sum(bool)<0.2*length(r)
            % power and Cp calculation 
            cc=ones(length(r),1);
            for i=1:length(r)
                cc(i)=cordalin(i)*(w(i)^2)*(interp1(angle,cl(:,i),betainf(i))*sin(phi(i))-interp1(angle,cd(:,i),betainf(i))*cos(phi(i)))*r(i)*dr*(R^2);
            end
            CpCD(j)=Nb*omega_d*sum(cc)/(pi*R^2*vCD(j)^3);
            PowerCD(j)=(ro_aria*Nb*omega_d*sum(cc)/2)/1000; % kW
            cc2=ones(length(r),1);
            dx2=dr*omega_d*R/vCD(j);
            lambda=omega_d*R/vCD(j);
            for i=1:length(r)
                cc2(i)=F_CD(i)*apCD(i)*(1-aCD(i))*(xCDj(i)^3)*dx2;
            end
            CpCD2(j)=8*sum(cc2)/(lambda^2);
            PowerCD2(j)=(4*ro_aria*pi*(R^2)*(vCD(j)^3)*sum(cc2)/(lambda^2))/1000;% kW
         else
            betac=betac+daB; %for every iteration, i will add the control pitch in order to avoid stall condition.
         end
    end
end
%% Range D-E
betac2=betac;
RatedPower=PowerCD(end);
tolPower=1e-2;%l'ultimo valore nel tratto CD è la massima potenza estraibile prima del pitch to feather
for j=1:length(vDE)
    xDEj = (omega_d*(r*R)/vDE(j));% non faccio xDE(j) perchè non serve salvarsi matrice, sovrascrive j-esimo vettore
    x_h = xDEj(r==r_H);
    x_p = xDEj(r==r_P);
    x_t = xDEj(r==r_T);
    PowerDE(j)=1;
    while abs((PowerDE(j)-RatedPower)/PowerDE(j))>tolPower
        %% d-e: H.P.T
        %hub
        ah=ad(r==r_H);
        aph=apd(r==r_H);
        aold=1;
        apold=1;
        nh=0;
        while abs(aph-apold)>tol 
            nh=nh+1;
            apold=aph;
            while abs(ah-aold)>tol
                aold=ah;
                phi=atan((1-ah)/((1+aph)*x_h));
                w=(1-ah)*vDE(j)/sin(phi);
                f=Nb*(R-r_H*R)/(2*r_H*R*sin(phi));
                F=2/pi*acos(exp(-f));
                Reh=w*cordalin(xDEj==x_h)/ni;
                    if Reh<=1e6
                        RE=Reh;
                    else
                        RE=1e6;
                    end
                    if RE>=5e5
                        RE=RE;
                    else
                        RE=5e5;
                    end
                medh=Re_medio(CLhub1,CLhub2,CDhub1,CDhub2,RE,1e6,5e5);
                betainf=phi*180/pi-betac2(r==r_H);
                ah=sol(xDEj==x_h)*w*(interp1(angle,medh(1,:),betainf)/tan(phi)+interp1(angle,medh(2,:),betainf))/(4*vDE(j)*F);
                ah=e*ah+(1-e)*aold;
            end
            aph=sol(xDEj==x_h)*w*(interp1(angle,medh(1,:),betainf)-interp1(angle,medh(2,:),betainf)/tan(phi))/(4*omega_d*r_H*R*F);
            aph=e*aph+(1-e)*apold;
        end 
        %primary
        ap=ad(r==r_P);
        app=apd(r==r_P);
        aold=1;
        apold=1;
        np=0;
        while abs(app-apold)>tol 
            np=np+1;
            apold=app;
            while abs(ap-aold)>tol
                aold=ap;
                phi=atan((1-ap)/((1+app)*x_p));
                w=(1-ap)*vDE(j)/sin(phi);
                f=Nb*(R-r_P*R)/(2*r_P*R*sin(phi));
                F=2/pi*acos(exp(-f));
                Rep=w*cordalin(xDEj==x_p)/ni;
                    if Rep<=1e6
                        RE=Rep;
                    else
                        RE=1e6;
                    end
                    if RE>=5e5
                        RE=RE;
                    else
                        RE=5e5;
                    end
                medp=Re_medio(CLpri1,CLpri2,CDpri1,CDpri2,RE,1e6,5e5);
                betainf=phi*180/pi-betac2(r==r_P);
                ap=sol(xDEj==x_p)*w*(interp1(angle,medp(1,:),betainf)/tan(phi)+interp1(angle,medp(2,:),betainf))/(4*vDE(j)*F);
                ap=e*ap+(1-e)*aold;
            end
            app=sol(xDEj==x_p)*w*(interp1(angle,medp(1,:),betainf)-interp1(angle,medp(2,:),betainf)/tan(phi))/(4*omega_d*r_P*R*F);
            app=e*app+(1-e)*apold;
        end
        %tip
        at=ad(r==r_T);
        apt=apd(r==r_T);
        aold=1;
        apold=1;
        nt=0;
        while abs(apt-apold)>tol 
            nt=nt+1;
            apold=apt;
            while abs(at-aold)>tol
                aold=at;
                phi=atan((1-at)/((1+apt)*x_t));
                w=(1-at)*vDE(j)/sin(phi);
                f=Nb*(R-r_T*R)/(2*r_T*R*sin(phi));
                F=2/pi*acos(exp(-f));
                Ret=w*cordalin(xDEj==x_t)/ni;
                    if Ret<=1e6
                        RE=Ret;
                    else
                        RE=1e6;
                    end
                    if RE>=5e5
                        RE=RE;
                    else
                        RE=5e5;
                    end
                medt=Re_medio(CLtip1,CLtip2,CDtip1,CDtip2,RE,1e6,5e5);
                betainf=phi*180/pi-betac2(r==r_T);
                at=sol(xDEj==x_t)*w*(interp1(angle,medt(1,:),betainf)/tan(phi)+interp1(angle,medt(2,:),betainf))/(4*vDE(j)*F);
                at=e*at+(1-e)*aold;
            end
            apt=sol(xDEj==x_t)*w*(interp1(angle,medt(1,:),betainf)-interp1(angle,medt(2,:),betainf)/tan(phi))/(4*omega_d*r_T*R*F);
            apt=e*apt+(1-e)*apold;
        end
        %% d-e: BLADESPAN
        % blade distribution of cl and cd
        [cl, cd]=blade(medh(1,:),medh(2,:),medp(1,:),medp(2,:),medt(1,:),medt(2,:),r_H,r_P,r_T,r);
        eff=cl./cd;
        % calculation of a and a' for all blade's parts
        aDE=ones(1,length(r));
        n=aDE;apDE=aDE;n2=aDE;
        for i=1:length(r)
            aDE(i)=ad(i);
            apDE(i)=apd(i);
            aold=1;
            apold=1;
            n(i)=0;
            n2(i)=0;
            while abs(apDE(i)-apold)>tol 
                n(i)=n(i)+1;
                apold=apDE(i);
                while abs(aDE(i)-aold)>tol
                    aold=aDE(i);
                    phi(i)=atan((1-aDE(i))/((1+apDE(i))*xDEj(i)));
                    w(i)=(1-aDE(i))*vDE(j)/sin(phi(i));
                    f=Nb*(R-r(i)*R)/(2*r(i)*R*sin(phi(i)));
                    F_DE(i)=2/pi*acos(exp(-f));
                    betainf(i)=phi(i)*180/pi-betac2(i);
                    aDE(i)=sol(i)*w(i)*(interp1(angle,cl(:,i),betainf(i))/tan(phi(i))+interp1(angle,cd(:,i),betainf(i)))/(4*vDE(j)*F_DE(i));
                    aDE(i)=e*aDE(i)+(1-e)*aold;
                end
                apDE(i)=sol(i)*w(i)*(interp1(angle,cl(:,i),betainf(i))-interp1(angle,cd(:,i),betainf(i))/tan(phi(i)))/(4*omega_d*r(i)*R*F_DE(i));
                apDE(i)=e*apDE(i)+(1-e)*apold;
            end
        end
             cc=ones(length(r),1);   
             for i=1:length(r)
                cc(i)=cordalin(i)*(w(i)^2)*(interp1(angle,cl(:,i),betainf(i))*sin(phi(i))-interp1(angle,cd(:,i),betainf(i))*cos(phi(i)))*r(i)*dr*(R^2);
            end
            CpDE(j)=Nb*omega_d*sum(cc)/(pi*R^2*vDE(j)^3);
            PowerDE(j)=(ro_aria*Nb*omega_d*sum(cc)/2)/1000; % kW
            
            cc2=ones(length(r),1);
            dx2=dr*omega_d*R/vDE(j);
            lambda=omega_d*R/vDE(j);
            for i=1:length(r)
                cc2(i)=F_DE(i)*apDE(i)*(1-aDE(i))*(xDEj(i)^3)*dx2;
            end
            CpDE2(j)=8*sum(cc2)/(lambda^2);
            PowerDE2(j)=(4*ro_aria*pi*(R^2)*(vDE(j)^3)*sum(cc2)/(lambda^2))/1000;% kW
            if abs((PowerDE(j)-RatedPower)/PowerDE(j))>tolPower
                betac2=betac2+daB;
            else
                betac2=betac2;
            end
    end
end 


Area = pi*R^2;
Pow_cub = [max(CpAC2)*ones(1,(length(CpAC2)+length(CpAC2)))].*(0.5*ro_aria*Area.*([vAC,vCD].^3))/1000;

SRO = max(PowerDE)/Area %kW/m2

%% PLOTS

figure()
plot(vAC,PowerAC./1000,'o',vCD,PowerCD,'o',vDE,PowerDE,'o');hold on;
plot([vAC,vCD],Pow_cub);hold off;
title('Power production');
legend('variable rpm','reaching Pmax','p.to feather','@ max Cp');
xlabel('Wind speed [m/s]'); 
ylabel('Power [kW]');

figure()
plot(vAC,CpAC,'o',vCD,CpCD,'o',vDE,CpDE,'o');
title('Coefficient of performance');    
xlabel('Wind speed [m/s]'); ylabel('Cp');
    
    














