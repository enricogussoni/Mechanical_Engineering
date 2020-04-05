%% blade shape estimation
%TBN: far eseguire "project1_1_main.m" prima per importare i dati
close all
clc
%HUB NREL S818, PRIM NREL S827, TIP NREL S828
load NREL_S818_profile.txt;
load NREL_S828_profile.txt;
load NREL_S827_profile.txt;

xh=NREL_S818_profile(:,1);yh=NREL_S818_profile(:,2);
xp=NREL_S827_profile(:,1);yp=NREL_S827_profile(:,2);
xt=NREL_S828_profile(:,1);yt=NREL_S828_profile(:,2);

L=min([length(xh),length(xp),length(xt)]);
xh=interp1(1:length(xh),xh,1:L); yh=interp1(1:length(yh),yh,1:L);
xp=interp1(1:length(xp),xp,1:L); yp=interp1(1:length(yp),yp,1:L);
xt=interp1(1:length(xt),xt,1:L); yt=interp1(1:length(yt),yt,1:L);

for i = 1:2:length(r)
    corda_i = cordalin(r==r(i));
    if r(i)<=r_P
        xi = interp1([r_H;r_P],[xh;xp],r(i)) ; 
        yi = interp1([r_H;r_P],[yh;yp],r(i)) ; 
    else
        xi = interp1([r_P;r_T],[xp;xt],r(i)) ;
        yi = interp1([r_P;r_T],[yp;yt],r(i)) ; 
    end
    xi=xi*corda_i; 
    yi=yi*corda_i; 
    zi=r(i)*R*ones(1,L);
    bi=beta_ci(i)*pi/180;
    Rz_i = [cos(bi) -sin(bi) 0;
            sin(bi)  cos(bi) 0;
            0        0       1];
    xyzi = Rz_i*[xi;yi;zi]; xir=xyzi(1,:); yir=xyzi(2,:); zir=xyzi(3,:);
    xyz(i,:,:)=xyzi;
   
    figure(1)
        plot(xir,yir);axis equal;hold on; xlabel('[m]');ylabel('[m]');
    
    figure(2)
        plot3(xir,yir,zir);axis equal; hold on;xlabel('[m]');zlabel('[m]');        
end

% N=6;
% for j=1:length(r)
%     Pup(j,:) = polyfit(xyz(j,1,1:end/2),xyz(j,2,1:end/2),N); % necessari per CAD
%     Pdw(j,:) = polyfit(xyz(j,1,(end/2+1):end),xyz(j,2,(end/2+1):end),N);
%     fprintf('\ni=%f\nPu: %f *1m*(x/1m)^6 + (%f)*1m*(x/1m)^5 + (%f)*1m*(x/1m)^4 + (%f)*1m*(x/1m)^3 + (%f)*1m*(x/1m)^2 + (%f)*1m*(x/1m) + (%f)*1m\n',...
%         j,Pup(j,1),Pup(j,2),Pup(j,3),Pup(j,4),Pup(j,5),Pup(j,6),Pup(j,7));
%     fprintf('xup = %f m\n',xyz(j,1,1));
%     fprintf('\ni=%f\nPd: %f *1m*(x/1m)^6 + (%f)*1m*(x/1m)^5 + (%f)*1m*(x/1m)^4 + (%f)*1m*(x/1m)^3 + (%f)*1m*(x/1m)^2 + (%f)*1m*(x/1m) + (%f)*1m\n',...
%         j,Pdw(j,1),Pdw(j,2),Pdw(j,3),Pdw(j,4),Pdw(j,5),Pdw(j,6),Pdw(j,7));
%     fprintf('xdw = %f m\n',xyz(j,1,end));
% end


