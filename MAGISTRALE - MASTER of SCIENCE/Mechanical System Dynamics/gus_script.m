clear
close all
clc

%% Structure properties 
IA = ; % [m^4]
E = 2.06e11;  % [N/m^2] -> Pa (non MPA)
rho = 7850;% [kg/m^3]
A = ;  % [m^2]
J = ;  % [m^4]
EJ = ; % [Nm^2]
EA = ; % [N]
m=rho*A; % [kg/m]

fmax = 15; 
SC = 2;
Ommax = 2*pi*SC*fmax;

Lmax = sqrt(pi^2/Ommax * sqrt(EJ/m))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE22062015');


%% Plot undeformed structure

figure
MeccFEM2_plotStructure(position,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')


%% Check IDB and ndof

MeccFEM2_DoFsTable(idf)


%% Assembly of Mass and Stiffness Matrices

ndof_total = 3*nnod;
[M,K] = MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_total);
ndoc=size(M,1)-ndof;

% Contribution of concentrated masses or springs
% M(,)=M(,) + [madd 0   0;
%              0   madd 0;
%              0    0  Jadd];

% Internal springs
%1)
% idof_k1=idf(3,1);
% K(idof_k1,idof_k1) = K(idof_k1,idof_k1)+k1;
% 
%2)
% K_k2 = [ k2 -k2 ;
%         -k2  k2 ];
%     
% idof_k2=idf([5,6],2);
% K(idof_k2,idof_k2)=K(idof_k2,idof_k2)+K_k2;
%
%3)
% k= ;
% Kloc = [k 0 0 -k 0 0;
%         0 0 0 0  0 0;
%         0 0 0 0  0 0;
%         -k 0 0 k 0 0;
%         0 0 0 0  0 0;
%         0 0 0 0  0 0];
% gammak = ;
% lambda = [cos(gammak) sin(gammak)  0;
%           -sin(gammak) cos(gammak) 0;
%            0            0        1];
% Lambda = [lambda zeros(3,3);
%           zeros(3,3) lambda];
% 
% idof_n2 = idf(2,:);
% idof_n4 = idf(4,:);
% idof_k3 = [idof_n2 idof_n4];
%
% Kglob = Lambda'*Kloc*Lambda;
% K(idof_k3,idof_k3) = K(idof_k3,idof_k3) + Kglob;

%% Partitioning of "FF" Mass and Stiffness Matrices
MFF = M(1:ndof,1:ndof);
KFF = K(1:ndof,1:ndof);


%% Analysis of Natural Frequencies and Mode Shapes

[modes, Om2] = eig(MFF\KFF);
om2 = diag(Om2);
frq = sqrt(om2)/2/pi;
om = sqrt(om2);

% Sort in ascending order frequencies and mode shapes
[frqord,ordmode] = sort(frq);
modes = modes(:,ordmode);

% Plot of mode shapes
scaleFactor = ;
for ii = 1:4
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  

%% Modal mass and modal stiffness

Phi = modes(:,1:3);
modalMFF_3f = Phi'*MFF*Phi; % for the first 3 mode shapes
modalKFF_3f = Phi'*KFF*Phi;
% freq_calc = sqrt(diag(Kmodal_3f)./diag(Mmodal_3f))./(2*pi);

%% Damping Matrix
om = 2*pi*frqord;

% A = [1/(2*om(1)) om(1)/2 ; 1/(2*om(2)) om(2)/2];
% h = [h1 h2 ...];
% ab = A^-1 * h;
alfah = ;
betah = ;

R = alfah*M + betah*K;

% Contribution of localized dampers
% rc = ; [Ns/m]

RFF = R(1:ndof,1:ndof);

modalRFF = phi'*RFF*phi;

%% Frequency Response Function

% Force applied in A
freq = ;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
modalf0 = phi'*f0;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
    
    modalA = -Om(ii)^2*modalMFF+sqrt(-1)*Om(ii)*modalRFF+modalKFF;
    modalxx(:,ii)=modalA\modalf0;
end

modalxx = phi*modalxx;

FRF_F_A = xx(idfAv,:);
idfBv = idf(9,2);
FRF_F_B = xx(idfBv,:);

modalFRF_F_A = modalxx(idfAv,:);
idfBv = idf(9,2);
modalFRF_F_B = modalxx(idfBv,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_F_A),freq,abs(modalFRF_F_A))
subplot(2,1,2)
plot(freq,angle(FRF_F_A),freq,angle(modalFRF_F_A))
figure
subplot(2,1,1)
plot(freq,abs(FRF_F_B),freq,abs(modalFRF_F_B))
subplot(2,1,2)
plot(freq,angle(FRF_F_B),freq,angle(modalFRF_F_B))

%% Response to a displacement
% 
% % Partitioning of "FC" Matrices
% MFC = M(1:ndof,ndof+1:end);
% RFC = R(1:ndof,ndof+1:end);
% KFC = K(1:ndof,ndof+1:end);
% 
% % freq=linspace(0,fmax,1000);
% % Om=2*pi*freq;
% 
% xc = zeros(ndoc,1); 
% 
% idof_theta=idf(1,3)-ndof; % riga della matrice -FC che corrisponde alla rotazione in 1
%                           % (bisogna sottrarre al numero della riga dato da
%                           % idf il numero di gradi di libertà.
% xc(idof_theta,1) = 1; % Imposed rotation in clamp joint
% 
% for ii=1:length(freq)    
%     A = -Om(ii)^2*MFF  +1i*Om(ii)*RFF  +KFF;
%     B = -Om(ii)^2*MFC  +1i*Om(ii)*RFC  +KFC;
%     xf(:,ii)= inv(A)*-B*xc;
% end
% 
% idof_YA = idf(4,2);
% 
% FRF_YA_Theta = xf(idof_YA,:);
% 
% figure
% subplot(2,1,1)
% plot(freq,abs(FRF_YA_Theta))
% grid on
% title('FRF: Vertical Displacement of A vs Imposed rotation clamp joint')
% ylabel(['|Y_A/\theta| [m/rad]'])
% subplot(2,1,2)
% plot(freq,angle(FRF_YA_Theta))
% ylabel(['\Psi [rad]'])
% xlabel('Freq [Hz]')
% grid on

%% Reaction forces
% 
% % Partitioning of "CF" Matrices
% 
% MCF=M(ndof+1:end,1:ndof);
% RCF=R(ndof+1:end,1:ndof);
% KCF=K(ndof+1:end,1:ndof);
% 
% for ii=1:length(freq)
%     A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
%     xx(:,ii)= A\f0;
%     
%     B = -Om(ii)^2*MCF + 1i*Om(ii)*RCF + KCF;
%     rr(:,ii)=B*xx(:,ii);
% end

% (1,3) 
% ngdv_F_1 = idf(1,3)-ndof; % the first index is the order of the node, the second the dof.
%                           % -ndof is needed to take care only of CC rows
% 
% FRF_F_C1 = rr(ngdv_F_1,:);
% 
% figure
% subplot(211)
% plot(freq,abs(FRF_F_C1))
% grid on
% title('FRF: Clamp joint Moment (node 1) vs Vertical Force node 4')
% ylabel('|M_1/F| [m]')
% subplot(212)
% plot(freq,angle(FRF_F_C1))
% ylabel('\Psi [rad]')
% xlabel('Freq [Hz]')
% grid on

%% Time response

T = ;
Fmax = ; %[N]
n_harmonics = ;
fsample=1000;
dt = 1/fsample;
t = dt:dt:T;
N = length(t);

Force = ;
figure
plot(t,Force);grid
xlabel('Time [s]');ylabel('Periodic Force [N]')
axis([0 t(end) -Fmax*1.05 Fmax*1.05])

% Fourier analysis with fft
fForce = fft(Force);
nmax = N/2;
delta_f = 1/T;
freq = 0:delta_f:delta_f*(nmax-1);
mod_Force(1) = fForce(1)/N;
mod_Force(2:nmax) = 2*abs(fForce(2:nmax))/N;
phase_Force(1) = 0;
phase_Force(2:nmax) = angle(fForce(2:nmax));

figure
subplot(211)
bar(freq,mod_Force)
subplot(212)
plot(freq,phase_Force,'o')

% Frequency response 
index_harm = [???]; 
freq_harm = freq(index_harm);
om_harm = 2*pi*freq_harm;
Force_harm = mod_Force(index_harm);

F0 = zeros(ndof,1);
idfBo = idf(6,1);

% superposition of the effects of the 3 harmonics
for iHarm = 1:nharmonics
    F0(idfBo) = Force_harm(iHarm);
    x(:,iHarm) = inv(-om_harm(iHarm)^2*MFF +1i*om_harm(iHarm)*RFF +KFF)*F0; % didplacement
%     xdotdot(:,iHarm) = -om_harm(iHarm)^2.*x(:,iHarm);
end

modx = abs(xdotdot);
phasex = angle(xdotdot);

% Time domain
td = [0:0.001:2];
x = zeros(ndof,length(td));

clear iHarm
for iHarm = 1:3
    for it = 1:length(td)
        A = modx(:,iHarm);
        variation = cos(om_harm(iHarm)*td(it)+phasex(:,iHarm));
        x(:,it) = x(:,it)-A.*variation;  % A * cos(om*t + phi)      
    end
end

xB = x(idfBo,:);
figure
plot(td,xB)

%% Internal actions

dfu_dx=1/l(10); % 10 è l'elemento finito di interesse
idof_n3 = idf(3,:);
idof_n11 = idf(11,:);
lambda = [cos(gamma(10)) sin(gamma(10)) 0; -sin(gamma(10)) cos(gamma(10)) 0; 0 0 1];
X3L = lambda*xx(idof_n3,:);
X11L = lambda*xx(idof_n11-ndof,:);
load_ass_n3 = EA2*dfu_dx*(X3L(1,:)-X11L(1,:)); %K*Dl

% alternativa
kL_ax = EA(2)/l(2)* [ 1 0 0 -1 0 0 % da MeccFEM2_beam
                0 0 0  0 0 0 
                0 0 0  0 0 0 
                -1 0 0  1 0 0 
                0 0 0  0 0 0 
                0 0 0  0 0 0 ] ;
            
load_ass_n3 = kL_ax*X3L + kL_ax*X11L;
mod_load_ass_n3 = abs();
phase_load_ass_n3 = angle();

% oppure (considero la K come EA/Ltot
dfu_dx=1/(2*l(10)); % 10 è l'elemento finito di interesse

load_ass_n3_opzione = EA2*dfu_dx*X3L(1,:); 

figure
subplot(211)
plot(freq,abs(load_ass_n3),freq,abs(load_ass_n3_opzione),'r')
subplot(212)
plot(freq,angle(load_ass_n3),freq,angle(load_ass_n3_opzione),'r')

%% Distributed load
% 
% p = 9400;
% % g1 = gamma(5)
% % g2 = gamma(15)
% 
% g = gamma(6);
% po = -p*sin(g);
% pv = p*cos(g);
% 
% Floc = zeros(ndof,1);
% Fglo = zeros(ndof,1);
% 
% Floc(idb(6,:),1) = [po*l(5)/2 -pv*l(5)/2 -pv*l(5)^2/12];
% lambda = [cos(g) sin(g) 0; -sin(g) cos(g) 0; 0 0 1];
% Fglo(idb(6,:),1) = Floc(idb(6,:),1)'*lambda;
% 
% Floc(idb(7,:),1) = [po*(l(5)+l(6))/2, -pv*(l(5)+l(6))/2, pv*(l(5)^2-l(6)^2)/12];
% Fglo(idb(7,:),1) = Floc(idb(7,:),1)'*lambda;
% 
% Floc(idb(8,:),1) = [po*(l(6)+l(7))/2, -pv*(l(6)+l(7))/2, pv*(l(6)^2-l(7)^2)/12];
% Fglo(idb(8,:),1) = Floc(idb(8,:),1)'*lambda;
% 
% Floc(idb(9,:),1) = [po*(l(7)+l(8))/2, -pv*(l(7)+l(8))/2, pv*(l(7)^2-l(8)^2)/12];
% Fglo(idb(9,:),1) = Floc(idb(9,:),1)'*lambda;
% 
% Floc(idb(10,:),1) = [po*l(8)/2, -pv*l(8)/2, pv*l(8)^2/12];
% Fglo(idb(10,:),1) = Floc(idb(10,:),1)'*lambda;
% 
% g = gamma(17);
% 
% Floc(idb(6,:),1) = [po*l(14)/2, pv*l(14)/2, pv*l(14)^2/12];
% lambda = [cos(g) sin(g) 0; -sin(g) cos(g) 0; 0 0 1];
% Fglo(idb(6,:),1) = Fglo(idb(6,:),1)' + Floc(idb(6,:),1)'*lambda;
% 
% Floc(idb(15,:),1) = [po*(l(14)+l(15))/2, pv*(l(14)+l(15))/2, pv*(l(15)^2-l(14)^2)/12];
% Fglo(idb(15,:),1) = Floc(idb(15,:),1)'*lambda;
% 
% Floc(idb(16,:),1) = [po*(l(15)+l(16))/2, pv*(l(15)+l(16))/2, pv*(l(16)^2-l(15)^2)/12];
% Fglo(idb(16,:),1) = Floc(idb(16,:),1)'*lambda;
% 
% Floc(idb(17,:),1) = [po*(l(16)+l(17))/2, pv*(l(16)+l(17))/2, pv*(l(17)^2-l(16)^2)/12];
% Fglo(idb(17,:),1) = Floc(idb(17,:),1)'*lambda;
% 
% Floc(idb(18,:),1) = [po*l(17)/2, pv*l(17)/2, pv*l(17)^2/12];
% Fglo(idb(18,:),1) = Floc(idb(18,:),1)'*lambda;
% 
% xstat = inv(KFF)*(Fglo);
% vertdispC = xstat(idb(6,2))
% hordispC = xstat(idb(6,1))