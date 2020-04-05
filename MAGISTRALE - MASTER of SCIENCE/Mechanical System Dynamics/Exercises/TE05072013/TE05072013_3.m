clear
close all
clc

%% Structure properties
E=2e11;
rho=7800;
A1 = 156e-4;
A2 = 156e-4;
A3 = 45.9e-4;
J1 = 92080e-8;
J2 = 92080e-8;
J3 = 5790e-8;
EA1=E*A1;
m1 = rho*A1;
EJ1=E*J1;
EA2=E*A2;
m2=rho*A2;
EJ2=E*J2;
EA3=E*A3;
m3=rho*A3;
EJ3=E*J3;

fmax = 30; 
SC = 2;
Ommax = 2*pi*SC*fmax;

Lmax1 = sqrt(pi^2/Ommax * sqrt(EJ1/m1))
Lmax2 = sqrt(pi^2/Ommax * sqrt(EJ2/m2))
Lmax3 = sqrt(pi^2/Ommax * sqrt(EJ1/m3))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE05072013_3');

% find(l>Lmax1)

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
scaleFactor = 2;
for ii = 1:3
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  
 

%% Damping Matrix
om = 2*pi*frqord;

alfah = 0.5;
betah = 5e-5;

R = alfah*M + betah*K;

% Contribution of localized dampers

RFF = R(1:ndof,1:ndof);

% Contribution of localized dampers

%% Frequency Response Function

% Force applied in C
freq = 0:0.05:30;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
idfCv = idf(10,2);
f0(idfCv)=1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
end
idfAv=idf(4,2);
FRF_F_A = xx(idfAv,:);
idfBv=idf(2,2);
FRF_F_B = xx(idfBv,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_F_A))
subplot(2,1,2)
plot(freq,angle(FRF_F_A))
figure
subplot(2,1,1)
plot(freq,abs(FRF_F_B))
subplot(2,1,2)
plot(freq,angle(FRF_F_B))

%% Modal mass and modal stiffness

Phi = modes(:,1:3);
Mmodal_3f = Phi'*MFF*Phi; % for the first 3 mode shapes
Kmodal_3f = Phi'*KFF*Phi;
% freq_calc = sqrt(diag(Kmodal_3f)./diag(Mmodal_3f))./(2*pi);

%% Response to a displacement
% Vertical displacement of B due to antiphase motion of O1 and O2

% Partitioning of "FC" Matrices
MFC = M(1:ndof,ndof+1:end);
RFC = R(1:ndof,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);

freq=linspace(0,fmax,1000);
Om=2*pi*freq;


xc = zeros(ndoc,1); 
idfO1v=idf(1,2)-ndof; 
idfO2v=idf(7,2)-ndof;

xc(idfO1v) = 1; 
xc(idfO2v) = -1;

clear ii A
for ii=1:length(freq)    
    A = -Om(ii)^2*MFF  +1i*Om(ii)*RFF  +KFF;
    B = -Om(ii)^2*MFC  +1i*Om(ii)*RFC  +KFC;
    xf(:,ii)= inv(A)*-B*xc;
end

FRF_O12_Bv = xf(idfBv,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_O12_Bv))
grid on
title('FRF: Vertical Displacement of A vs Imposed rotation clamp joint')
ylabel(['|Y_A/\theta| [m/rad]'])
subplot(2,1,2)
plot(freq,angle(FRF_O12_Bv))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on

%% Internal reactions on the right of A
dfu_dx = l/l(4);
idfA = idf(4,:);
idf5 = idf(5,:);
lambda = [cos(gamma(4)) sin(gamma(4)) 0;
          -sin(gamma(4)) cos(gamma(4)) 0;
          0               0            1];
x3 = lambda*xx(idfA,:);
x5 = lambda*xx(idf5,:);
load_ass_n3 = EA1*dfu_dx'*(x3(1,:)-x5(1,:));

figure
subplot(211)
plot(freq,abs(load_ass_n3))
subplot(212)
plot(freq,angle(load_ass_n3))

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

