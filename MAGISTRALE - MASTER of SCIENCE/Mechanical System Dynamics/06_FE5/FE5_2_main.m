clear all
close all
clc
addpath('..\..\Codice\MeccFEM2\')

%% Structure properties

m  = 59;    %[kg/m]

M1 = 200;        %[kg]

EA = 1.58e9;   %[N]
EJ = 6.3e6;   %[Nm^2]

k1 = 50e6;      %[N/m]
k2 = 20e6;      %[N/m]

% Maximum length of finite element
fmax=1500;
coef=2;
om=coef*(fmax*2*pi)
Lmax_fe = sqrt(pi^2/om*sqrt(EJ/m))

%% Load Structure Data

[file_i,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam]=MeccFEM2_loadstructure('FE5_2');

% Plot undeformed structure
figure
MeccFEM2_plotStructure(position,xy,m,EA,EJ);
xlabel('x [m]'); ylabel('y [m]')

% Check idf and ndof
MeccFEM2_DoFsTable(idf)
fprintf('Number of DoFs: %i\n',ndof)

%% Assembly of Mass and Stiffness Matrices

ndof_tot = 3*nnod;
[M,K]=MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_tot);


%Contribution due to lumped masses
idof_ver_n6 = idf(6,2);
idof_ver_n7 = idf(7,2);
M(idof_ver_n6,idof_ver_n6)=M1;
M(idof_ver_n7,idof_ver_n7)=M1;

K_k1 = [k1 -k1; -k1 k1];
K_k2 = [k2 -k2; -k2 k2];

idof_ver_n2 = idf(2,2);
idof_ver_n4 = idf(4,2);
idof_ver_n8 = idf(8,2);
idof_ver_n9 = idf(9,2);

ind_n2_6 = [idof_ver_n2 idof_ver_n6];
ind_n4_7 = [idof_ver_n4 idof_ver_n7];
ind_n6_8 = [idof_ver_n6 idof_ver_n8];
ind_n7_9 = [idof_ver_n7 idof_ver_n9];

K(ind_n2_6,ind_n2_6) = K(ind_n2_6,ind_n2_6) + K_k1;
K(ind_n4_7,ind_n4_7) = K(ind_n4_7,ind_n4_7) + K_k1;
K(ind_n6_8,ind_n6_8) = K(ind_n6_8,ind_n6_8) + K_k2;
K(ind_n7_9,ind_n7_9) = K(ind_n7_9,ind_n7_9) + K_k2;

%% Partitioning of "FF" Mass and Stiffness Matrices
MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);


%% Analysis of Natural Frequencies and Mode Shapes


[modes,Om2]=eig(inv(MFF)*KFF);
om2=diag(Om2);
frq=sqrt(om2)/2/pi;

% Sort in ascending order frequencies and mode shapes
[frq,ordmode]=sort(frq);
modes = modes(:,ordmode);

% Plot of mode shapes
fscale=1;
for i=1:4
    mode=modes(:,i);
    figure;
    MeccFEM2_plotDeformedStructure(mode,fscale,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(i) ': Freq [Hz]=' num2str(frq(i))])
end

%% Modal mass and modal stiffness
Phi = modes(:,1:3);
Mmodal_3f = Phi'*MFF*Phi; % for the first 3 mode shapes
Kmodal_3f = Phi'*KFF*Phi; % for the first 3 mode shapes

fprintf('Modal Mass      (1st to 3rd mode)  : [%d %d %d]\n', diag(Mmodal_3f));
fprintf('Modal Stiffness (1st to 3rd mode)  : [%d %d %d]\n', diag(Kmodal_3f));

freq_calc = sqrt(diag(Kmodal_3f)./diag(Mmodal_3f))./(2*pi);
fprintf('Natural frequency (1st to 3rd mode)  : [%d %d %d]\n', frq(1:3));
fprintf('Natural frequency with modal approach (1st to 3rd mode)  : [%d %d %d]\n', freq_calc);


%% Damping Matrix

alpha = 0;      %[s^-1]
beta = 0.1e-5;    %[s]
C=alpha*M+beta*K;

CFF=C(1:ndof,1:ndof);



%% Frequency Response Function

freq=[0:1:fmax].';
Om=2*pi*freq;
F0 = zeros(ndof,1);
idof_C_vert = idf(3,2);
idof_B_vert = idf(2,2);
F0(idof_C_vert) = 1;

for ii=1:length(freq) 
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*CFF+KFF;
    xx(:,ii)=inv(A)*F0;
end

FRF_YC_F = xx(idof_C_vert,:);
FRF_YB_F = xx(idof_B_vert,:);

figure
subplot(211)
plot(freq,abs(FRF_YC_F))
grid on
title('FRF: Vertical Displacement of C vs Force')
ylabel(['|Y_C/F| [m/N]'])
subplot(212)
plot(freq,angle(FRF_YC_F))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on

figure
subplot(211)
plot(freq,abs(FRF_YB_F))
grid on
title('FRF: Vertical Displacement of B vs Force')
ylabel(['|Y_B/F| [m/N]'])
subplot(212)
plot(freq,angle(FRF_YB_F))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on

%% Partitioning of "CF" Matrices

MCF = M(ndof+1:end,1:ndof);
CCF = C(ndof+1:end,1:ndof);
KCF = K(ndof+1:end,1:ndof);

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*CFF+KFF;
    B = -Om(ii)^2*MCF+sqrt(-1)*Om(ii)*CCF+KCF;
    rr(:,ii)=(B)*inv(A)*F0;
end

idof_vert_n8 = idf(8,2)-ndof;

FRF_VF_F = rr(idof_vert_n8,:);

figure
subplot(211)
semilogy(freq,abs(FRF_VF_F))
grid on
title('FRF: Vertical Constraint Force in F vs Force')
ylabel(['|V_F/F| [m/N]'])
subplot(212)
plot(freq,angle(FRF_VF_F))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on

%% Partitioning of "FC" Matrices

MFC = M(1:ndof,ndof+1:end);
CFC = C(1:ndof,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);

ndof_constraint = ndof_tot - ndof;
Xc = zeros(ndof_constraint,1);

Xc(idof_vert_n8,1) = 1; % Imposed vertical displacement in F (node 8)

for ii=1:length(freq)  
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*CFF+KFF;
    C = -Om(ii)^2*MFC+sqrt(-1)*Om(ii)*CFC+KFC;
    xx_imp(:,ii)=inv(A)*(-C*Xc);
end

idof_vert_n4 = idf(4,2);
FRF_YD_YF = xx_imp(idof_vert_n4,:);

figure
subplot(211)
plot(freq,abs(FRF_YD_YF))
grid on
title('FRF: Vertical Displacement of D vs Imposed vertical displacement in F')
ylabel(['|Y_D/Y_F| [m/m]'])
subplot(212)
plot(freq,angle(FRF_YD_YF))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on



