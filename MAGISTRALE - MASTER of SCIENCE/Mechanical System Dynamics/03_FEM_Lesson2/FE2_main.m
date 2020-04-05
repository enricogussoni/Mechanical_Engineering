clear all
close all
clc

%% Structure properties

EJ = 1.34e4;    % [Nm2]
m = 9.75;       % [kg/m]
EA = 2.57e7;    % [N]
Mc = 10;         % [kg]
Jc = 1;          % [kgm^2]
kx = 2e6;       % [N/m]
ky = 3e6;       % [N/m]


% Maximum frequency
% fre_ef=(pi/L)^2*sqrt(EJ/m)/2/pi

% Maximum length of finite element
fmax=100;
coef=1.5;
om=coef*(fmax*2*pi);
Lmax_fe=sqrt(pi^2/om*sqrt(EJ/m));
disp('Maximum length')
disp(Lmax_fe)

%% Load Structure Data

[file_i,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam]=MeccFEM2_loadstructure('FE2');

% Plot undeformed structure
figure
MeccFEM2_plotStructure(position,xy)
xlabel('x [m]'); ylabel('y [m]')

% Check IDB and ndof
MeccFEM2_DoFsTable(idf)

%% Assembly of Mass and Stiffness Matrices

ndof_total = 3*nnod;

[M,K]=MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_total);

idof_mc = idf(5,:);
Mcc = [Mc 0 0; 0 Mc 0; 0 0 Jc];
M(idof_mc,idof_mc) = M(idof_mc,idof_mc)+Mcc; 

idof_kc=idf(5,1:2);
Kcc=[kx 0;
     0 ky];
K(idof_kc,idof_kc)=K(idof_kc,idof_kc)+Kcc;

%% Partitioning of "FF" Mass and Stiffness Matrices

MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);

%% Analysis of Natural Frequencies and Mode Shapes

[modes,Om2]=eig(inv(MFF)*KFF);
om2=diag(Om2);
frq=sqrt(om2)/2/pi;

% Sort in ascending order frequencies and mode shapes
[frqord,ordmode]=sort(frq);


% Plot of mode shapes
scale=2;
for ii=1:3
    mode=modes(:,ordmode(ii));
    figure
    MeccFEM2_plotDeformedStructure(mode,scale,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end


%% Damping Matrix

h1 = 0.01;
h2 = 0.015;

AA = [1/(2*2*pi*frqord(1)) (2*pi*frqord(1)/2); 1/(2*2*pi*frqord(2)) (2*pi*frqord(2)/2)];
XX = AA^-1*[h1;h2];
alpha = XX(1);
beta = XX(2);

C=alpha*M+beta*K;
CFF=C(1:ndof,1:ndof);

%% Frequency Response Function

% Force applied in node 4 (vertical)
freq=0:0.1:fmax;
Om=2*pi*freq;
f0=zeros(ndof,1);
ngdl_F_4 = idf(4,2); % the first index is the order of the node the second one is the direction
f0(ngdl_F_4)=1; % Vertical Force in node 4

for ii=1:length(freq)    
    xx(:,ii)=inv(-Om(ii)^2*MFF+sqrt(-1)*Om(ii)*CFF+KFF)*f0;
end
FRF_F_Y4 = xx(ngdl_F_4,:);
figure
subplot(211)
plot(freq,abs(FRF_F_Y4))
grid on
title('FRF: Vertical Displacement node 4 vs Vertical Force node 4')
ylabel('|Y_4/F| [m/N]')
subplot(212)
plot(freq,angle(FRF_F_Y4))
ylabel('\Psi [rad]')
xlabel('Freq [Hz]')
grid on

ngdl_F_Y5 = idf(5,2);
FRF_F_Y5 = xx(ngdl_F_Y5,:);
figure
subplot(211)
plot(freq,abs(FRF_F_Y5))
grid on
title('FRF: Vertical Displacement node 5 vs Vertical Force node 4')
ylabel('|Y_5/F| [m/N]')
subplot(212)
plot(freq,angle(FRF_F_Y5))
ylabel('\Psi [rad]')
xlabel('Freq [Hz]')
grid on


%% Frequency Response Function - Constraint Force

% Partitioning of "CF" Matrices

MCF=M(ndof+1:end,1:ndof);
CCF=C(ndof+1:end,1:ndof);
KCF=K(ndof+1:end,1:ndof);

for ii=1:length(freq)
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*CFF+KFF;
    xx(:,ii)= A\f0;
    
    B = -Om(ii)^2*MCF + 1i*Om(ii)*CCF + KCF;
    rr(:,ii)=B*xx(:,ii);
end

ngdv_F_1 = idf(1,3)-ndof; % the first index is the order of the node, the second the dof

FRF_F_C1 = rr(ngdv_F_1,:);

figure
subplot(211)
plot(freq,abs(FRF_F_C1))
grid on
title('FRF: Clamp joint Moment (node 1) vs Vertical Force node 4')
ylabel('|M_1/F| [m]')
subplot(212)
plot(freq,angle(FRF_F_C1))
ylabel('\Psi [rad]')
xlabel('Freq [Hz]')
grid on

