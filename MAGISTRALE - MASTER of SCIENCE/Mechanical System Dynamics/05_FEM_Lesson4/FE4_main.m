clear all
close all
clc

addpath('..\Codice\MeccFEM2')

%% Structure properties

E = 2.06e11;    %[N/m^2]
rho = 7800;     %[kg/m^3]
k = 1e6;        %[N/m]

A1 = 20.1e-4;   %[m^2]
J1 = 869e-8;    %[m^4]

D2 = 20e-3;      %[m]
A2 = pi*D2^2/4;  %[m^2]
J2 = pi*D2^4/64; %[m^4]

EJ1 = E*J1;      %[Nm^2]
EA1 = E*A1;      %[N]
m1 = rho*A1;     %[kg/m]

EJ2 = E*J2;      %[Nm^2]
EA2 = E*A2;      %[N]
m2 = rho*A2;     %[kg/m]



% Maximum length of finite element
fmax=300;
coef=2;
om=coef*(fmax*2*pi);
Lmax_fe_1 = sqrt(pi^2/om*sqrt(EJ1/m1));
Lmax_fe_2 = sqrt(pi^2/om*sqrt(EJ2/m2));


%% Load Structure Data

[file_i,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam]=MeccFEM2_loadstructure('FE4');

% Plot undeformed structure
figure
MeccFEM2_plotStructure(position,xy);
xlabel('x [m]'); ylabel('y [m]')

% Check IDB and ndof
MeccFEM2_DoFsTable(idf)
fprintf('Number of DoFs: %i\n',ndof)


%% Assembly of Mass and Stiffness Matrices

ndof_tot = 3*nnod;
[M,K]=MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_tot);

%Contribution due to k
K(17,17) = K(17,17) + k;

%% Partitioning of "FF" Mass and Stiffness Matrices
MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);


%% Analysis of Natural Frequencies and Mode Shapes

[modes,Om2]=eig(MFF\KFF);
om2=diag(Om2);
frq=sqrt(om2)/(2*pi);

% Sort in ascending order frequencies and mode shapes
[frq,ordmode]=sort(frq);
modes = modes(:,ordmode);

% Plot of mode shapes
fscala=1;
for iHarm=1:5
    mode=modes(:,iHarm);
    figure
%     diseg2(mode,fscala,incidence,l,gamma,position,idb,xy);
    MeccFEM2_plotDeformedStructure(mode,fscala,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(iHarm) ': Freq [Hz]=' num2str(frq(iHarm))])
end


%% Damping Matrix

h1 = 0.01;
h2 = 0.01;
h3 = 0.01;
h4 = 0.01;

frq_damp=frq(1:4);
coeffs = [1./(2*2*pi*frq_damp) (2*pi*frq_damp/2)];
h=[h1;h2;h3;h4];

alpha_beta = coeffs\h;

alpha = alpha_beta(1);
beta = alpha_beta(2);

C=alpha*M+beta*K;

CFF=C(1:ndof,1:ndof);

% Plot of h=h(f)
h_f_spot = alpha./(2*2*pi*frq_damp)+ beta.*(2*pi*frq_damp./2);

ffs=[1:0.1:fmax]';
h_f_fun = alpha./(2*2*pi*ffs)+ beta.*(2*pi*ffs./2);

figure
plot(ffs,h_f_fun,'b',frq_damp,h_f_spot,'ro',frq_damp,h,'go')
title('Non-dimensional damping ratio')
xlabel('Freq. [Hz]')
grid on
ylim([0 0.03])


%% Frequency Response Function

freq=[0:0.1:fmax].';
Om=2*pi*freq;
F0 = zeros(ndof,1);
F0(idf(1,3)) = 1;

for iHarm=1:length(freq)    
    x(:,iHarm)=inv(-Om(iHarm)^2*MFF+sqrt(-1)*Om(iHarm)*CFF+KFF)*F0;
end
idof_XA = idf(4,1);
idof_XB = idf(3,1);
FRF_XA_C = x(idof_XA,:);
FRF_XB_C = x(idof_XB,:);

figure
subplot(211)
plot(freq,abs(FRF_XA_C))
grid on
title('FRF: Horizontal Displacement of A vs Torque')
ylabel(['|X_A/C| [m/Nm]'])
subplot(212)
plot(freq,angle(FRF_XA_C))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on

figure
subplot(211)
plot(freq,abs(FRF_XB_C))
grid on
title('FRF: Horizontal Displacement of B vs Torque')
ylabel(['|X_B/C| [m/Nm]'])
subplot(212)
plot(freq,angle(FRF_XB_C))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on


%% Case 2 with concentrated springs

clearvars -except alpha beta k EA2
kc = 1e5;   %  EA2/l !!!! [N/m]
l = 0.5*sqrt(2);
kc = EA2/l;

%% Load Structure Data

[file_i,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam]=MeccFEM2_loadstructure('FE4_spring');

% Plot undeformed structure
figure

MeccFEM2_plotStructure(position,xy);
% dis_stru(posiz,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')


%% Assembly of Mass and Stiffness Matrices

ndof_tot = 3*nnod;
[M,K]=MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_tot);

%Contribution due to k
K(17,17) = K(17,17)+k;

% Contribution due to kc
K_kc_local = [1 0 0 -1 0 0]'*kc*[1 0 0 -1 0 0];
g1 = pi/4;
g2 = 3/4*pi;
lambda_kc1 = [cos(g1) sin(g1) 0; -sin(g1) cos(g1) 0; 0 0 1];
Lambda_kc1 = [lambda_kc1 zeros(3,3); zeros(3,3) lambda_kc1];
K_kc1_global = Lambda_kc1'*K_kc_local*Lambda_kc1;

lambda_kc2 = [cos(g2) sin(g2) 0; -sin(g2) cos(g2) 0; 0 0 1];
Lambda_kc2 = [lambda_kc2 zeros(3,3); zeros(3,3) lambda_kc2];
K_kc2_global = Lambda_kc2'*K_kc_local*Lambda_kc2;

idof_n3 = idf(3,:);
idof_n5 = idf(5,:);
idof_kc1 = [idof_n3 idof_n5];
K(idof_kc1,idof_kc1)=K(idof_kc1,idof_kc1)+K_kc1_global;

idof_n6 = idf(6,:);
idof_n8 = idf(8,:);
idof_kc2 = [idof_n6 idof_n8];
K(idof_kc2,idof_kc2)=K(idof_kc2,idof_kc2)+K_kc2_global;


%% Partitioning of "FF" Mass and Stiffness Matrices

MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);
CFF=alpha*MFF+beta*KFF;

%% Analysis of Natural Frequencies and Mode Shapes

[modes,Om2]=eig(MFF\KFF);
om2=diag(Om2);
frq=sqrt(om2)/2/pi;

% Sort in ascending order frequencies and mode shapes
[frq,ordmode]=sort(frq);
modes = modes(:,ordmode);

% Plot of mode shapes
fscala=1;
for iHarm=1:3
    mode=modes(:,iHarm);
    figure
    MeccFEM2_plotDeformedStructure(mode,fscala,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(iHarm) ': Freq [Hz]=' num2str(frq(iHarm))])
end



%% Periodic External Torque

T = 3;                  %[s]
Cmax = 500;             %[Nm]
fsample = 1000;         %[Hz]
dt = 1/fsample;         % Time step
time_vec = [dt:dt:T];   % Time vector
N = length(time_vec);
Torque = zeros(1,N);
Torque(1:end/2) = Cmax;
Torque(end/2+1:end) = -Cmax;
figure
plot(time_vec,Torque);
grid
xlabel('Time [s]');ylabel('Periodic Torque [Nm]')
axis([0 T -Cmax*1.05 Cmax*1.05])

%Fourier Analysis with fft
Torque_f = fft(Torque);
nmax = N/2;
delta_freq = 1/T;                   %Frequency step 
freq_vec = delta_freq*[0:(nmax-1)]; %Frequency vector
Torque_mod(1) = Torque_f(1)/N;      %Torque modulus at 0 Hz 
Torque_mod(2:nmax)=2*abs(Torque_f(2:nmax))/N;   %Torque modulus for frequencies > 0 Hz
Torque_phase(1) = 0;
Torque_phase(2:nmax) = angle(Torque_f(2:nmax));

figure
h1 = subplot(2,1,1);
stem(freq_vec,Torque_mod);
set(gca,'Xlim',[0 500]);
h2 = subplot(212);
plot(freq_vec,Torque_phase,'o')
linkaxes([h1,h2],'x')

%Frequency response
index_3_contributions = [2 4 6]; %even harmonic contributions and mean value are 0!
Omega_3_conributions = 2*pi*freq_vec(index_3_contributions);
Torque_f_3contributions = Torque_mod(index_3_contributions).*exp(1i*Torque_phase(index_3_contributions));

F0 = zeros(ndof,1);

for iHarm=1:3
    F0(idf(1,3)) = Torque_f_3contributions(iHarm);
    x(:,iHarm) = inv(-Omega_3_conributions(iHarm)^2 * MFF + 1i * Omega_3_conributions(iHarm) *CFF + KFF)*F0;
end

modx = abs(x);
phasex = angle(x);

% Time domain

t = [0:0.01:10];
xtot = zeros(ndof,length(t));

for iHarm = 1:3
    for it=1:length(t)
        xtot(:,it) = xtot(:,it) +...
                      modx(:,iHarm).*cos(Omega_3_conributions(iHarm)*t(it)+phasex(:,iHarm));                       % A * cos(om*t + phi)
                      ... real(modx(:,iHarm).* exp(1i * (Omega_3_conributions(iHarm)*t(it)+phasex(:,iHarm))));     % RE( A * e^i(om*t + phi))
    end
end

idof_XA = idf(4,1);
XA = xtot(idof_XA,:);
figure
plot(t,XA)
grid
xlabel('[s]')
ylabel('X_A[m]')
title('Steady-state response of X_A')


