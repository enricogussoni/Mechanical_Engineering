clear all
close all
clc
addpath('..\..\Codice\MeccFEM2\')


%% Structure properties


m1  = 33.50;    %[kg/m]
m2  = 10.30;    %[kg/m]

M1 = 20;        %[kg]
M2 = 50;        %[kg]

J1 = 0.04;      %[kgm^2]
J2 = 0.1;       %[kgm^2]

EA1 = 8.85e8;   %[N]
EA2 = 2.72e8;   %[N]
EJ1 = 3.10e6;   %[Nm^2]
EJ2 = 6.55e5;   %[Nm^2]

% Maximum length of finite element
fmax=200;
coef=2;
om=coef*(fmax*2*pi);
Lmax_fe_1 = sqrt(pi^2/om*sqrt(EJ1/m1))
Lmax_fe_2 = sqrt(pi^2/om*sqrt(EJ2/m2))


%% Load Structure Data

[file_i,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam]=MeccFEM2_loadstructure('FE5_1mod');

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
idof_n4 = idf(4,:);
idof_n5 = idf(5,:);
idof_n6 = idf(6,:);
MM1 = diag([M1 M1 J1]);
MM2 = diag([M2 M2 J2]);

M(idof_n4,idof_n4) = M(idof_n4,idof_n4)+MM1;
M(idof_n5,idof_n5) = M(idof_n5,idof_n5)+MM2;
M(idof_n6,idof_n6) = M(idof_n6,idof_n6)+MM1;

% Partitioning of "FF" Mass and Stiffness Matrices
MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);


%% Analysis of Natural Frequencies and Mode Shapes

[modes,Om2]=eig(MFF\KFF);
om2=diag(Om2);
frq=sqrt(om2)/(2*pi);

% Sort in ascending order frequencies and mode shapes
[frq,ordmode]=sort(frq);
modes = modes(:,ordmode);

nMod = sum(frq<fmax);

% Plot of mode shapes
fscale=1;
for i=1:3
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
fprintf('Natural frequency (1st to 3rd mode)  : [%i %i %i]\n', frq(1:3));
fprintf('Natural frequency with modal approach (1st to 3rd mode)  : [%i %i %i]\n', freq_calc);

%% Damping Matrix
alpha = 1;      %[s^-1]
beta = 1e-5;    %[s]
C=alpha*M+beta*K;

CFF=C(1:ndof,1:ndof);


%% Frequency Response Function

freq=[0:1:fmax].';
Om=2*pi*freq;
F0 = zeros(ndof,1);
idof_C_vert = idf(5,2);
idof_E_hori = idf(7,1);
F0(idof_C_vert) = 1;
for ii=1:length(freq)    
    xx(:,ii)=inv(-Om(ii)^2*MFF+sqrt(-1)*Om(ii)*CFF+KFF)*F0;
end

FRF_YC_F = xx(idof_C_vert,:);
FRF_XE_F = xx(idof_E_hori,:);

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
plot(freq,abs(FRF_XE_F))
grid on
title('FRF: Horizontal Displacement of E vs Force')
ylabel(['|X_E/F| [m/N]'])
subplot(212)
plot(freq,angle(FRF_XE_F))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on


%% Periodic External Torque

T = 0.1;                %[s]
Fmax = 1e5;             %[N]
fsample = 1000;         %[Hz]
dt = 1/fsample;         % Time step
time_vec = [dt:dt:T];   % Time vector
N = length(time_vec);

t_ramp1 = dt:dt:T/4;
t_ramp2 = T/4+dt:dt:3/4*T;
t_ramp3 = 3/4*T+dt:dt:T;

Force = [interp1([0,T/4],[0,Fmax],t_ramp1)  interp1([T/4, 3/4*T],[Fmax,-Fmax],t_ramp2)   interp1([3*T/4, T],[-Fmax, 0],t_ramp3)];
figure
plot(time_vec,Force);grid
xlabel('Time [s]');ylabel('Periodic Force [N]')
axis([0 time_vec(end) -Fmax*1.05 Fmax*1.05])

%Fourier Analysis with fft
Force_f = fft(Force,N);
nmax = N/2;
delta_freq = 1/T;                   %Frequency step 
freq_vec = delta_freq*[0:(nmax-1)]; %Frequency vector
Force_mod(1) = Force_f(1)/N;        %Torque modulus at 0 Hz 
Force_mod(2:nmax)=2/N*abs(Force_f(2:nmax));   %Torque modulus for frequencies > 0 Hz
Force_phase(1) = 0;
Force_phase(2:nmax) = angle(Force_f(2:nmax));          %Torque phase in radiant

figure
subplot(211)
bar(freq_vec,Force_mod,0.1);
set(gca,'Xlim',[0 500]);
subplot(212)
plot(freq_vec,Force_phase,'o')

%Frequency response
index_3_contributions = [2 4 6]; %even harmonic contributions and mean value are 0!
Ome_vec = 2*pi*freq_vec(index_3_contributions);
Force_f_3c = Force_mod(index_3_contributions).*exp(sqrt(-1)*Force_phase(index_3_contributions));

F0 = zeros(ndof,1);
idof_XA = idf(4,1);
for ii=1:length(Ome_vec)
    F0(idof_XA) = Force_f_3c(ii);
    xxp(:,ii)=inv(-Ome_vec(ii)^2*MFF+sqrt(-1)*Ome_vec(ii)*CFF+KFF)*F0;
end

% Time domain
t = [0:0.001:2];
xtot = zeros(ndof,length(t));

for jjj = 1:length(Ome_vec)
    for iii=1:length(t)
        xtot(:,iii)=xtot(:,iii)-Ome_vec(jjj)^2*abs(xxp(:,jjj)).*cos(Ome_vec(jjj)*t(iii)+angle(xxp(:,jjj)));
    end
end
idof_XB = idf(6,1);
figure
plot(t,xtot(idof_XB,:))
grid
xlabel('[s]')
ylabel('$\ddot{X}_B  [m/s^2]$','interpreter','latex')
title('Steady-state response of $\ddot{X}_B$','interpreter','latex')

%% 5
dfu_dx=1/l(10); % 10 è l'elemento finito di interesse (FE5_1.inp!!!!)
idof_n3 = idf(3,:);
idof_n11 = idf(11,:);
lambda = [cos(gamma(10)) sin(gamma(10)) 0; -sin(gamma(10)) cos(gamma(10)) 0; 0 0 1];
X3L = lambda*xx(idof_n3,:);
X11L = lambda*xx(idof_n11,:);
load_ass_n3 = EA2*dfu_dx*(X3L(1,:)-X11L(1,:)); %K*Dl

kL_ax = EA(2)/l(2)* [ 1 0 0 -1 0 0 % da MeccFEM2_beam
                0 0 0  0 0 0 
                0 0 0  0 0 0 
                -1 0 0  1 0 0 
                0 0 0  0 0 0 
                0 0 0  0 0 0 ] ; 

%oppure (considero la K come EA/Ltot
dfu_dx=1/(2*l(10)); % 10 è l'elemento finito di interesse (FE5_1.inp!!!!)

load_ass_n3_opzione = EA2*dfu_dx*X3L(1,:); 

figure
subplot(211)
plot(freq,abs(load_ass_n3),freq,abs(load_ass_n3_opzione),'r')
grid on
title('FRF: Assial load in D vs Force')
ylabel(['|Ax Load/F| [N/N]'])
subplot(212)
plot(freq,angle(load_ass_n3),freq,angle(load_ass_n3_opzione),'r')
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on


