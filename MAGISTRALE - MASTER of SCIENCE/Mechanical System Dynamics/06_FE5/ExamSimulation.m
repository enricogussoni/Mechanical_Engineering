%% Exam Simulation

clear
close all
clc

%% Structure properties
m1 = 33.50;
EA1 = 8.85e8;
EJ1 = 3.10e6;
m2 = 10.30;
EA2 = 2.72e8;
EJ2 = 6.55e5;

fmin = 0;
fmax = 200;
sc = 2;
Ommax = sc*2*pi*fmax;

Lmax1 = sqrt(pi^2/Ommax * sqrt(EJ1/m1))
Lmax2 = sqrt(pi^2/Ommax * sqrt(EJ2/m2))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('ExamSimulationINP');


%% Plot undeformed structure

figure
MeccFEM2_plotStructure(position,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')


%% Check IDB and ndof

MeccFEM2_DoFsTable(idf)


%% Assembly of Mass and Stiffness Matrices

ndof_total = 3*nnod;
[M,K] = MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_total);

% Local masses contribution
M1 = 20;
I1 = 0.04;
M2 = 50;
I2 = 0.01;

idfA = idf(4,:);
idfC = idf(5,:);
idfB = idf(6,:);

M(idfA,idfA)= M(idfA,idfA)+[M1 0 0;
               0 M1 0;
               0  0 I1];
M(idfC,idfC)= M(idfC,idfC)+[M2 0 0;
               0 M2 0;
               0  0 I2];
M(idfB,idfB)= M(idfB,idfB)+[M1 0 0;
               0 M1 0;
               0  0 I1];

%% Partitioning of "FF" Mass and Stiffness Matrices
MFF = M(1:ndof,1:ndof);
KFF = K(1:ndof,1:ndof);


%% Analysis of Natural Frequencies and Mode Shapes

[modes, Om2] = eig(MFF\KFF);
om2 = diag(Om2);
frq = sqrt(om2)/2/pi;

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
 
% Modal mass and stiffness

Kstar = modes'*KFF*modes;
Mstar = modes'*MFF*modes;

Dm=diag(Mstar);
disp('First 3 modal mass:')
Dm(1:3)

Dk=diag(Kstar);
disp('First 3 modal stiffness:')
Dk(1:3)

%% Damping Matrix

alfah = 1;
betah = 1e-5;

R = alfah*M + betah*K;
RFF = R(1:ndof,1:ndof);


%% Frequency Response Function

% Vertical force in C
fres = 1;
freq = 0:fres:fmax;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
idfCv = idfC(2);
f0(idfCv) = 1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
end
FRF_F_C = xx(idfCv,:);

idfEo = idf(7,1);
FRF_F_E = xx(idfEo,:);

figure(1)
subplot(2,1,1)
plot(freq,abs(FRF_F_C));
subplot(2,1,2)
plot(freq,angle(FRF_F_C));
figure(2)
subplot(2,1,1)
plot(freq,abs(FRF_F_E));
subplot(2,1,2)
plot(freq,angle(FRF_F_E));

%% Time response
T = 0.1;
Fmax = 1e5; %[N]
n_harmonics = 3;
fsample=1000;
dt = 1/fsample;
t = dt:dt:T;
N = length(t);

% F1 = t(1:N/4)*Fmax/(N/4);
% F2 = -t(N/4+1:3*N/4)*Fmax/(N/4)+F1(end)*2;
% F3 = t(3*N/4+1:end)*Fmax/(N/4)+F2(end)*2;
% 
% figure
% plot(t(1:N/4),F1)
% hold on
% plot(t(N/4+1:3*N/4),F2)
% plot(t(3*N/4+1:end),F3)

t_ramp1 = dt:dt:T/4;
t_ramp2 = T/4+dt:dt:3/4*T;
t_ramp3 = 3/4*T+dt:dt:T;

Force = [interp1([0,T/4],[0,Fmax],t_ramp1)  interp1([T/4, 3/4*T],[Fmax,-Fmax],t_ramp2)   interp1([3*T/4, T],[-Fmax, 0],t_ramp3)];
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

% Frequency response (1st, 2nd and 3rd harmonics)
index_harm = [2 4 6]; %even harmonic contributions and mean value are 0!
freq_harm = freq(index_harm);
om_harm = 2*pi*freq_harm;
Force_harm = mod_Force(index_harm);

F0 = zeros(ndof,1);
idfBo = idf(6,1);

% superposition of the effects of the 3 harmonics
for iHarm = 1:3
    F0(idfBo) = Force_harm(iHarm);
    x(:,iHarm) = inv(-om_harm(iHarm)^2*MFF +1i*om_harm(iHarm)*RFF +KFF)*F0; % didplacement
    xdotdot(:,iHarm) = -om_harm(iHarm)^2.*x(:,iHarm);
end

modxdotdot = abs(xdotdot);
phasexdotdot = angle(xdotdot);

% Time domain
td = [0:0.001:2];
xdotdottot = zeros(ndof,length(td));

for iHarm = 1:3
    for it = 1:length(td)
        xdotdottot(:,it) = xdotdottot(:,it)+modxdotdot(:,iHarm).*cos(om_harm(iHarm)*td(it)+phasexdotdot(:,iHarm));
        % A * cos(om*t + phi)
    end
end

xBdotdot = xdotdottot(idfBo,:);
figure
plot(td,xBdotdot)

%% Frequency Response Function of the axial load of the upper end of O3D 
%  beam due to a horizontal force applied in C.














