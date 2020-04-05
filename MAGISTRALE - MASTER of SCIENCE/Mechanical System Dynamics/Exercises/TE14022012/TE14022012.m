clear
close all
clc

%% Structure properties
m1= 50;
m2= 15;
m3= 20;
EA1=4e8;
EA2=2e8;
EA3=2e8;
EJ1=1e8;
EJ2=6e7;
EJ3=4.9e11;

fmax = 8; 
SC = 2;
Ommax = 2*pi*SC*fmax;

Lmax1 = sqrt(pi^2/Ommax * sqrt(EJ1/m1))
Lmax2 = sqrt(pi^2/Ommax * sqrt(EJ2/m2))
Lmax3 = sqrt(pi^2/Ommax * sqrt(EJ3/m3))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE14022012');


%% Plot undeformed structure

figure
MeccFEM2_plotStructure(position,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')


%% Check IDB and ndof

MeccFEM2_DoFsTable(idf)


%% Assembly of Mass and Stiffness Matrices

ndof_total = 3*nnod;
[M,K] = MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_total);

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
scaleFactor = 4;
for ii = 1:3
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  
 

%% Damping Matrix

h1 = 0.01;
h2 = 0.005;
h3 = h2;
h = [h1 h2 h3];

A = [1/(2*om(1)) om(1)/2 ; 1/(2*om(2)) om(2)/2;  1/(2*om(3)) om(3)/2];
ab = A\h';
alfah = ab(1);
betah = ab(2);

R = alfah*M + betah*K;

RFF = R(1:ndof,1:ndof);

%% Frequency Response Function

% Force applied in A and C
freq = 0:0.01:8;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
idfAv = idf(8,2);
idfDv = idf(16,2);
idfCv = idf(19,2);
f0(idfAv) = 1;
f0(idfCv) = -1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
end
FRF_F_A = xx(idfAv,:);
FRF_F_D = xx(idfDv,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_F_A))
subplot(2,1,2)
plot(freq,angle(FRF_F_A))
figure
subplot(2,1,1)
plot(freq,abs(FRF_F_D))
subplot(2,1,2)
plot(freq,angle(FRF_F_D))

%% Response to a displacement
% 
% % Partitioning of "FC" Matrices
% MFC = M(1:ndof,ndof+1:end);
% RFC = R(1:ndof,ndof+1:end);
% KFC = K(1:ndof,ndof+1:end);
% 
% freq=linspace(0,fmax,1000);
% Om=2*pi*freq;
% 
% xc = zeros(5,1); % 5 = doc
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

for iHarm = 1:3
    for it = 1:length(td)
        x(:,it) = x(:,it)+modx(:,iHarm).*cos(om_harm(iHarm)*td(it)+phasex(:,iHarm));
        % A * cos(om*t + phi)
    end
end

xB = x(idfBo,:);
figure
plot(td,xB)