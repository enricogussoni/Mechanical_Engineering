clear
close all
clc

%% Structure properties
rho=7850;
E = 2.06e11;
A = 3.657e-4;  % [m^2]
J = 11.260e-8;  % [m^4]
EJ = E*J ; % [Nm^2]
EA = E*A; % [N]
m=rho*A; % [kg/m]

fmax = 15; 
SC = 2;
Ommax = 2*pi*SC*fmax;

Lmax = sqrt(pi^2/Ommax * sqrt(EJ/m))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE06092011');


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
madd = 0.5;
Iadd = 3e-4;
idfM = [idf(2,:) idf(6,:)];
MM = [madd 0   0;
      0   madd 0;
      0    0  Iadd];
M(idfM,idfM)=M(idfM,idfM) + [MM zeros(3,3);zeros(3,3) MM];

% Internal springs
kx = 7.5e6;
ky = 3e6;

idof_kx1=idf(2,1);
K(idof_kx1,idof_kx1) = K(idof_kx1,idof_kx1)+kx;
idof_ky1=idf(2,2);
K(idof_ky1,idof_ky1) = K(idof_ky1,idof_ky1)+ky;
idof_kx6=idf(6,1);
K(idof_kx6,idof_kx6) = K(idof_kx6,idof_kx6)+kx;
idof_ky6=idf(6,2);
K(idof_ky6,idof_ky6) = K(idof_ky6,idof_ky6)+ky;

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

% A = [1/(2*om(1)) om(1)/2 ; 1/(2*om(2)) om(2)/2];
% h = [h1 h2 ...];
% ab = A^-1 * h;
alfah = 0.3;
betah = 6e-6;

R = alfah*M + betah*K;

RFF = R(1:ndof,1:ndof);

%% Frequency Response Function

% Force applied in E
freq = 0:0.1:400;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
idfEv = idf(4,2);
f0(idfEv)=1;

% Vertical displacement of D
for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
end
idfDv = idf(7,2);
FRF_F_Dv = xx(idfDv,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_F_Dv))
subplot(2,1,2)
plot(freq,angle(FRF_F_Dv))

clear ii
% Reaction in F
MFC = M(1:ndof,ndof+1:end);
RFC = R(1:ndof,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xf(:,ii) = A\f0;
    B = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    r(:,ii) = B*xf(:,ii);
end

idfFv=idf(8,2)-ndof;
FRF_FE_Fv = r(idfFv,:);

figure
subplot(211)
plot(f,abs(FRF_FE_Fv));
subplot(212)
plot(f,phase(FRF_FE_Fv));

%% Modal mass and modal stiffness

Phi = modes(:,1:3);
Mmodal_3f = Phi'*MFF*Phi; % for the first 3 mode shapes
Kmodal_3f = Phi'*KFF*Phi;
% freq_calc = sqrt(diag(Kmodal_3f)./diag(Mmodal_3f))./(2*pi);

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

for iHarm = 1:3
    for it = 1:length(td)
        x(:,it) = x(:,it)+modx(:,iHarm).*cos(om_harm(iHarm)*td(it)+phasex(:,iHarm));
        % A * cos(om*t + phi)
    end
end

xB = x(idfBo,:);
figure
plot(td,xB)