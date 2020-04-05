clear
close all
clc

%% Structure properties
m = 59;
EA = 1.58e9;
EJ = 6.3e6;
M1 = 200;
k1 = 50e6;
k2 = 20e6;
L = 1.2;

fmax = 1500; 
SC = 2;
Ommax = 2*pi*SC*fmax;

Lmax = sqrt(pi^2/Ommax * sqrt(EJ/m))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('ExamSimulation2INP')


%% Plot undeformed structure

figure
MeccFEM2_plotStructure(position,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')


%% Check IDB and ndof

MeccFEM2_DoFsTable(idf)


%% Assembly of Mass and Stiffness Matrices

ndof_total = 3*nnod;
[M,K] = MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_total);

% Contribution of concentrated masses or springs
idof_ver_n6 = idf(6,2);
idof_ver_n7 = idf(7,2);
M(idof_ver_n6,idof_ver_n6)=M1;
M(idof_ver_n7,idof_ver_n7)=M1;

% Internal springs
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
for ii = 1:4
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  
 

%% Damping Matrix

% A = [1/(2*om(1)) om(1)/2 ; 1/(2*om(2)) om(2)/2];
% h = [h1 h2 ...];
% ab = A^-1 * h;
alfah = 0;
betah = 0.1e-5;

R = alfah*M + betah*K;

% Contribution of localized dampers

RFF = R(1:ndof,1:ndof);

% Contribution of localized dampers

%% Frequency Response Function

% Force applied in C (vertically)
fres = 1;
freq = 0:fres:fmax;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
idfCv = idf(3,2);
f0(idfCv) = 1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
end
idfBv = idf(2,2);
FRF_F_B = xx(idfBv,:);
idfCv = idf(3,2);
FRF_F_C = xx(idfCv,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_F_B))
subplot(2,1,2)
plot(freq,angle(FRF_F_B))
figure
subplot(2,1,1)
plot(freq,abs(FRF_F_C))
subplot(2,1,2)
plot(freq,angle(FRF_F_C))

%% Response to a displacement

% Partitioning of "FC" Matrices
MFC = M(1:ndof,ndof+1:end);
RFC = R(1:ndof,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);

freq=linspace(0,fmax,1000);
Om=2*pi*freq;

xc = zeros(size(M,1)-ndof,1); 

idfFv=idf(8,2)-ndof; 
                        
xc(idfFv,1) = 1; 

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF  +1i*Om(ii)*RFF  +KFF;
    B = -Om(ii)^2*MFC  +1i*Om(ii)*RFC  +KFC;
    xf(:,ii)= inv(A)*-B*xc;
end

idfDv = idf(4,2);

FRF_Dv_Fv = xf(idfDv,:);

figure
subplot(2,1,1)
semilogy(freq,abs(FRF_Dv_Fv))
subplot(2,1,2)
plot(freq,angle(FRF_Dv_Fv))

%% Reaction forces

% Partitioning of "CF" Matrices

MCF=M(ndof+1:end,1:ndof);
RCF=R(ndof+1:end,1:ndof);
KCF=K(ndof+1:end,1:ndof);

for ii=1:length(freq)
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)= A\f0;
    
    B = -Om(ii)^2*MCF + 1i*Om(ii)*RCF + KCF;
    rr(:,ii)=B*xx(:,ii);
end

idfFv = idf(8,2)-ndof; % the first index is the order of the node, the second the dof.
                       % -ndof is needed to take care only of CC rows
FRF_F_C = rr(idfFv,:);

figure
subplot(211)
plot(freq,abs(FRF_F_C))
subplot(212)
plot(freq,angle(FRF_F_C))

%% Spring reaction (B due to vertical displ. of F)
% Nodi interessati: 6 e 2 (2=B)
% Fmolla_k1 = k1*delta_L_k1
% delta_L_k1 = disp_2 + disp_6

idfBv=idf(2,2); 
idfk1v=idf(6,2);

xc(idfFv,1) = 1; 

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF  +1i*Om(ii)*RFF  +KFF;
    B = -Om(ii)^2*MFC  +1i*Om(ii)*RFC  +KFC;
    xf(:,ii)= inv(A)*-B*xc;
end

FRF_Bv_Fv = xf(idfBv,:);
FRF_k1v_Fv = xf(idfk1v,:);

delta_L_k1 = FRF_Bv_Fv + FRF_k1v_Fv;
Fspring = k1.*delta_L_k1;

figure
subplot(2,1,1)
semilogy(freq,abs(Fspring))
subplot(2,1,2)
plot(freq,angle(Fspring))

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