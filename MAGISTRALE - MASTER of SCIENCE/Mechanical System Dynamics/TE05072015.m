%% TE05072013
clear
close all
clc

%% Structure properties
L = 30; %[m]
H = 5.625; %[m]
E = 2.06e11; %[N/m^2]
rho = 7800; %[km/m^3]

A1=15600e-4; %[m^2]
A2=15600e-4; %[m^2]
A3=45.9e-2;
J1=92080e-8; %[m^4]
J2=92080e-8; %[m^4]
J3=5790e-8;
EA1 = E*A1;
EA2 = E*A2;
EA3 = E*A3;
EJ1 = E*J1;
EJ2 = E*J2;
EJ3 = E*J3;

SC = 2;
fmax = 30; %[Hz]
Ommax = 2*pi*fmax*SC;
Lmax1 = sqrt(pi^2/Ommax * sqrt(EJ1/rho))
Lmax2 = Lmax1
Lmax3 = sqrt(pi^2/Ommax * sqrt(EJ2/rho))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE05072013_2');


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

% Sort in ascending order frequencies and mode shapes
[frqord,ordmode] = sort(frq);
modes = modes(:,ordmode);

% Plot of mode shapes
scaleFactor = 3;
for ii = 1:3
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  
 

%% Damping Matrix

alpha = 0.5;
beta = 5e-5;

R = alpha*M + beta*K;
RFF = R(1:ndof,1:ndof);


%% Frequency Response Function

% Force applied in 10, vertical
fres = 0.05;
freq = 0:fres:fmax;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
f0(idf(10,2))=1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
end
FRF_F_A = xx(idf(4,2),:);


figure
subplot(2,1,1)
plot(freq, abs(FRF_F_A))
subplot(2,1,2)
plot(freq, angle(FRF_F_A))

FRF_F_B = xx(idf(2,2),:);


figure
subplot(2,1,1)
plot(freq, abs(FRF_F_B))
subplot(2,1,2)
plot(freq, angle(FRF_F_B))

%% Constrain displacement
MFC = M(1:ndof,ndof+1:end);
RFC = R(1:ndof,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);

doc = size(M)-ndof;
xc = zeros(doc(1),1);
xc(idf(1,2)-ndof) = 1;
xc(idf(7,2)-ndof) = -1;

for kk=1:length(freq)
    B = -Om(kk)^2*MFF+sqrt(-1)*Om(kk)*RFF+KFF;
    C = -Om(kk)^2*MFC+sqrt(-1)*Om(kk)*RFC+KFC;
    xx(:,kk) = B^-1 * -C*xc;
end

idfBv=idf(2,2);
FRF_Yc_B = xx(idfBv,:);

figure
subplot(2, 1, 1)
plot(freq, abs(FRF_Yc_B))
subplot(2, 1, 2)
plot(freq, angle(FRF_Yc_B))

%% Reaction forces
MCF = M(ndof+1:end,1:ndof);
RCF = R(ndof+1:end,1:ndof);
KCF = K(ndof+1:end,1:ndof);

for cc=1:length(freq)
    xf(:,cc) = (-Om(cc)^2*MFF+sqrt(-1)*Om(cc)*RFF+KFF)^-1 * f0;
    fr(:,cc) = (-Om(cc)^2*MCF+sqrt(-1)*Om(cc)*RCF+KCF)*xf(:,cc);
end

idfFlecA = idf(4,3);
FRF_F_MA = xx(idfFlecA,:);

figure
subplot(2, 1, 1)
plot(freq, abs(FRF_F_MA))
subplot(2, 1, 2)
plot(freq, angle(FRF_F_MA))