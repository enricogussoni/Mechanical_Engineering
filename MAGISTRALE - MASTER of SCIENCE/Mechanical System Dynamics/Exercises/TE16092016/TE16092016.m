clear
close all
clc

%% Structure properties 
mo = 2;
EAo = 4e7;
EJo = 8e3;
mi = 1;
EAi = 1e5;
EJi = 5e3;

fmax = 3; 
SC = 1.5;
Ommax = 2*pi*SC*fmax;

Lmaxo = sqrt(pi^2/Ommax * sqrt(EJo/mo))
Lmaxi = sqrt(pi^2/Ommax * sqrt(EJi/mi))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('EG489');


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
scaleFactor = 5;
for ii = 1:3
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  

%% Modal mass and modal stiffness

Phi = modes(:,1:2);
modalMFF = Phi'*MFF*Phi; 
modalKFF = Phi'*KFF*Phi;

%% Damping Matrix
% om = 2*pi*frqord;

h1 = 0.01;
h2 = h1;

% A = [1/(2*om(1)) om(1)/2 ; 1/(2*om(2)) om(2)/2];
A = [1/(2*2*pi*frqord(1)) (2*pi*frqord(1)/2); 1/(2*2*pi*frqord(2)) (2*pi*frqord(2)/2)];
h = [h1 h2];
ab = A^-1 * h';
alfah = ab(1);
betah = ab(2);

R = alfah*M + betah*K;

RFF = R(1:ndof,1:ndof);

modalRFF = Phi'*RFF*Phi;

%% Frequency Response Function

% Force applied in Ax
freq = 0:0.01:3;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
idfAx = idf(10,1);
f0(idfAx)=1;
modalf0 = Phi'*f0;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
    
    modalA = -Om(ii)^2*modalMFF+sqrt(-1)*Om(ii)*modalRFF+modalKFF;
    modalxx(:,ii)=modalA\modalf0;
end

modalxx = Phi*modalxx;

idfBy = idf(12,2);
FRF_F_By = xx(idfBy,:);
idfBx = idf(12,1);
FRF_F_Bx = xx(idfBx,:);

modalFRF_F_By = modalxx(idfBy,:);
modalFRF_F_Bx = modalxx(idfBx,:);

figure
subplot(2,1,1)
title('FRF B vertical')
plot(freq,abs(FRF_F_By),freq,abs(modalFRF_F_By))
subplot(2,1,2)
plot(freq,angle(FRF_F_By),freq,angle(modalFRF_F_By))
legend('Usaual method','Modal coordinates')
figure
title('FRF B horizontal')
subplot(2,1,1)
plot(freq,abs(FRF_F_Bx),freq,abs(modalFRF_F_Bx))
subplot(2,1,2)
plot(freq,angle(FRF_F_Bx),freq,angle(modalFRF_F_Bx))
legend('Usaual method','Modal coordinates')

%% Response to a displacement

% Partitioning of "FC" Matrices
MFC = M(1:ndof,ndof+1:end);
RFC = R(1:ndof,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);

xc = zeros(ndoc,1); 

idof1x=idf(1,1)-ndof; 
idof19x=idf(19,1)-ndof; 
xc(idof1x,1) = 1;
xc(idof19x,1) = 1; 

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF  +1i*Om(ii)*RFF  +KFF;
    B = -Om(ii)^2*MFC  +1i*Om(ii)*RFC  +KFC;
    xf(:,ii)= inv(A)*-B*xc;
end

idofAx = idf(10,1);

FRF_eq_Ax = xf(idofAx,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_eq_Ax))
subplot(2,1,2)
plot(freq,angle(FRF_eq_Ax))

%% Distributed load
