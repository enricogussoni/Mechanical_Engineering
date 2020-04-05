clear
close all
clc

%% Structure properties
mA = 1000;
IA = 60;
E = 2.06e11;
rho = 7800;
A = 2.848e-3;
J = 1.9432e-5;
EJ = E*J;
EA = E*A;
m=rho*A;

fmax = 15; 
SC = 2;
Ommax = 2*pi*SC*fmax;

Lmax = sqrt(pi^2/Ommax * sqrt(EJ/m))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE22062015');


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
idfA = idf(6,:);
M(idfA,idfA)=M(idfA,idfA) + [mA  0  0;
                             0   mA 0;
                             0   0  IA];

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

% A = [1/(2*om(1)) om(1)/2 ; 1/(2*om(2)) om(2)/2];
% h = [h1 h2 ...];
% ab = A^-1 * h;
alfah = 0.4;
betah = 6e-5;

R = alfah*M + betah*K;

% Contribution of localized dampers

RFF = R(1:ndof,1:ndof);

%% Modal mass, stifness and damping
phi = modes(:,1:3);
Mmodal = phi'*MFF*phi;
Kmodal = phi'*KFF*phi;
Rmodal = phi'*RFF*phi;

%% Frequency Response Function

% Force applied in A
% Calculate the structure frequency response functions which relate the input
% force in node A (perpendicular to the axis of the beam connecting points A 
% and B) to the output vertical displacements evaluated in nodes A and B
freq = 0:0.01:15;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
idfAo = idf(6,1);
idfAv = idf(6,2);
f0(idfAo)=cos(pi/3);
f0(idfAv)=sin(pi/3);

modalf0 = phi'*f0;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
    
    modalA = -Om(ii)^2*Mmodal+sqrt(-1)*Om(ii)*Rmodal+Kmodal;
    modalxx(:,ii)=modalA\modalf0;
end

modalxx = phi*modalxx;              % passaggio dal significato ignoto

FRF_F_A = xx(idfAv,:);
idfBv = idf(9,2);
FRF_F_B = xx(idfBv,:);

modalFRF_F_A = modalxx(idfAv,:);
idfBv = idf(9,2);
modalFRF_F_B = modalxx(idfBv,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_F_A),freq,abs(modalFRF_F_A))
subplot(2,1,2)
plot(freq,angle(FRF_F_A),freq,angle(modalFRF_F_A))
figure
subplot(2,1,1)
plot(freq,abs(FRF_F_B),freq,abs(modalFRF_F_B))
subplot(2,1,2)
plot(freq,angle(FRF_F_B),freq,angle(modalFRF_F_B))
