clear
close all
clc

%% Structure properties
ml = 6600;
EAl = 1.6e11;
EJl = 1.0e13;
mt = 12000;
EAt = 3.6e11;
EJt = 1.2e12;

fmax = 5;
SC = 2;
Ommax = 2*pi*SC*fmax;

Lmaxl = sqrt(pi^2/Ommax * sqrt(EJl/ml))
Lmaxt = sqrt(pi^2/Ommax * sqrt(EJt/mt))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE10072015');


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
scaleFactor = ;
for ii = 1:4
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  
 

%% Damping Matrix

alfah = ;
betah = ;

R = alfah*M + betah*K;
RFF = R(1:ndof,1:ndof);


%% Frequency Response Function

% Force applied in A
freq = ;
Om = 2*pi*freq;

f0 = ;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
end
FRF_F_A = ;


figure
subplot(2,1,1)
subplot(2,1,2)