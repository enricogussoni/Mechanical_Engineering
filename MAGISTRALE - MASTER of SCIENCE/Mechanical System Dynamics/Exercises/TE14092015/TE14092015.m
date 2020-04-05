clear
close all
clc

%% Structure properties
m1 = 10;
EJ1 = 5e6;
EA1 = 2.6e8;
m2 = 50;
EJ2 = 1e8;
EA2 = 1.3e9;

fmax = 20;
SC = 1.5;
Ommax = 2*pi*SC*fmax;

Lmax1 = sqrt(pi^2/Ommax * sqrt(EJ1/m1))
Lmax2 = sqrt(pi^2/Ommax * sqrt(EJ2/m2))


%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE14092015');


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
scaleFactor = 2 ;
for ii = 1:3
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