clear
close all
clc

%% Structure properties
me = 4;
EAe = 4e7;
EJe = 4e3;
mi = 1;
EAi = 1e5;
EJi = 1e3;

fmax = 10; 
SC = 1;
Ommax = 2*pi*SC*fmax;

Lmaxe = sqrt(pi^2/Ommax * sqrt(EJe/me))
Lmaxi = sqrt(pi^2/Ommax * sqrt(EJi/mi))

%% Load Structure Data

% [file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE12022016');
[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('Exam_20110620');

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
betah = 0;

R = alfah*M + betah*K;

% Contribution of localized dampers

RFF = R(1:ndof,1:ndof);

% Contribution of localized dampers

%% Frequency Response Function

% Force applied in 11 and orizzontally
freq = 0:0.01:10;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
f0(idf(11,1))= 1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+KFF;
    xx(:,ii)=A\f0;
end
FRF_F_xj = xx(idf(13,1),:);

figure
subplot(2,1,1)
semilogy(freq,abs(FRF_F_xj))
subplot(2,1,2)
plot(freq,angle(FRF_F_xj))
hold on

% Modal coordinates model
phi = modes(:,1:2);
Mmodal = phi'*MFF*phi;
Kmodal = phi'*KFF*phi;
mmodal = diag(Mmodal);
kmodal = diag(Kmodal);

freq_modal = sqrt(kmodal./mmodal)/2/pi;

for ii=1:length(freq)    
    A = -Om(ii)^2*Mmodal+Kmodal;
%     xxmodal(:,ii)=A\f0;
end
FRF_F_xjmodal = xxmodal(idf(13,1),:);

figure
subplot(2,1,1)
semilogy(freq,abs(FRF_F_xjmodal))
subplot(2,1,2)
plot(freq,angle(FRF_F_xjmodal))
hold off