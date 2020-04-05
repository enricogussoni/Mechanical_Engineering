clear
close all
clc

%% Structure properties
fmax = 5; % [Hz]
Ommax = 2*pi*fmax;
ml= 6600;
mt = 12000;
EAl=1.6e11; % legs
EAt=3.6e11; % transverse
EJl=1e13;
EJt=1.2e12;

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
scaleFactor = 10;
for ii = 1:4
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  
 

%% Damping Matrix

alfah = 0.2;
betah = 2;

R = alfah*M + betah*K;
RFF = R(1:ndof,1:ndof);


%% Frequency Response Function

% Force applied in A
fres=0.01;
freq = 0:fres:fmax;
Om = 2*pi*freq;
idfA=idf(17,1);

f0 = zeros(1,ndof);
f0(idfA)=1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0';
end
FRF_F_A =xx(idfA,:);


figure
subplot(2,1,1)
plot(freq,abs(FRF_F_A));
subplot(2,1,2)
plot(freq,angle(FRF_F_A));

%..............................................
% Fourier decomposition of a periodic signal
% period of signal to be assigned manually
T=;

% values of signal at equally spaced discrete time values to be assigned
% manually
vect_F=[];

fftout=fft(vect_F);
N=length(vect_F);
df=1/T;
fmax=(N/2-1)*df;
vett_freq=0:df:fmax;
modf(1)=1/N*abs(fftout(1));
modf(2:N/2)=2/N*abs(fftout(2:N/2));
fasf(1:N/2)=angle(fftout(1:N/2));

figure
subplot 211;bar(vett_freq,modf);grid;
subplot 212;bar(vett_freq,fasf);grid



