clear
close all
clc

%% Structure properties
fmax = 300;
Ommax = 2*pi*fmax;
SC=2;
E = 2.06e11;    %[N/m^2]
rho = 7800;     %[kg/m^3]
k = 1e6;        %[N/m]
L=0.5;          %[m]

A1 = 20.1e-4;   %[m^2]
J1 = 869e-8;    %[m^4]

D2 = 20e-3;      %[m]
A2 = pi*D2^2/4;  %[m^2]
J2 = pi*D2^4/64; %[m^4]

EJ1 = E*J1;      %[Nm^2]
EA1 = E*A1;      %[N]
m1 = rho*A1;     %[kg/m]

EJ2 = E*J2;      %[Nm^2]
EA2 = E*A2;      %[N]
m2 = rho*A2;     %[kg/m]

Lmax1 = sqrt(pi^2/(2*Ommax) + sqrt(EJ1/rho));
Lmax2 = sqrt(pi^2/(2*Ommax) + sqrt(EJ2/rho));
Lmax = max(Lmax1,Lmax2)

%% Load Structure Data

% [file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('nodi');

% Modellizations of type 2 beams with springs
[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('nodi2springs');


%% Plot undeformed structure

figure
MeccFEM2_plotStructure(position,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')


%% Check IDB and ndof

MeccFEM2_DoFsTable(idf)


%% Assembly of Mass and Stiffness Matrices

ndof_total = 3*nnod;
[M,K] = MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_total);

%Contribution due to k
K(17,17) = K(17,17) + k;

%% Modellizations of type 2 beams with springs 

kadd = E*A2/(sqrt(2)*L);
gammaL = pi/4;
lambdaL = [cos(gammaL) sin(gammaL) 0;
          -sin(gammaL) cos(gammaL) 0;
          0 0 1];
LambdaL = [lambdaL zeros(3,3);
           zeros(3,3) lambdaL];
      
gammaR = 3*pi/4;
lambdaR = [cos(gammaR) sin(gammaR) 0;
          -sin(gammaR) cos(gammaR) 0;
          0 0 1];
LambdaR = [lambdaR zeros(3,3);
           zeros(3,3) lambdaR];
      
kL= [kadd  0 0 -kadd 0 0;
     0     0 0 0     0 0;
     0     0 0 0     0 0;
     -kadd 0 0  kadd 0 0;
     0     0 0 0     0 0;
     0     0 0 0     0 0];
 
kLG = LambdaL'*kL*LambdaL; % OK
 
idof_n3 = idf(3,:);
idof_n5 = idf(5,:);
idof_kL = [idof_n3 idof_n5];
K(idof_kL,idof_kL)=K(idof_kL,idof_kL)+kLG;

kR= [kadd  0 0 -kadd 0 0;
     0     0 0 0     0 0;
     0     0 0 0     0 0;
     -kadd 0 0  kadd 0 0;
     0     0 0 0     0 0;
     0     0 0 0     0 0];
 
kRG = LambdaR'*kR*LambdaR; % OK
 
idof_n6 = idf(6,:);
idof_n8 = idf(8,:);
idof_kR = [idof_n6 idof_n8];
K(idof_kR,idof_kR)=K(idof_kR,idof_kR)+kRG;

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
for ii = 1:3 %5
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  
 

%% Damping Matrix
h1 = 0.01;
h2 = h1;
h3 = h1;
h4 = h1;

h = [h1 h2 h3 h4];
% AA = [1/2/om(1) om(1)/2 ; 1/2/om(2) om(2)/2 ; 1/2/om(3) om(3)/2 ; 1/2/om(4) om(4)/2 ];
frq_damp=frq(1:4);
AA = [1./(2*2*pi*frq_damp) (2*pi*frq_damp/2)];
ab = AA\h';

alfah = ab(1);
betah = ab(2);

R = alfah*M + betah*K;
RFF = R(1:ndof,1:ndof);

% freq_damp = frq(1:length(h));
% h_damp = alfah./(2*pi*freq_damp)+ betah.*(2*pi*freq_damp)/2;
h_damp = alfah./(2*2*pi*frq_damp)+ betah.*(2*pi*frq_damp./2);
fxplot = [1:0.1:fmax]';
hxplot = alfah./(2*2*pi*fxplot) + betah.*(2*pi*fxplot./2);

% figure
% plot(fxplot,hxplot,'b',freq_damp,h_damp,'o');
% legend('Experimental damping', 'Extimated damping')
% xlabel('Frequency'); ylabel('Damping ratio')
figure
plot(fxplot,hxplot,frq_damp,h_damp,'o');
ylim([0 0.03])
% legend('Experimental damping', 'Extimated damping')
% xlabel('Frequency'); ylabel('Damping ratio')

%% Frequency Response Function

% Force applied in A
fres=0.1;
freq =0:fres:fmax ;
Om = 2*pi*freq;

f0 = zeros(1,ndof);
idf1or = idf(1,3);
f0(idf1or) = 1;

clear ii
clear xx
for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=(A^-1)*f0';
end
idfAor=idf(4,1);
FRF_F_A = xx(idfAor,:);
idfBor=idf(3,1);
FRF_F_B = xx(idfBor,:);

figure
subplot(2,2,1)
plot(freq,abs(FRF_F_A ))
xlabel('Frequency'); ylabel('|FRF_F_A |');
subplot(2,2,2)
plot(freq,abs(FRF_F_B ))
xlabel('Frequency'); ylabel('|FRF_F_B |');
subplot(2,2,3)
plot(freq,angle(FRF_F_A ))
xlabel('Frequency'); ylabel('<FRF_F_A ');
subplot(2,2,4)
plot(freq,angle(FRF_F_B ))
xlabel('Frequency'); ylabel('<FRF_F_B ');