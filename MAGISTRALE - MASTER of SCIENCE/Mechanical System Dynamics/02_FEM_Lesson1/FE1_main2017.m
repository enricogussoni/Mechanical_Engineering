clear
close all
clc

%% Structure properties
% parte indipendente dal resto che serve a controllare che nessun elemento
% vada in risonanza -> fissa la lunghezza massima dell'elemento
% SLIDE 3/22

EJ = 1.34e4;    % [Nm2]
m = 9.75;       % [kg/m]
EA = 2.57e7;    % [N}
kx = 2e6;       % [N/m]
ky = 3e6;       % [N/m]
Mc = 10;        % [kg]
J  = 1;         % [kgm^2]

% Maximum length of finite element
fmax=100;       % [Hz]
coef=1.5; % coefficiente di sicurezza
om=coef*(fmax*2*pi);
Lmax_fe=sqrt(pi^2/om*sqrt(EJ/m));

disp('Maximum length:')
disp(Lmax_fe) % informazione a schermo


%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam]=MeccFEM2_loadstructure('FE1');


%% Plot undeformed structure

figure
MeccFEM2_plotStructure(position,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')


%% Check IDB and ndof

MeccFEM2_DoFsTable(idf)


%% Assembly of Mass and Stiffness Matrices

ndof_total = 3*nnod;
[M,K]=MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_total);

Mx = [Mc 0 0;
      0 Mc 0;
      0  0 J];
  
M(10:12,10:12) = M(10:12,10:12)+Mx;

Kx= [kx 0  0;
      0 ky 0;
      0  0 0];
  
K(10:12,10:12) = K(10:12,10:12)+Kx;

%% Partitioning of "FF" Mass and Stiffness Matrices
% SLIDE 13/22

MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);


%% Analysis of Natural Frequencies and Mode Shapes

[modes,Om2]=eig(inv(MFF)*KFF); % diagonale con autovalori sulla diagonale
%[modes,Om2]=eig(MFF\KFF);
om2=diag(Om2);
frq=sqrt(om2)/2/pi;

% Sort in ascending order frequencies and mode shapes
[frqord,ordmode]=sort(frq);


% Plot of mode shapes
scaleFactor=2;        % a scelta -> ampiezza del disegno 
for ii=1:3
    mode=modes(:,ordmode(ii));
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  
 

%% Damping Matrix

alfah=0.2;
betah=1e-4;

R=alfah*M+betah*K;
RFF=R(1:ndof,1:ndof);


%% Frequency Response Function

% Force applied in A (node 11)
freq=0:0.1:250;
Om=2*pi*freq;
f0=zeros(ndof,1);     % Force vector
ngdl_F_A = idf(11,2); % Vertical displacement in the 11th node
f0(ngdl_F_A)=1;       % Vertical force in A (node 11)

for ii=1:length(freq)    
    xx(:,ii)=inv(-Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF)*f0;
end
FRF_F_A = xx(ngdl_F_A,:);

figure(11)
subplot(211)
semilogy(freq,abs(FRF_F_A))
grid on
title('FRF: Vertical Displacement in A vs Vertical Force in A')
ylabel(['|Y_A/Fy_A|'])
subplot(212)
plot(freq,angle(FRF_F_A))
ylabel(['[rad]'])
xlabel('Freq [Hz]')
grid on

