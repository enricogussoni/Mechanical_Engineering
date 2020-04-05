clear all
close all
clc

%-------------------------------------------------------------------------
% Structure properties
%-------------------------------------------------------------------------

mp = 4;  %[kg/m]
mr = 1;  %[kg/m]
EAp = 4e7;%[N]
EAr = 4e3;%[N]
EJp = 1e5;%[Nm^2]
EJr = 1e3;%[Nm^2]

% control on finite elements length
sf = 2;
fmax = 20; %{Hz]
ommax = 2*pi*fmax*2;

lmaxp = sqrt(pi^2/ommax * sqrt(EJp/mp))
lmaxr = sqrt(pi^2/ommax * sqrt(EJr/mr))

%-------------------------------------------------------------------------
% Load Structure Data
%-------------------------------------------------------------------------
[file_i,xy,nnod,sizew,idb,ngdl,incidenze,l,gamma,m,EA,EJ,posiz,nbeam]=loadstructure(IMP04092015)

% Plot undeformed structure
figure
dis_stru(posiz,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')

% Check IDB and ndof
disp('Node  Degree-of-freedom')
fprintf('Node  %i:\t%i\t%i\t%i\n',[[1:size(idb,1)]' idb])
ngdl

disp('Press any key')
pause

%-------------------------------------------------------------------------
% Assembly of Mass and Stiffness Matrices
%-------------------------------------------------------------------------
[M,K]=assem(incidenze,l,m,EA,EJ,gamma,idb);

% Partitioning of "FF" Mass and Stiffness Matrices
MFF=M(1:ngdl,1:ngdl);
KFF=K(1:ngdl,1:ngdl);

%-------------------------------------------------------------------------
% Analysis of Natural Frequencies and Mode Shapes
%-------------------------------------------------------------------------

[modes,Om2]=eig(inv(MFF)*KFF);
om2=diag(Om2)
frq=sqrt(om2)/2/pi;

% Sort in ascending order frequencies and mode shapes
[frqord,ordmode]=sort(frq);


% Plot of mode shapes
fscala=2;
for ii=1:4
    mode=modes(:,ordmode(ii));
    figure
    diseg2(mode,fscala,incidenze,l,gamma,posiz,idb,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  
 
 



%---------------------------------------------------------------------
% Damping Matrix
%---------------------------------------------------------------------
alfah=0.2;
betah=1e-4;

R=alfah*M+betah*K;
RFF=R(1:ngdl,1:ngdl);

disp('Press any key')
pause

%---------------------------------------------------------------------
% Frequency Response Function
%---------------------------------------------------------------------
% Force applied in A
freq=[0:0.1:250].';
Om=2*pi*freq;
f0=zeros(ngdl,1);
ngdl_F_A = idb(11,2)
f0(ngdl_F_A)=1; % Vertical Force in A (node 11)

for ii=1:length(freq)    
    xx(:,ii)=inv(-Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF)*f0;
end
FRF_F_A = xx(ngdl_F_A,:);


figure
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