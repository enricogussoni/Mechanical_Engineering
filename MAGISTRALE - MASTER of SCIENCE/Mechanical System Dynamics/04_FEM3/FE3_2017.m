clear 
clc

addpath('MeccFEM2')

%% Structure properties

EJ = 1.10e5;    % [Nm2]
m = 19.5;       % [kg/m]
EA = 5.15e8;    % [N]
mc = 20;         % [kg]
k1 = 2e6;       % [N/m]
k2 = 3e6;       % [N/m]
k3 = 2e6;       % [N/m]
c2 = 310;       % [Ns/m]


% Maximum frequency
% fre_ef=(pi/L)^2*sqrt(EJ/m)/2/pi

%% Maximum length of finite element
fmax=200;
coef=2;
om=coef*(fmax*2*pi);
Lmax_fe=sqrt(pi^2/om*sqrt(EJ/m));
disp('Maximum length')
disp(Lmax_fe)


%% Load Structure Data

[file_i,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam]=MeccFEM2_loadstructure('FE3');


% Plot undeformed structure
figure
MeccFEM2_plotStructure(position,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')

% Check IDB and ndof
MeccFEM2_DoFsTable(idf)
fprintf('Number of DoFs: %i\n',ndof)


%% Assembly of Mass and Stiffness Matrices

ndof_tot = 3*nnod;
[M,K]=MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_tot);

idof_m=idf(6,2);
M(idof_m,idof_m)=M(idof_m,idof_m)+mc;  % introduction of the lumped mass in the Mass Matrix

%% Contribution due to k1
idof_k1=idf(3,1);
K(idof_k1,idof_k1) = K(idof_k1,idof_k1)+k1;

%% Contribution due to k2
K_k2 = [ k2 -k2 ;
        -k2  k2 ];
    
idof_k2=idf([5,6],2);
K(idof_k2,idof_k2)=K(idof_k2,idof_k2)+K_k2;

%% Contribution due to k3
K_k3_local = [ k3 0 0 -k3 0 0;
               0  0 0  0  0 0;
               0  0 0  0  0 0;
              -k3 0 0  k3 0 0;
               0  0 0  0  0 0;
               0  0 0  0  0 0;];
alpha_k3 = 3/4*pi;
lambda_k3 = [ cos(alpha_k3) sin(alpha_k3) 0 ; 
             -sin(alpha_k3) cos(alpha_k3) 0 ;
              0      0      1 ];
          
Lambda_k3 = [lambda_k3  zeros(3,3) ; 
             zeros(3,3) lambda_k3  ];
         
K_k3_global = Lambda_k3'*K_k3_local*Lambda_k3;

idof_n2 = idf(2,:);
idof_n4 = idf(4,:);
idof_k3 = [idof_n2 idof_n4];

K(idof_k3,idof_k3)=K(idof_k3,idof_k3)+K_k3_global;

%% Partitioning of "FF" Mass and Stiffness Matrices
MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);


%% Analysis of Natural Frequencies and Mode Shapes


[modes,Om2]=eig(MFF\KFF);
om2=diag(Om2);
frq=sqrt(om2)/2/pi;

% Sort in ascending order frequencies and mode shapes
[frq,sort_index]=sort(frq);
modes=modes(:,sort_index);

% Plot of mode shapes
fscale=1;
for ii=1:3
    figure
    MeccFEM2_plotDeformedStructure(modes(:,ii),fscale,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frq(ii))])
end


%% Damping Matrix

% Structural isteretic damping
h1 = 0.01;
h2 = 0.015;
h3 = 0.018;

h=[h1;h2;h3];

frq_damp=frq(1:3);
coeffs = [1./(2*2*pi*frq_damp) (2*pi*frq_damp/2)];

alpha_beta = (coeffs'*coeffs)^-1*coeffs'*[h1;h2;h3];
% alpha_beta = coeffs\h;

alpha = alpha_beta(1);
beta = alpha_beta(2);

C=alpha*M+beta*K;

% Contribution due to c2
C_c2 = [ c2 -c2 ; 
        -c2 c2  ];
    
idof_c2 = idf([5,6],2);
C(idof_c2,idof_c2)=C(idof_c2,idof_c2)+C_c2;

% Extract the C_FF matrix
CFF=C(1:ndof,1:ndof);

% Plot of h=h(f)
h_computed=alpha./(2*2*pi*frq_damp) + beta.*(2*pi*frq_damp)/2;

ffs=linspace(0,fmax,1000);
h_fun=alpha./(2*2*pi*ffs) + beta.*(2*pi*ffs)/2;

figure
plot(ffs,h_fun,'b',frq_damp,h_computed,'ro')
title('Non-dimensional damping ratio')
xlabel('Freq. [Hz]')
grid on


%% Frequency Response Function

% Partitioning of "FC" Matrices
MFC = M(1:ndof,ndof+1:end);
CFC = C(1:ndof,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);

freq=linspace(0,fmax,1000);
Om=2*pi*freq;

xc = zeros(5,1); % 5 = doc

idof_theta=idf(1,3)-ndof; % riga della matrice -FC che corrisponde alla rotazone in 1
                          % (bisogna sottrarre al numero della riga dato da
                          % idf il numero di gradi di libertà.
xc(idof_theta,1) = 1; % Imposed rotation in clamp joint

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF  +1i*Om(ii)*CFF  +KFF;
    B = -Om(ii)^2*MFC  +1i*Om(ii)*CFC  +KFC;
    xf(:,ii)= inv(A)*-B*xc;
end

idof_YA = idf(4,2);

FRF_YA_Theta = xf(idof_YA,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_YA_Theta))
grid on
title('FRF: Vertical Displacement of A vs Imposed rotation clamp joint')
ylabel(['|Y_A/\theta| [m/rad]'])
subplot(2,1,2)
plot(freq,angle(FRF_YA_Theta))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on


idof_YB = idf(5,2);
FRF_YB_Theta = xf(idof_YB,:);

figure
subplot(211)
plot(freq,abs(FRF_YB_Theta))
grid on
title('FRF: Vertical Displacement of B vs Imposed rotation clamp joint')
ylabel(['|Y_B/\theta| [m/rad]'])
subplot(212)
plot(freq,angle(FRF_YB_Theta))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on


%% Frequency Response Function - Distributed force

f0 = zeros(ndof,1);
idof_y3 = idf(3,2);
idof_theta3 = idf(3,3);
idof_y4 = idf(4,2);
idof_theta4 = idf(4,3);
idof_y5 = idf(5,2);
idof_theta5 = idf(5,3);

l_fe = 0.5;
p0 = 1;

f0(idof_y3,1) = -p0*l_fe/2;
f0(idof_theta3,1) = -p0*l_fe^2/12;
f0(idof_y4,1) = -p0*l_fe/2 *2;
f0(idof_theta4,1) = 0; % =-p0*l_fe^2/12 + p0*l_fe^2/12
f0(idof_y5,1) = -p0*l_fe/2;
f0(idof_theta5,1) = p0*l_fe^2/12;

for ii=1:length(freq)   
    A = -Om(ii)^2*MFF  +1i*Om(ii)*CFF   +KFF;
    xf(:,ii)=A\f0;
end

FRF_YA_p = xf(idof_YA,:);

figure
subplot(211)
plot(freq,abs(FRF_YA_p))
grid on
title('FRF: Vertical Displacement of A vs Distributed force')
ylabel(['|Y_A/p_0| [m/N/m]'])
subplot(212)
plot(freq,angle(FRF_YA_p))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on

