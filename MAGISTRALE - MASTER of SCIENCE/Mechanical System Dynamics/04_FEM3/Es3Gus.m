clear all
close all

% Data
m = 19.5;     % Mass per unit length [kg]
EJ = 1.1e5;   % Flexural rigidity [Nm^2]
EA = 5.15e10; % Axial stiffness [N]
mc = 20;      % Lumped mass [kg]
k1 = 2e6;
k2 = 3e6;
k3 = 2e6;
c2 = 310;     % Damper [Ns/m]

% Maximum length of finite elements
fmax = 200;   
safety_coef = 2;
om=safety_coef*(fmax*2*pi);
Lmax_fe=sqrt(pi^2/om*sqrt(EJ/m));

disp('Maximum length:')
disp(Lmax_fe)

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam]=MeccFEM2_loadstructure('FE3gus');

% Undeformed structure
% figure
% MeccFEM2_plotStructure(position,l,gamma,xy)
% xlabel('x [m]'); ylabel('y [m]')

%% Check IDB and ndof
MeccFEM2_DoFsTable(idf)

% Assembling of elements from .inp file 
ndof_total = 3*nnod;
[M,K]=MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_total);

%% Assembling of lumped parameter elements to be done manually

% Contribution due to m
idof_m = idf(6,2);
M(idof_m,idof_m) = M(idof_m,idof_m) + mc;

% " " " k1
idof_k1=idf(3,1);
K(idof_k1,idof_k1) = K(idof_k1,idof_k1) + k1;

% " " " k2
K_k2 = [ k2 -k2;
        -k2  k2];

idof_k2 = idf([5,6],2);
K(idof_k2,idof_k2) = K(idof_k2,idof_k2)+K_k2;

% " " " k3
K_k3_local = [k3 0 0 -k3 0 0;
              0  0 0  0  0 0;
              0  0 0  0  0 0;
             -k3 0 0  k3 0 0;
              0  0 0  0  0 0;
              0  0 0  0  0 0];
    
alpha_k3 = 3/4*pi;
lambda_k3 = [cos(alpha_k3)  sin(alpha_k3) 0;
             -sin(alpha_k3) cos(alpha_k3) 0;
             0                 0          1]; 

Lambda_k3 = [lambda_k3     zeros(3,3);
             zeros(3,3)    lambda_k3];
             
K_k3_global = Lambda_k3'*K_k3_local*Lambda_k3;

idof_k3_2 = idf(2,:);
idof_k3_4 = idf(4,:);
idof_k3 = [idof_k3_2, idof_k3_4];

K(idof_k3,idof_k3) = K(idof_k3,idof_k3)+K_k3_global;

%% Partitioning of matrices 
MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);

%% Natural frequencies and modes of vibration
[modes,Om2]=eig(inv(MFF)*KFF);
%[modes,Om2]=eig(MFF\KFF);
nfreq=sqrt(diag(Om2))/2/pi;
[frqord,ind]=sort(nfreq);

% Plot of mode shapes
scaleFactor=2;        % a scelta -> ampiezza del disegno 
for ii=1:3
    mode=modes(:,ind(ii));
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  


%% alfa and beta values for structural damping to be assigned manually
Om = sqrt(Om2);

h1 = 0.01;
h2 = 0.015;
h3 = 0.018;

h = [h1 h2 h3]';
A = [1/(2*Om(1)) Om(1)/2;
     1/(2*Om(2)) Om(2)/2;
     1/(2*Om(3)) Om(3)/2];
 
x = A\h;

alfa=x(1);
beta=x(2);

C=alfa*M+beta*K;

% Contribution of c2 (like for k2)

C_c2 = [ c2 -c2;
        -c2  c2];

idof_c2 = idf([5,6],2);
C(idof_c2,idof_c2) = K(idof_c2,idof_k2)+C_c2;

CFF = C(1:ndof,1:ndof);

%% Number of modes to be plotted to be assigned manually
nmodes=3;
modes=ind(1:nmodes);

% vector of frequency to be defined
vett_f=linspace(0,fmax,1000);
omega=vett_f*2*pi;

% % Forcing vector 
% F0=zeros(1,length(x)); % no forces applied
% for k=1:length(vett_f)
%     A=-omega(k)^2*MFF+1i*omega(k)*CFF+KFF;
%     x=A\F0;
%     %outputs to be defined, script prepared for two outputs 
%     yg1=idf(4,2); % vertical displacement of A
%     yg2=idf(5,2); % " " " B
%     mod1(k)=abs(yg1);
%     fas1(k)=angle(yg1);
%     mod2(k)=abs(yg2);
%     fas2(k)=angle(yg2);
% end

% figure
% subplot 211;plot(vett_f,mod1);grid;
% subplot 212;plot(vett_f,fas1);grid
% 
% figure
% subplot 211;plot(vett_f,mod2);grid;
% subplot 212;plot(vett_f,fas2);grid


%..............................................
% Fourier decomposition of a periodic signal
% period of signal to be assigned manually
% T=;
% 
% % values of signal at equally spaced discrete time values to be assigned
% % manually
% vect_F=[];
% 
% fftout=fft(vect_F);
% N=length(vect_F);
% df=1/T;
% fmax=(N/2-1)*df;
% vett_freq=0:df:fmax;
% modf(1)=1/N*abs(fftout(1));
% modf(2:N/2)=2/N*abs(fftout(2:N/2));
% fasf(1:N/2)=angle(fftout(1:N/2));
% 
% figure
% subplot 211;bar(vett_freq,modf);grid;
% subplot 212;bar(vett_freq,fasf);grid

%% Distributed load
% extraction of all the dof affected by the distributed load
freq=linspace(1,fmax,1000);
p0 = 1;
L=0.5;
p=p0;
Pvect=zeros(1,ndof);
idof5 = idf(5,:);
idof3 = idf(3,:);
Pvect(1,idof5) = [p , p*L/2 , -p*L^2 /12]; 
Pvect(1,idof3) = [p , p*L/2 , p*L^2 /12];

for k=1:length(freq)
    A=-omega(k)^2*MFF+1i*omega(k)*CFF+KFF;
    xf(:,k)=A\Pvect';
end

idofA = idf(4,2);
y4 = xf(idofA,:);

ampl = abs(y4);
phase = angle(y4);

